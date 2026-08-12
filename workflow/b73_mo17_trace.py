#!/usr/bin/env python3
"""Prepare and replay compact B73/Mo17 inter-anchor alignment traces.

The trace source is an AnchorWave ``-b`` method BED.  ``prepare`` performs a
deterministic stratified sample, extracts the corresponding oriented sequence
pairs from B73 and Mo17 FASTA files, and writes a portable manifest.  ``replay``
runs the public ``anchorwave ali`` entry point without changing production
selector code.  ``summarize-pipeline`` records a bounded ``genoAli`` run, and
``compare`` compares two replay or bounded-pipeline result directories.

Important limitation: ``ali`` does not install the process-wide ``-M`` memory
scheduler and has no ``-t`` option.  Replay therefore uses independent ``ali``
processes and a conservative external concurrency limit.  Use a bounded
``genoAli`` run for authoritative ``-t/-M`` scheduler measurements.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import mmap
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Dict, List, Mapping, Optional, Sequence, Tuple


GIB = 1024 ** 3
MIB = 1024 ** 2
MANIFEST_FIELDS = (
    "task_id",
    "source_line",
    "ref_chr",
    "ref_start",
    "ref_end",
    "query_chr",
    "query_start",
    "query_end",
    "strand",
    "source_method",
    "source_score",
    "ref_length",
    "query_length",
    "max_length",
    "length_ratio",
    "length_bin",
    "ratio_bin",
    "ref_fasta",
    "query_fasta",
    "ref_sha256",
    "query_sha256",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def length_bin(length: int) -> str:
    if length <= 1_000:
        return "01_le_1kb"
    if length <= 10_000:
        return "02_1_10kb"
    if length <= 50_000:
        return "03_10_50kb"
    if length <= 250_000:
        return "04_50_250kb"
    return "05_gt_250kb"


def ratio_value(first: int, second: int) -> float:
    shorter = min(first, second)
    if shorter == 0:
        return math.inf if max(first, second) else 1.0
    return max(first, second) / shorter


def ratio_bin(value: float) -> str:
    if value <= 1.25:
        return "01_balanced"
    if value <= 4.0:
        return "02_imbalanced"
    return "03_extreme"


@dataclass(frozen=True)
class BedRecord:
    source_line: int
    ref_chr: str
    ref_start: int
    ref_end: int
    query_chr: str
    query_start: int
    query_end: int
    strand: str
    source_method: str
    source_score: str

    @property
    def ref_length(self) -> int:
        return self.ref_end - self.ref_start

    @property
    def query_length(self) -> int:
        return self.query_end - self.query_start

    @property
    def max_length(self) -> int:
        return max(self.ref_length, self.query_length)

    @property
    def length_ratio(self) -> float:
        return ratio_value(self.ref_length, self.query_length)

    @property
    def stratum(self) -> Tuple[str, str, str, str]:
        return (
            self.source_method,
            length_bin(self.max_length),
            ratio_bin(self.length_ratio),
            self.strand,
        )

    def stable_text(self) -> str:
        return "|".join(
            (
                str(self.source_line),
                self.ref_chr,
                str(self.ref_start),
                str(self.ref_end),
                self.query_chr,
                str(self.query_start),
                str(self.query_end),
                self.strand,
                self.source_method,
            )
        )


def read_method_bed(
    path: Path, include_filling: bool, maximum_sequence_bases: int
) -> Tuple[List[BedRecord], Counter]:
    records: List[BedRecord] = []
    excluded: Counter = Counter()
    with path.open("rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                excluded["malformed"] += 1
                continue
            try:
                ref_start = int(fields[1])
                ref_end = int(fields[2])
                query_start = int(fields[7])
                query_end = int(fields[8])
            except ValueError:
                excluded["invalid_coordinate"] += 1
                continue
            if min(ref_start, ref_end, query_start, query_end) < 0:
                excluded["negative_coordinate"] += 1
                continue
            if ref_end < ref_start or query_end < query_start:
                excluded["reversed_coordinate"] += 1
                continue
            method = fields[3]
            if method == "FILLING" and not include_filling:
                excluded["filling"] += 1
                continue
            record = BedRecord(
                source_line=line_number,
                ref_chr=fields[0],
                ref_start=ref_start,
                ref_end=ref_end,
                query_chr=fields[6],
                query_start=query_start,
                query_end=query_end,
                strand=fields[5],
                source_method=method,
                source_score=fields[4],
            )
            if record.strand not in ("+", "-"):
                excluded["invalid_strand"] += 1
                continue
            if maximum_sequence_bases > 0 and (
                record.ref_length > maximum_sequence_bases
                or record.query_length > maximum_sequence_bases
            ):
                excluded["over_maximum_sequence_bases"] += 1
                continue
            records.append(record)
    return records, excluded


def choose_records(records: Sequence[BedRecord], sample_size: int, seed: int) -> List[BedRecord]:
    if sample_size <= 0 or sample_size >= len(records):
        return sorted(records, key=lambda record: record.stable_text())

    strata: Dict[Tuple[str, str, str, str], List[BedRecord]] = defaultdict(list)
    for record in records:
        strata[record.stratum].append(record)
    seed_text = str(seed)
    for stratum in strata:
        strata[stratum].sort(
            key=lambda record: hashlib.sha256(
                (seed_text + "|" + record.stable_text()).encode("utf-8")
            ).digest()
        )

    # Balance the top-level source methods first, then cycle through the
    # length/ratio/strand strata within each method.  A flat lexicographic
    # stratum walk would omit later method names whenever sample_size is
    # smaller than the number of populated strata.
    method_strata: Dict[str, List[Tuple[str, str, str, str]]] = defaultdict(list)
    for stratum in strata:
        method_strata[stratum[0]].append(stratum)
    for method in method_strata:
        method_strata[method].sort(
            key=lambda stratum: hashlib.sha256(
                (seed_text + "|stratum|" + "|".join(stratum)).encode("utf-8")
            ).digest()
        )

    selected: List[BedRecord] = []
    offsets = {stratum: 0 for stratum in strata}
    stratum_offsets = {method: 0 for method in method_strata}
    ordered_methods = sorted(method_strata)
    while len(selected) < sample_size:
        made_progress = False
        for method in ordered_methods:
            candidate_strata = method_strata[method]
            for _ in range(len(candidate_strata)):
                stratum_index = stratum_offsets[method] % len(candidate_strata)
                stratum_offsets[method] += 1
                stratum = candidate_strata[stratum_index]
                offset = offsets[stratum]
                if offset >= len(strata[stratum]):
                    continue
                selected.append(strata[stratum][offset])
                offsets[stratum] += 1
                made_progress = True
                break
            if len(selected) == sample_size:
                break
        if not made_progress:
            break
    return sorted(selected, key=lambda record: record.stable_text())


@dataclass(frozen=True)
class FastaEntry:
    sequence_start: int
    sequence_end: int
    bases_per_line: int
    bytes_per_line: int
    sequence_length: int


class IndexedFasta:
    """Small mmap-based FASTA index that also handles chromosome-long lines."""

    def __init__(self, path: Path):
        self.path = path
        self.handle = path.open("rb")
        self.mapping = mmap.mmap(self.handle.fileno(), 0, access=mmap.ACCESS_READ)
        self.entries = self._index()

    def close(self) -> None:
        self.mapping.close()
        self.handle.close()

    def __enter__(self) -> "IndexedFasta":
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.close()

    def _index(self) -> Dict[str, FastaEntry]:
        entries: Dict[str, FastaEntry] = {}
        size = len(self.mapping)
        position = 0
        while position < size:
            if self.mapping[position : position + 1] != b">":
                marker = self.mapping.find(b"\n>", position)
                if marker < 0:
                    break
                position = marker + 1
            header_end = self.mapping.find(b"\n", position)
            if header_end < 0:
                raise ValueError(f"FASTA header lacks newline: {self.path}")
            header = self.mapping[position + 1 : header_end].decode(
                "utf-8", errors="replace"
            )
            name = header.split()[0]
            if not name:
                raise ValueError(f"empty FASTA record name: {self.path}")
            next_marker = self.mapping.find(b"\n>", header_end + 1)
            sequence_end = size if next_marker < 0 else next_marker + 1
            sequence_start = header_end + 1
            first_newline = self.mapping.find(b"\n", sequence_start, sequence_end)
            first_line_end = sequence_end if first_newline < 0 else first_newline
            first_line = self.mapping[sequence_start:first_line_end]
            bases_per_line = len(first_line.rstrip(b"\r"))
            if bases_per_line <= 0:
                raise ValueError(f"empty FASTA sequence: {name} in {self.path}")
            bytes_per_line = first_line_end - sequence_start
            if first_newline >= 0:
                bytes_per_line += 1

            content_end = sequence_end
            while content_end > sequence_start and self.mapping[content_end - 1] in (10, 13):
                content_end -= 1
            last_newline = self.mapping.rfind(b"\n", sequence_start, content_end)
            last_line_start = sequence_start if last_newline < 0 else last_newline + 1
            completed_lines = (last_line_start - sequence_start) // bytes_per_line
            last_bases = len(self.mapping[last_line_start:content_end].rstrip(b"\r"))
            sequence_length = completed_lines * bases_per_line + last_bases
            if name in entries:
                raise ValueError(f"duplicate FASTA record {name!r}: {self.path}")
            entries[name] = FastaEntry(
                sequence_start=sequence_start,
                sequence_end=sequence_end,
                bases_per_line=bases_per_line,
                bytes_per_line=bytes_per_line,
                sequence_length=sequence_length,
            )
            if next_marker < 0:
                break
            position = next_marker + 1
        if not entries:
            raise ValueError(f"no FASTA records found: {self.path}")
        return entries

    def extract(self, name: str, start: int, end: int) -> bytes:
        if name not in self.entries:
            raise KeyError(f"FASTA record {name!r} is absent from {self.path}")
        entry = self.entries[name]
        if start < 0 or end < start or end > entry.sequence_length:
            raise ValueError(
                f"invalid {name}:{start}-{end}; sequence length is {entry.sequence_length}"
            )
        if start == end:
            return b""

        def file_offset(coordinate: int) -> int:
            return (
                entry.sequence_start
                + (coordinate // entry.bases_per_line) * entry.bytes_per_line
                + coordinate % entry.bases_per_line
            )

        raw = self.mapping[file_offset(start) : file_offset(end)]
        sequence = raw.replace(b"\n", b"").replace(b"\r", b"")
        expected = end - start
        if len(sequence) != expected:
            raise ValueError(
                f"extracted {len(sequence)} rather than {expected} bases from "
                f"{name}:{start}-{end} in {self.path}"
            )
        return bytes(sequence).upper()


COMPLEMENT = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")


def reverse_complement(sequence: bytes) -> bytes:
    return sequence.translate(COMPLEMENT)[::-1]


def write_single_fasta(path: Path, name: str, sequence: bytes) -> None:
    with path.open("wb") as handle:
        handle.write(b">" + name.encode("utf-8") + b"\n")
        for offset in range(0, len(sequence), 80):
            handle.write(sequence[offset : offset + 80] + b"\n")


def task_identifier(record: BedRecord) -> str:
    return "gap_" + hashlib.sha256(record.stable_text().encode("utf-8")).hexdigest()[:16]


def prepare_trace(args: argparse.Namespace) -> int:
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(
            f"output directory is not empty: {output}; choose a new directory"
        )
    output.mkdir(parents=True, exist_ok=True)
    pair_dir = output / "pairs"
    pair_dir.mkdir()

    records, excluded = read_method_bed(
        args.methods.resolve(), args.include_filling, args.max_sequence_bases
    )
    selected = choose_records(records, args.sample_size, args.seed)
    if not selected:
        raise ValueError("no eligible method-BED records were selected")

    rows: List[Dict[str, str]] = []
    with IndexedFasta(args.reference.resolve()) as reference, IndexedFasta(
        args.query.resolve()
    ) as query:
        for record in selected:
            task_id = task_identifier(record)
            ref_sequence = reference.extract(
                record.ref_chr, record.ref_start, record.ref_end
            )
            query_sequence = query.extract(
                record.query_chr, record.query_start, record.query_end
            )
            if record.strand == "-":
                query_sequence = reverse_complement(query_sequence)
            ref_relative = Path("pairs") / f"{task_id}.ref.fa"
            query_relative = Path("pairs") / f"{task_id}.query.fa"
            write_single_fasta(
                output / ref_relative,
                f"{task_id}|B73|{record.ref_chr}:{record.ref_start}-{record.ref_end}",
                ref_sequence,
            )
            write_single_fasta(
                output / query_relative,
                f"{task_id}|Mo17|{record.query_chr}:{record.query_start}-{record.query_end}|{record.strand}",
                query_sequence,
            )
            ratio = record.length_ratio
            rows.append(
                {
                    "task_id": task_id,
                    "source_line": str(record.source_line),
                    "ref_chr": record.ref_chr,
                    "ref_start": str(record.ref_start),
                    "ref_end": str(record.ref_end),
                    "query_chr": record.query_chr,
                    "query_start": str(record.query_start),
                    "query_end": str(record.query_end),
                    "strand": record.strand,
                    "source_method": record.source_method,
                    "source_score": record.source_score,
                    "ref_length": str(record.ref_length),
                    "query_length": str(record.query_length),
                    "max_length": str(record.max_length),
                    "length_ratio": "inf" if math.isinf(ratio) else f"{ratio:.8f}",
                    "length_bin": length_bin(record.max_length),
                    "ratio_bin": ratio_bin(ratio),
                    "ref_fasta": str(ref_relative),
                    "query_fasta": str(query_relative),
                    "ref_sha256": sha256_bytes(ref_sequence),
                    "query_sha256": sha256_bytes(query_sequence),
                }
            )

    manifest = output / "manifest.tsv"
    with manifest.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    method_counts = Counter(row["source_method"] for row in rows)
    method_bases = Counter()
    stratum_counts = Counter()
    for row in rows:
        method_bases[row["source_method"]] += int(row["max_length"])
        stratum_counts["|".join((row["source_method"], row["length_bin"], row["ratio_bin"], row["strand"]))] += 1
    summary = {
        "format_version": 1,
        "kind": "b73_mo17_alignment_trace",
        "created_utc": utc_now(),
        "source_methods_bed": str(args.methods.resolve()),
        "source_reference_fasta": str(args.reference.resolve()),
        "source_query_fasta": str(args.query.resolve()),
        "seed": args.seed,
        "requested_sample_size": args.sample_size,
        "maximum_sequence_bases": args.max_sequence_bases,
        "eligible_records": len(records),
        "selected_records": len(rows),
        "excluded_records": dict(sorted(excluded.items())),
        "source_method_counts": dict(sorted(method_counts.items())),
        "source_method_max_bases": dict(sorted(method_bases.items())),
        "stratum_counts": dict(sorted(stratum_counts.items())),
        "manifest_sha256": sha256_file(manifest),
        "note": "source_method is a label from the source BED, not a forced replay algorithm",
    }
    write_json(output / "trace.json", summary)
    print(json.dumps(summary, sort_keys=True))
    return 0


def read_manifest(trace_dir: Path) -> List[Dict[str, str]]:
    manifest = trace_dir / "manifest.tsv"
    with manifest.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(MANIFEST_FIELDS) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"manifest lacks columns: {sorted(missing)}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"empty trace manifest: {manifest}")
    return rows


def parse_gnu_time(path: Path) -> Dict[str, Optional[float]]:
    values: Dict[str, Optional[float]] = {
        "user_seconds": None,
        "system_seconds": None,
        "max_rss_kib": None,
        "elapsed_wall_seconds": None,
    }
    if not path.exists():
        return values
    text = path.read_text(encoding="utf-8", errors="replace")
    patterns = {
        "user_seconds": r"User time \(seconds\):\s*([0-9.]+)",
        "system_seconds": r"System time \(seconds\):\s*([0-9.]+)",
        "max_rss_kib": r"Maximum resident set size \(kbytes\):\s*([0-9]+)",
    }
    for key, pattern in patterns.items():
        match = re.search(pattern, text)
        if match:
            values[key] = float(match.group(1))
    elapsed = re.search(
        r"Elapsed \(wall clock\) time.*?:\s*([0-9]+(?::[0-9]+)*(?:\.[0-9]+)?)",
        text,
    )
    if elapsed:
        pieces = [float(piece) for piece in elapsed.group(1).split(":")]
        seconds = 0.0
        for piece in pieces:
            seconds = seconds * 60.0 + piece
        values["elapsed_wall_seconds"] = seconds
    return values


def validate_ali_output(
    output_path: Path, expected_reference_sha: str, expected_query_sha: str
) -> bool:
    names: List[str] = []
    sequences: List[bytearray] = []
    with output_path.open("rb") as handle:
        for raw in handle:
            if raw.startswith(b">"):
                names.append(raw[1:].strip().decode("utf-8", errors="replace"))
                sequences.append(bytearray())
            elif sequences:
                sequences[-1].extend(raw.strip())
    if len(sequences) != 2:
        return False
    reference_ungapped = bytes(sequences[0]).replace(b"-", b"").upper()
    query_ungapped = bytes(sequences[1]).replace(b"-", b"").upper()
    return (
        sha256_bytes(reference_ungapped) == expected_reference_sha
        and sha256_bytes(query_ungapped) == expected_query_sha
    )


REPLAY_FIELDS = (
    "task_id",
    "source_method",
    "ref_length",
    "query_length",
    "status",
    "return_code",
    "wall_seconds",
    "user_seconds",
    "system_seconds",
    "max_rss_kib",
    "output_sha256",
    "reconstruction_ok",
)


def run_replay_task(
    row: Mapping[str, str],
    trace_dir: Path,
    output: Path,
    anchorwave: Path,
    gnu_time: Path,
    window_width: int,
    bt_minutes: float,
    timeout_seconds: float,
    keep_task_files: bool,
) -> Dict[str, str]:
    task_id = row["task_id"]
    task_dir = output / "task-files"
    stdout_path = task_dir / f"{task_id}.aligned.fa"
    stderr_path = task_dir / f"{task_id}.stderr.txt"
    time_path = task_dir / f"{task_id}.time.txt"
    command = [
        str(gnu_time),
        "-v",
        "-o",
        str(time_path),
        str(anchorwave),
        "ali",
        "-r",
        str(trace_dir / row["ref_fasta"]),
        "-s",
        str(trace_dir / row["query_fasta"]),
        "-w",
        str(window_width),
        "-bt",
        str(bt_minutes),
    ]
    started = time.monotonic()
    return_code = -1
    status = "error"
    with stdout_path.open("wb") as stdout_handle, stderr_path.open("wb") as stderr_handle:
        process = subprocess.Popen(
            command,
            stdout=stdout_handle,
            stderr=stderr_handle,
            start_new_session=True,
        )
        try:
            return_code = process.wait(timeout=timeout_seconds)
            status = "ok" if return_code == 0 else "failed"
        except subprocess.TimeoutExpired:
            status = "timeout"
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            return_code = 124
    wall_seconds = time.monotonic() - started
    timing = parse_gnu_time(time_path)
    output_sha = sha256_file(stdout_path) if stdout_path.exists() else ""
    reconstruction_ok = False
    if status == "ok":
        reconstruction_ok = validate_ali_output(
            stdout_path, row["ref_sha256"], row["query_sha256"]
        )
        if not reconstruction_ok:
            status = "invalid_output"
    result = {
        "task_id": task_id,
        "source_method": row["source_method"],
        "ref_length": row["ref_length"],
        "query_length": row["query_length"],
        "status": status,
        "return_code": str(return_code),
        "wall_seconds": f"{wall_seconds:.6f}",
        "user_seconds": "" if timing["user_seconds"] is None else str(timing["user_seconds"]),
        "system_seconds": "" if timing["system_seconds"] is None else str(timing["system_seconds"]),
        "max_rss_kib": "" if timing["max_rss_kib"] is None else str(int(timing["max_rss_kib"])),
        "output_sha256": output_sha,
        "reconstruction_ok": "yes" if reconstruction_ok else "no",
    }
    if not keep_task_files and status == "ok":
        stdout_path.unlink(missing_ok=True)
        stderr_path.unlink(missing_ok=True)
        time_path.unlink(missing_ok=True)
    return result


def replay_trace(args: argparse.Namespace) -> int:
    trace_dir = args.trace_dir.resolve()
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(
            f"output directory is not empty: {output}; comparisons require immutable runs"
        )
    output.mkdir(parents=True, exist_ok=True)
    (output / "task-files").mkdir()
    rows = read_manifest(trace_dir)
    if args.limit > 0:
        rows = rows[: args.limit]
    anchorwave = args.anchorwave.resolve()
    gnu_time = args.gnu_time.resolve()
    if not os.access(anchorwave, os.X_OK):
        raise ValueError(f"AnchorWave binary is not executable: {anchorwave}")
    if not os.access(gnu_time, os.X_OK):
        raise ValueError(f"GNU time is not executable: {gnu_time}")

    maximum_pair_bases = max(
        int(row["ref_length"]) + int(row["query_length"]) for row in rows
    )
    estimated_per_process = (
        args.window_width * args.window_width
        + 256 * MIB
        + 4 * maximum_pair_bases
    )
    effective_workers = args.threads
    reserve_bytes = 0
    if args.max_memory_gib > 0:
        maximum_bytes = int(args.max_memory_gib * GIB)
        reserve_bytes = max(GIB, int(maximum_bytes * 0.05))
        usable = max(0, maximum_bytes - reserve_bytes)
        memory_workers = usable // estimated_per_process
        if memory_workers < 1:
            raise ValueError(
                "the conservative ali replay envelope cannot admit one task; "
                "increase -M or reduce -w"
            )
        effective_workers = min(effective_workers, int(memory_workers))
    if effective_workers < 1:
        raise ValueError("threads must be positive")

    binary_sha = sha256_file(anchorwave)
    scenario = {
        "format_version": 1,
        "kind": "b73_mo17_ali_replay",
        "created_utc": utc_now(),
        "trace_dir": str(trace_dir),
        "trace_manifest_sha256": sha256_file(trace_dir / "manifest.tsv"),
        "anchorwave": str(anchorwave),
        "anchorwave_sha256": binary_sha,
        "requested_threads": args.threads,
        "requested_max_memory_gib": args.max_memory_gib,
        "effective_external_workers": effective_workers,
        "window_width": args.window_width,
        "bt_minutes": args.bt_minutes,
        "timeout_seconds_per_task": args.timeout_seconds,
        "requested_task_limit": args.limit,
        "estimated_external_per_process_bytes": estimated_per_process,
        "external_reserve_bytes": reserve_bytes,
        "limitation": (
            "anchorwave ali has neither -t nor -M and does not install the production "
            "process-wide scheduler; -t/-M here control conservative external process "
            "concurrency only. Use bounded genoAli for authoritative scheduler results."
        ),
    }
    write_json(output / "scenario.json", scenario)

    started = time.monotonic()
    results: List[Dict[str, str]] = []
    with ThreadPoolExecutor(max_workers=effective_workers) as executor:
        futures = [
            executor.submit(
                run_replay_task,
                row,
                trace_dir,
                output,
                anchorwave,
                gnu_time,
                args.window_width,
                args.bt_minutes,
                args.timeout_seconds,
                args.keep_task_files,
            )
            for row in rows
        ]
        for future in as_completed(futures):
            results.append(future.result())
    benchmark_wall = time.monotonic() - started
    results.sort(key=lambda result: result["task_id"])
    tasks_path = output / "tasks.tsv"
    with tasks_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=REPLAY_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)

    statuses = Counter(result["status"] for result in results)
    source_methods = Counter(result["source_method"] for result in results)
    source_method_wall = Counter()
    source_method_rss: Dict[str, float] = defaultdict(float)
    combined = hashlib.sha256()
    maximum_task_rss = 0
    for result in results:
        source_method_wall[result["source_method"]] += float(result["wall_seconds"])
        if result["max_rss_kib"]:
            rss = int(result["max_rss_kib"])
            maximum_task_rss = max(maximum_task_rss, rss)
            source_method_rss[result["source_method"]] = max(
                source_method_rss[result["source_method"]], rss
            )
        combined.update(
            (result["task_id"] + "\t" + result["output_sha256"] + "\n").encode(
                "utf-8"
            )
        )
    summary = {
        **scenario,
        "task_count": len(results),
        "statuses": dict(sorted(statuses.items())),
        "benchmark_wall_seconds": benchmark_wall,
        "sum_task_wall_seconds": sum(float(result["wall_seconds"]) for result in results),
        "maximum_single_ali_rss_kib": maximum_task_rss,
        "combined_output_sha256": combined.hexdigest(),
        "tasks_tsv_sha256": sha256_file(tasks_path),
        "source_method_counts": dict(sorted(source_methods.items())),
        "source_method_wall_seconds": dict(sorted(source_method_wall.items())),
        "source_method_max_single_ali_rss_kib": dict(sorted(source_method_rss.items())),
        "observed_replay_method_distribution": None,
        "rss_limitation": (
            "maximum_single_ali_rss_kib is not aggregate concurrent RSS; bounded "
            "genoAli GNU-time RSS is required for process-wide -M validation"
        ),
    }
    write_json(output / "summary.json", summary)
    print(json.dumps(summary, sort_keys=True))
    return 0 if statuses == Counter({"ok": len(results)}) else 2


def method_distribution(path: Path) -> Tuple[Dict[str, int], Dict[str, int]]:
    counts: Counter = Counter()
    bases: Counter = Counter()
    with path.open("rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            method = fields[3]
            counts[method] += 1
            try:
                bases[method] += max(int(fields[2]) - int(fields[1]), int(fields[8]) - int(fields[7]))
            except ValueError:
                pass
    return dict(sorted(counts.items())), dict(sorted(bases.items()))


def semantic_anchor_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        handle.readline()
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def summarize_pipeline(args: argparse.Namespace) -> int:
    result_dir = args.result_dir.resolve()
    prefix = result_dir / args.prefix
    required = {
        "anchors": prefix.with_suffix(".anchors"),
        "maf": Path(str(prefix) + ".maf"),
        "fragment_maf": Path(str(prefix) + ".fragment.maf"),
        "methods_bed": Path(str(prefix) + ".methods.bed"),
        "time": Path(str(prefix) + ".time.txt"),
        "scenario": result_dir / "scenario.json",
    }
    missing = [str(path) for path in required.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"bounded pipeline outputs are missing: {missing}")
    scenario = json.loads(required["scenario"].read_text(encoding="utf-8"))
    counts, bases = method_distribution(required["methods_bed"])
    timing = parse_gnu_time(required["time"])
    summary = {
        **scenario,
        "method_counts": counts,
        "method_max_bases": bases,
        "user_seconds": timing["user_seconds"],
        "system_seconds": timing["system_seconds"],
        "elapsed_wall_seconds": timing["elapsed_wall_seconds"],
        "maximum_process_rss_kib": timing["max_rss_kib"],
        "hashes": {
            "anchors_raw": sha256_file(required["anchors"]),
            "anchors_semantic_without_command_header": semantic_anchor_sha256(required["anchors"]),
            "maf": sha256_file(required["maf"]),
            "fragment_maf": sha256_file(required["fragment_maf"]),
            "methods_bed": sha256_file(required["methods_bed"]),
        },
    }
    write_json(result_dir / "summary.json", summary)
    print(json.dumps(summary, sort_keys=True))
    return 0


def numeric_delta(first: object, second: object) -> Optional[float]:
    if isinstance(first, (int, float)) and isinstance(second, (int, float)):
        return float(second) - float(first)
    return None


def compare_results(args: argparse.Namespace) -> int:
    left_dir = args.left.resolve()
    right_dir = args.right.resolve()
    left = json.loads((left_dir / "summary.json").read_text(encoding="utf-8"))
    right = json.loads((right_dir / "summary.json").read_text(encoding="utf-8"))
    comparison: Dict[str, object] = {
        "format_version": 1,
        "kind": "b73_mo17_benchmark_comparison",
        "created_utc": utc_now(),
        "left": str(left_dir),
        "right": str(right_dir),
        "left_kind": left.get("kind"),
        "right_kind": right.get("kind"),
    }
    if left.get("kind") == "b73_mo17_ali_replay" and right.get("kind") == "b73_mo17_ali_replay":
        def load_tasks(directory: Path) -> Dict[str, Dict[str, str]]:
            with (directory / "tasks.tsv").open("rt", encoding="utf-8", newline="") as handle:
                return {row["task_id"]: row for row in csv.DictReader(handle, delimiter="\t")}

        left_tasks = load_tasks(left_dir)
        right_tasks = load_tasks(right_dir)
        common = sorted(set(left_tasks) & set(right_tasks))
        hash_matches = sum(
            left_tasks[task]["output_sha256"] == right_tasks[task]["output_sha256"]
            for task in common
        )
        comparison.update(
            {
                "common_tasks": len(common),
                "left_only_tasks": len(set(left_tasks) - set(right_tasks)),
                "right_only_tasks": len(set(right_tasks) - set(left_tasks)),
                "byte_identical_common_outputs": hash_matches,
                "different_common_outputs": len(common) - hash_matches,
                "benchmark_wall_seconds": {
                    "left": left.get("benchmark_wall_seconds"),
                    "right": right.get("benchmark_wall_seconds"),
                    "delta_right_minus_left": numeric_delta(
                        left.get("benchmark_wall_seconds"), right.get("benchmark_wall_seconds")
                    ),
                },
                "maximum_single_ali_rss_kib": {
                    "left": left.get("maximum_single_ali_rss_kib"),
                    "right": right.get("maximum_single_ali_rss_kib"),
                    "delta_right_minus_left": numeric_delta(
                        left.get("maximum_single_ali_rss_kib"),
                        right.get("maximum_single_ali_rss_kib"),
                    ),
                },
                "source_method_counts": {
                    "left": left.get("source_method_counts"),
                    "right": right.get("source_method_counts"),
                },
            }
        )
    else:
        comparison.update(
            {
                "hashes_equal": left.get("hashes") == right.get("hashes"),
                "left_hashes": left.get("hashes"),
                "right_hashes": right.get("hashes"),
                "method_counts": {
                    "left": left.get("method_counts"),
                    "right": right.get("method_counts"),
                },
                "maximum_process_rss_kib": {
                    "left": left.get("maximum_process_rss_kib"),
                    "right": right.get("maximum_process_rss_kib"),
                    "delta_right_minus_left": numeric_delta(
                        left.get("maximum_process_rss_kib"),
                        right.get("maximum_process_rss_kib"),
                    ),
                },
                "elapsed_wall_seconds": {
                    "left": left.get("elapsed_wall_seconds"),
                    "right": right.get("elapsed_wall_seconds"),
                    "delta_right_minus_left": numeric_delta(
                        left.get("elapsed_wall_seconds"),
                        right.get("elapsed_wall_seconds"),
                    ),
                },
            }
        )
    output = args.output.resolve() if args.output else right_dir / "comparison-to-left.json"
    write_json(output, comparison)
    print(json.dumps(comparison, sort_keys=True))
    return 0


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(description=__doc__)
    subparsers = root.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser("prepare", help="extract a deterministic B73/Mo17 trace")
    prepare.add_argument("--methods", type=Path, required=True)
    prepare.add_argument("--reference", type=Path, required=True)
    prepare.add_argument("--query", type=Path, required=True)
    prepare.add_argument("--output", type=Path, required=True)
    prepare.add_argument("--sample-size", type=int, default=0, help="0 keeps all eligible records")
    prepare.add_argument("--seed", type=int, default=20260809)
    prepare.add_argument("--max-sequence-bases", type=int, default=500_000, help="0 disables this compactness bound")
    prepare.add_argument("--include-filling", action="store_true")
    prepare.set_defaults(function=prepare_trace)

    replay = subparsers.add_parser("replay", help="replay a trace through anchorwave ali")
    replay.add_argument("--trace-dir", type=Path, required=True)
    replay.add_argument("--output", type=Path, required=True)
    replay.add_argument("--anchorwave", type=Path, required=True)
    replay.add_argument("--gnu-time", type=Path, required=True)
    replay.add_argument("-t", "--threads", type=int, default=1)
    replay.add_argument("-M", "--max-memory-gib", type=float, default=0.0)
    replay.add_argument("-w", "--window-width", type=int, default=100_000)
    replay.add_argument("-bt", "--bt-minutes", type=float, default=30.0)
    replay.add_argument("--timeout-seconds", type=float, default=300.0)
    replay.add_argument("--limit", type=int, default=0, help="run only the first N manifest tasks; 0 runs all")
    replay.add_argument("--keep-task-files", action="store_true")
    replay.set_defaults(function=replay_trace)

    summarize = subparsers.add_parser("summarize-pipeline", help="summarize a bounded genoAli run")
    summarize.add_argument("--result-dir", type=Path, required=True)
    summarize.add_argument("--prefix", default="run")
    summarize.set_defaults(function=summarize_pipeline)

    compare = subparsers.add_parser("compare", help="compare two replay or bounded result directories")
    compare.add_argument("--left", type=Path, required=True)
    compare.add_argument("--right", type=Path, required=True)
    compare.add_argument("--output", type=Path)
    compare.set_defaults(function=compare_results)
    return root


def main(argv: Optional[Sequence[str]] = None) -> int:
    arguments = parser().parse_args(argv)
    try:
        return arguments.function(arguments)
    except (FileNotFoundError, FileExistsError, KeyError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
