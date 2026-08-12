#!/usr/bin/env python3
"""Sequential, paired exact-aligner benchmark for a real gap corpus.

Every candidate is launched in a fresh process and no two candidates overlap.
The runner is deliberately conservative: KSW2-full is skipped when its
analytic traceback allocation cannot fit the configured per-candidate budget,
and timeouts are recorded as right-censored observations rather than timings.
Results are resumable through immutable per-attempt JSON receipts.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


GIB = 1024 ** 3
MIB = 1024 ** 2
ALGORITHMS = (
    "ksw2-score-certified",
    "ksw2-full",
    "singletrack",
    "wfa-high",
    "wfa-medium",
    "wfa-low",
    "biwfa",
)
RESULT_FIELDS = (
    "task_id",
    "algorithm",
    "ref_length",
    "query_length",
    "max_length",
    "length_bin",
    "ratio_bin",
    "strand",
    "status",
    "return_code",
    "timed_out",
    "censored",
    "skip_reason",
    "score",
    "wall_seconds",
    "internal_seconds",
    "user_seconds",
    "system_seconds",
    "max_rss_kib",
    "wfa_memory_bytes",
    "wfa_memory_peak_bytes",
    "predicted_ksw_memory_bytes",
    "memory_budget_bytes",
    "reconstruction_ok",
    "alignment_hash_fnv1a64",
    "wfa_status",
    "configured_threads",
    "actual_threads",
    "wfa_parallel_support",
    "configured_biwfa_leaf_score",
    "biwfa_leaf_alignments",
    "biwfa_max_leaf_score",
    "estimated_score",
    "conservative_score",
    "profile_uncertainty",
    "sketch_jaccard",
    "length_ratio_predicted",
    "predicted_candidate_memory_bytes",
    "predicted_minutes_p50",
    "predicted_minutes_p90",
    "predicted_success_probability",
    "predicted_memory_feasible",
    "predicted_time_feasible",
    "certified_optimal_score",
    "certified_initial_band",
    "certified_maximum_band",
    "certified_final_band",
    "certified_traceback_attempts",
    "stdout_sha256",
    "stderr_sha256",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(MIB)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wt", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
    temporary.replace(path)


def parse_assignment(values: Sequence[str], allowed: Sequence[str], name: str) -> Dict[str, float]:
    result: Dict[str, float] = {}
    for item in values:
        if "=" not in item:
            raise ValueError(f"{name} must use ALGORITHM=VALUE: {item!r}")
        algorithm, raw_value = item.split("=", 1)
        if algorithm not in allowed:
            raise ValueError(f"unknown algorithm in {name}: {algorithm!r}")
        value = float(raw_value)
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"{name} value must be positive: {item!r}")
        result[algorithm] = value
    return result


def parse_algorithms(value: str) -> List[str]:
    algorithms: List[str] = []
    for algorithm in value.split(","):
        algorithm = algorithm.strip()
        if algorithm not in ALGORITHMS:
            raise ValueError(f"unknown algorithm: {algorithm!r}")
        if algorithm not in algorithms:
            algorithms.append(algorithm)
    if not algorithms:
        raise ValueError("at least one algorithm is required")
    return algorithms


def read_manifest(trace_dir: Path) -> Tuple[List[Dict[str, str]], List[str]]:
    path = trace_dir / "manifest.tsv"
    with path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or ())
        required = {
            "task_id", "ref_length", "query_length", "max_length",
            "length_bin", "ratio_bin", "strand", "ref_fasta", "query_fasta",
        }
        missing = required - set(fields)
        if missing:
            raise ValueError(f"manifest lacks columns: {sorted(missing)}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"empty manifest: {path}")
    return rows, fields


def stable_quantile_sample(rows: Sequence[Mapping[str, str]], sample_size: int) -> List[Mapping[str, str]]:
    """Balance strata and span their length range, including tail tasks."""
    if sample_size <= 0 or sample_size >= len(rows):
        return list(rows)
    strata: Dict[Tuple[str, str, str], List[Mapping[str, str]]] = defaultdict(list)
    for row in rows:
        strata[(row["length_bin"], row["ratio_bin"], row["strand"])].append(row)
    for values in strata.values():
        values.sort(key=lambda row: (int(row["max_length"]), row["task_id"]))

    ordered_strata = sorted(strata)
    quotas = {key: 0 for key in ordered_strata}
    remaining = sample_size
    while remaining:
        progressed = False
        for key in ordered_strata:
            if quotas[key] < len(strata[key]):
                quotas[key] += 1
                remaining -= 1
                progressed = True
                if remaining == 0:
                    break
        if not progressed:
            break

    selected: List[Mapping[str, str]] = []
    for key in ordered_strata:
        values = strata[key]
        count = quotas[key]
        if count == 0:
            continue
        if count == 1:
            indices = [len(values) - 1]
        else:
            indices = [round(index * (len(values) - 1) / (count - 1))
                       for index in range(count)]
        selected.extend(values[index] for index in indices)

    # Interleave strata so an interrupted run still contains broad coverage.
    selected.sort(key=lambda row: (
        hashlib.sha256((row["length_bin"] + "|" + row["ratio_bin"] + "|" +
                        row["strand"] + "|" + row["task_id"]).encode()).digest(),
        row["task_id"],
    ))
    return selected


def ksw2_full_memory_bytes(query_length: int, reference_length: int) -> int:
    # Keep this in sync with AlignmentAlgorithmSelector.cpp. KSW2 stores one
    # traceback byte per retained diagonal/column plus a fixed working set.
    diagonal_count = query_length + reference_length
    retained_columns = min(query_length, reference_length) + 160
    return 16 * MIB + diagonal_count * retained_columns


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
        seconds = 0.0
        for piece in elapsed.group(1).split(":"):
            seconds = seconds * 60.0 + float(piece)
        values["elapsed_wall_seconds"] = seconds
    return values


def parse_candidate_output(path: Path) -> Dict[str, str]:
    if not path.exists():
        return {}
    lines = [line for line in path.read_text(encoding="utf-8", errors="replace").splitlines()
             if line.strip()]
    if len(lines) != 2:
        return {}
    reader = csv.DictReader(lines, delimiter="\t")
    rows = list(reader)
    return rows[0] if len(rows) == 1 else {}


def blank_result(row: Mapping[str, str], algorithm: str, memory_bytes: int,
                 predicted_ksw: int) -> Dict[str, object]:
    result: Dict[str, object] = {field: "" for field in RESULT_FIELDS}
    for field in ("task_id", "ref_length", "query_length", "max_length",
                  "length_bin", "ratio_bin", "strand"):
        result[field] = row[field]
    result.update({
        "algorithm": algorithm,
        "memory_budget_bytes": memory_bytes,
        "predicted_ksw_memory_bytes": predicted_ksw,
        "timed_out": "no",
        "censored": "no",
    })
    return result


def terminate_process_group(process: subprocess.Popen, grace_seconds: float = 2.0) -> None:
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        process.wait(timeout=grace_seconds)
        return
    except subprocess.TimeoutExpired:
        pass
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    process.wait()


def run_attempt(row: Mapping[str, str], algorithm: str, args: argparse.Namespace,
                output: Path, memory_bytes: int, timeout_seconds: float) -> Dict[str, object]:
    predicted_ksw = ksw2_full_memory_bytes(
        int(row["query_length"]), int(row["ref_length"]))
    result = blank_result(row, algorithm, memory_bytes, predicted_ksw)
    receipt_dir = output / "receipts" / row["task_id"]
    receipt_dir.mkdir(parents=True, exist_ok=True)
    receipt_path = receipt_dir / f"{algorithm}.json"
    if receipt_path.exists():
        return json.loads(receipt_path.read_text(encoding="utf-8"))

    attempt_dir = output / "attempts" / row["task_id"]
    attempt_dir.mkdir(parents=True, exist_ok=True)
    stdout_path = attempt_dir / f"{algorithm}.stdout.tsv"
    stderr_path = attempt_dir / f"{algorithm}.stderr.txt"
    time_path = attempt_dir / f"{algorithm}.time.txt"
    if algorithm == "ksw2-full" and predicted_ksw > int(
            memory_bytes * args.ksw_admission_fraction):
        result.update({
            "status": "skipped_infeasible",
            "skip_reason": "analytic_ksw_memory_exceeds_budget",
            "censored": "yes",
        })
        atomic_json(receipt_path, result)
        return result

    command = [
        str(args.gnu_time), "-v", "-o", str(time_path),
        str(args.benchmark),
        "--reference", str((args.trace / row["ref_fasta"]).resolve()),
        "--query", str((args.trace / row["query_fasta"]).resolve()),
        "--algorithm", algorithm,
        "--memory-bytes", str(memory_bytes),
        "--wfa-threads", str(args.wfa_threads),
        "--min-offsets-per-thread", str(args.min_offsets_per_thread),
        "--biwfa-leaf-score", str(args.biwfa_leaf_score),
        "--memory-probe-score-interval", str(args.memory_probe_score_interval),
    ]
    started = time.monotonic()
    timed_out = False
    with stdout_path.open("wb") as stdout_handle, stderr_path.open("wb") as stderr_handle:
        process = subprocess.Popen(
            command, stdout=stdout_handle, stderr=stderr_handle,
            start_new_session=True,
        )
        try:
            process.wait(timeout=timeout_seconds)
        except subprocess.TimeoutExpired:
            timed_out = True
            terminate_process_group(process)
    wall_seconds = time.monotonic() - started
    time_values = parse_gnu_time(time_path)
    candidate = parse_candidate_output(stdout_path)
    result.update({
        "return_code": process.returncode,
        "wall_seconds": f"{wall_seconds:.6f}",
        "user_seconds": time_values["user_seconds"] if time_values["user_seconds"] is not None else "",
        "system_seconds": time_values["system_seconds"] if time_values["system_seconds"] is not None else "",
        "max_rss_kib": time_values["max_rss_kib"] if time_values["max_rss_kib"] is not None else "",
        "timed_out": "yes" if timed_out else "no",
        "censored": "yes" if timed_out else "no",
        "stdout_sha256": sha256_file(stdout_path),
        "stderr_sha256": sha256_file(stderr_path),
    })
    if timed_out:
        result["status"] = "timeout"
    elif not candidate:
        result["status"] = "process_error" if process.returncode else "malformed_output"
    else:
        result.update({
            "status": candidate.get("status", ""),
            "score": candidate.get("score", ""),
            "internal_seconds": candidate.get("seconds", ""),
            "wfa_memory_bytes": candidate.get("memory_bytes", ""),
            "wfa_memory_peak_bytes": candidate.get("memory_peak_bytes", ""),
            "reconstruction_ok": candidate.get("sequence_reconstruction", ""),
            "alignment_hash_fnv1a64": candidate.get(
                "alignment_hash_fnv1a64", ""),
            "wfa_status": candidate.get("wfa_status", ""),
            "configured_threads": candidate.get("configured_threads", ""),
            "actual_threads": candidate.get("actual_threads", ""),
            "wfa_parallel_support": candidate.get("wfa_parallel_support", ""),
            "configured_biwfa_leaf_score": candidate.get(
                "configured_biwfa_leaf_score", ""),
            "biwfa_leaf_alignments": candidate.get("biwfa_leaf_alignments", ""),
            "biwfa_max_leaf_score": candidate.get("biwfa_max_leaf_score", ""),
            "estimated_score": candidate.get("estimated_score", ""),
            "conservative_score": candidate.get("conservative_score", ""),
            "profile_uncertainty": candidate.get("profile_uncertainty", ""),
            "sketch_jaccard": candidate.get("sketch_jaccard", ""),
            "length_ratio_predicted": candidate.get("length_ratio", ""),
            "predicted_candidate_memory_bytes": candidate.get(
                "predicted_memory_bytes", ""),
            "predicted_minutes_p50": candidate.get("predicted_minutes_p50", ""),
            "predicted_minutes_p90": candidate.get("predicted_minutes_p90", ""),
            "predicted_success_probability": candidate.get(
                "predicted_success_probability", ""),
            "predicted_memory_feasible": candidate.get(
                "predicted_memory_feasible", ""),
            "predicted_time_feasible": candidate.get(
                "predicted_time_feasible", ""),
            "certified_optimal_score": candidate.get(
                "certified_optimal_score", ""),
            "certified_initial_band": candidate.get(
                "certified_initial_band", ""),
            "certified_maximum_band": candidate.get(
                "certified_maximum_band", ""),
            "certified_final_band": candidate.get(
                "certified_final_band", ""),
            "certified_traceback_attempts": candidate.get(
                "certified_traceback_attempts", ""),
        })
        if result["status"] == "memory_limit":
            result["censored"] = "yes"
        if process.returncode != 0:
            result["status"] = "process_error"
    atomic_json(receipt_path, result)
    return result


def write_results(output: Path, rows: Iterable[Mapping[str, object]]) -> None:
    path = output / "results.tsv"
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_FIELDS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def summarize(output: Path, results: Sequence[Mapping[str, object]], config: Mapping[str, object]) -> None:
    status_counts: Dict[str, int] = defaultdict(int)
    algorithm_counts: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    scores: Dict[str, List[Tuple[str, int]]] = defaultdict(list)
    paths: Dict[str, List[Tuple[str, str]]] = defaultdict(list)
    for row in results:
        status = str(row["status"])
        algorithm = str(row["algorithm"])
        status_counts[status] += 1
        algorithm_counts[algorithm][status] += 1
        if status == "completed" and row.get("score", "") != "":
            scores[str(row["task_id"])].append((algorithm, int(row["score"])))
            if row.get("alignment_hash_fnv1a64", ""):
                paths[str(row["task_id"])].append(
                    (algorithm, str(row["alignment_hash_fnv1a64"])))

    inconsistent: List[Dict[str, str]] = []
    for task_id, task_scores in sorted(scores.items()):
        unique = {score for _, score in task_scores}
        if len(unique) > 1:
            inconsistent.append({
                "task_id": task_id,
                "scores": ",".join(f"{algorithm}:{score}"
                                    for algorithm, score in sorted(task_scores)),
            })
    with (output / "score_inconsistencies.tsv").open(
            "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("task_id", "scores"), delimiter="\t")
        writer.writeheader()
        writer.writerows(inconsistent)

    path_differences: List[Dict[str, str]] = []
    for task_id, task_paths in sorted(paths.items()):
        if len({path_hash for _, path_hash in task_paths}) > 1:
            path_differences.append({
                "task_id": task_id,
                "alignment_hashes": ",".join(
                    f"{algorithm}:{path_hash}"
                    for algorithm, path_hash in sorted(task_paths)),
            })
    with (output / "alignment_path_differences.tsv").open(
            "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=("task_id", "alignment_hashes"), delimiter="\t")
        writer.writeheader()
        writer.writerows(path_differences)

    summary = {
        "format_version": 1,
        "updated_utc": utc_now(),
        "configuration": config,
        "attempts": len(results),
        "tasks_with_completed_candidate": len(scores),
        "status_counts": dict(sorted(status_counts.items())),
        "algorithm_status_counts": {
            algorithm: dict(sorted(counts.items()))
            for algorithm, counts in sorted(algorithm_counts.items())
        },
        "score_inconsistency_count": len(inconsistent),
        "alignment_path_difference_count": len(path_differences),
    }
    atomic_json(output / "summary.json", summary)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--trace", type=Path, required=True)
    result.add_argument("--benchmark", type=Path, required=True)
    result.add_argument("--output", type=Path, required=True)
    result.add_argument("--gnu-time", type=Path, default=Path("/opt/homebrew/bin/gtime"))
    result.add_argument("--algorithms", default=",".join(ALGORITHMS))
    result.add_argument("--sample-size", type=int, default=192,
                        help="0 uses every task; otherwise stratified quantile sample")
    result.add_argument("--task-id", action="append", default=[],
                        help="run an explicit manifest task; repeatable and overrides sampling")
    result.add_argument("--memory-gib", type=float, default=64.0)
    result.add_argument("--algorithm-memory-gib", action="append", default=[],
                        metavar="ALGORITHM=GIB")
    result.add_argument("--timeout-minutes", type=float, default=20.0)
    result.add_argument("--algorithm-timeout-minutes", action="append", default=[],
                        metavar="ALGORITHM=MINUTES")
    result.add_argument("--ksw-admission-fraction", type=float, default=0.90)
    result.add_argument("--wfa-threads", type=int, default=1)
    result.add_argument("--min-offsets-per-thread", type=int, default=500)
    result.add_argument("--biwfa-leaf-score", type=int, default=0)
    result.add_argument("--memory-probe-score-interval", type=int, default=3000)
    result.add_argument("--resume", action="store_true")
    return result


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parser().parse_args(argv)
    try:
        args.trace = args.trace.resolve()
        args.benchmark = args.benchmark.resolve()
        args.output = args.output.resolve()
        args.gnu_time = args.gnu_time.resolve()
        algorithms = parse_algorithms(args.algorithms)
        if args.sample_size < 0:
            raise ValueError("--sample-size cannot be negative")
        if args.memory_gib <= 0 or args.timeout_minutes <= 0:
            raise ValueError("memory and timeout must be positive")
        if not 0 < args.ksw_admission_fraction <= 1:
            raise ValueError("--ksw-admission-fraction must be in (0,1]")
        if args.wfa_threads <= 0 or args.min_offsets_per_thread <= 0:
            raise ValueError("WFA thread settings must be positive")
        if args.biwfa_leaf_score < 0 or args.memory_probe_score_interval <= 0:
            raise ValueError("invalid BiWFA leaf/probe setting")
        for path, label in ((args.benchmark, "benchmark"),
                            (args.gnu_time, "GNU time")):
            if not path.is_file() or not os.access(path, os.X_OK):
                raise FileNotFoundError(f"{label} is not executable: {path}")

        memory_overrides = parse_assignment(
            args.algorithm_memory_gib, ALGORITHMS, "--algorithm-memory-gib")
        timeout_overrides = parse_assignment(
            args.algorithm_timeout_minutes, ALGORITHMS,
            "--algorithm-timeout-minutes")
        rows, manifest_fields = read_manifest(args.trace)
        if args.task_id:
            requested = set(args.task_id)
            selected = [row for row in rows if row["task_id"] in requested]
            found = {row["task_id"] for row in selected}
            missing = requested - found
            if missing:
                raise ValueError(f"requested task IDs are absent: {sorted(missing)}")
            selected.sort(key=lambda row: row["task_id"])
        else:
            selected = stable_quantile_sample(rows, args.sample_size)
        config = {
            "format_version": 1,
            "trace": str(args.trace),
            "trace_manifest_sha256": sha256_file(args.trace / "manifest.tsv"),
            "benchmark": str(args.benchmark),
            "benchmark_sha256": sha256_file(args.benchmark),
            "gnu_time": str(args.gnu_time),
            "algorithms": algorithms,
            "sample_size_requested": args.sample_size,
            "explicit_task_ids": sorted(args.task_id),
            "selected_tasks": len(selected),
            "memory_gib": args.memory_gib,
            "algorithm_memory_gib": memory_overrides,
            "timeout_minutes": args.timeout_minutes,
            "algorithm_timeout_minutes": timeout_overrides,
            "ksw_admission_fraction": args.ksw_admission_fraction,
            "wfa_threads": args.wfa_threads,
            "min_offsets_per_thread": args.min_offsets_per_thread,
            "biwfa_leaf_score": args.biwfa_leaf_score,
            "memory_probe_score_interval": args.memory_probe_score_interval,
        }
        config_path = args.output / "configuration.json"
        if args.output.exists() and any(args.output.iterdir()):
            if not args.resume:
                raise FileExistsError(
                    f"output is not empty: {args.output}; use --resume or a new path")
            if not config_path.exists():
                raise ValueError("resume directory lacks configuration.json")
            previous = json.loads(config_path.read_text(encoding="utf-8"))
            if previous != config:
                raise ValueError("resume configuration differs from original run")
        else:
            args.output.mkdir(parents=True, exist_ok=True)
            atomic_json(config_path, config)
            with (args.output / "selected_manifest.tsv").open(
                    "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=manifest_fields, delimiter="\t")
                writer.writeheader()
                writer.writerows(selected)

        results: List[Dict[str, object]] = []
        total = len(selected) * len(algorithms)
        completed = 0
        for row in selected:
            for algorithm in algorithms:
                memory_gib = memory_overrides.get(algorithm, args.memory_gib)
                timeout_minutes = timeout_overrides.get(
                    algorithm, args.timeout_minutes)
                result_row = run_attempt(
                    row, algorithm, args, args.output,
                    int(memory_gib * GIB), timeout_minutes * 60.0)
                results.append(result_row)
                completed += 1
                print(
                    f"[{completed}/{total}] {row['task_id']} {algorithm} "
                    f"status={result_row['status']} wall={result_row['wall_seconds']}",
                    flush=True,
                )
                write_results(args.output, results)
                summarize(args.output, results, config)
        return 0
    except (FileExistsError, FileNotFoundError, OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
