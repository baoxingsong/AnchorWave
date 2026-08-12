#!/usr/bin/env python3
"""Stream-compare two pairwise MAF files without loading them into memory."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from itertools import zip_longest
from pathlib import Path
from typing import Iterator, List, Optional, Sequence, Tuple


@dataclass(frozen=True)
class MafRow:
    source: str
    start: int
    size: int
    strand: str
    source_size: int
    text: str

    @property
    def coordinates(self) -> Tuple[str, int, int, str, int]:
        return self.source, self.start, self.size, self.strand, self.source_size

    @property
    def ungapped(self) -> str:
        return self.text.replace("-", "")


@dataclass(frozen=True)
class MafBlock:
    attributes: Tuple[Tuple[str, str], ...]
    rows: Tuple[MafRow, ...]

    @property
    def score(self) -> Optional[str]:
        return dict(self.attributes).get("score")


@dataclass
class MafStats:
    blocks: int = 0
    rows: int = 0
    alignment_columns: int = 0
    ungapped_bases: int = 0
    gap_characters: int = 0
    invalid_blocks: int = 0


@dataclass(frozen=True)
class ScoringParameters:
    mismatch: int
    gap_open1: int
    gap_extend1: int
    gap_open2: int
    gap_extend2: int


def parse_attributes(line: str) -> Tuple[Tuple[str, str], ...]:
    attributes = []
    for field in line.split()[1:]:
        if "=" in field:
            key, value = field.split("=", 1)
            attributes.append((key, value))
    return tuple(attributes)


def parse_row(line: str, line_number: int, path: Path) -> MafRow:
    fields = line.split()
    if len(fields) != 7:
        raise ValueError(f"{path}:{line_number}: expected 7 fields in MAF s row")
    return MafRow(
        source=fields[1],
        start=int(fields[2]),
        size=int(fields[3]),
        strand=fields[4],
        source_size=int(fields[5]),
        text=fields[6],
    )


def read_maf(path: Path, stats: MafStats) -> Iterator[MafBlock]:
    attributes: Optional[Tuple[Tuple[str, str], ...]] = None
    rows: List[MafRow] = []

    def finish_block() -> Optional[MafBlock]:
        nonlocal attributes, rows
        if attributes is None and not rows:
            return None
        block = MafBlock(attributes or tuple(), tuple(rows))
        attributes = None
        rows = []
        stats.blocks += 1
        stats.rows += len(block.rows)
        if block.rows:
            stats.alignment_columns += len(block.rows[0].text)
        stats.ungapped_bases += sum(len(row.ungapped) for row in block.rows)
        stats.gap_characters += sum(row.text.count("-") for row in block.rows)
        if not validate_block(block):
            stats.invalid_blocks += 1
        return block

    with path.open("rt", encoding="ascii") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                block = finish_block()
                if block is not None:
                    yield block
                continue
            if line.startswith("#"):
                continue
            record_type = line.split(None, 1)[0]
            if record_type == "a":
                block = finish_block()
                if block is not None:
                    yield block
                attributes = parse_attributes(line)
                continue
            if record_type == "s":
                rows.append(parse_row(line, line_number, path))
                continue
            # Ignore optional MAF i/e/q records. AnchorWave pairwise MAF output
            # currently contains only a/s records, but accepting the optional
            # records keeps this comparator useful for converted files.
            if record_type in {"i", "e", "q"}:
                continue
            raise ValueError(f"{path}:{line_number}: unsupported MAF record: {line[:40]}")

    block = finish_block()
    if block is not None:
        yield block


def validate_block(block: MafBlock) -> bool:
    if not block.rows:
        return False
    column_count = len(block.rows[0].text)
    for row in block.rows:
        if row.strand not in {"+", "-"}:
            return False
        if row.start < 0 or row.size < 0 or row.source_size < 0:
            return False
        if row.start + row.size > row.source_size:
            return False
        if len(row.text) != column_count or len(row.ungapped) != row.size:
            return False
    return True


def score_pairwise_block(block: MafBlock,
                         scoring: ScoringParameters) -> Optional[int]:
    if len(block.rows) != 2 or not validate_block(block):
        return None
    first, second = (row.text for row in block.rows)
    score = 0
    column = 0
    while column < len(first):
        first_gap = first[column] == "-"
        second_gap = second[column] == "-"
        if first_gap and second_gap:
            return None
        if first_gap or second_gap:
            gap_in_first = first_gap
            end = column + 1
            while end < len(first):
                if ((first[end] == "-") != gap_in_first or
                        (second[end] == "-") == gap_in_first):
                    break
                end += 1
            length = end - column
            score -= min(
                scoring.gap_open1 + scoring.gap_extend1 * length,
                scoring.gap_open2 + scoring.gap_extend2 * length,
            )
            column = end
            continue
        if first[column] != second[column]:
            score -= scoring.mismatch
        column += 1
    return score


def compare_mafs(old_path: Path, new_path: Path,
                 scoring: Optional[ScoringParameters] = None) -> dict:
    old_stats = MafStats()
    new_stats = MafStats()
    comparison = {
        "paired_blocks": 0,
        "exact_blocks": 0,
        "coordinate_equivalent_blocks": 0,
        "score_equal_blocks": 0,
        "coordinate_mismatches": 0,
        "ungapped_sequence_mismatches": 0,
        "alignment_text_mismatches": 0,
        "score_mismatches": 0,
        "missing_old_blocks": 0,
        "missing_new_blocks": 0,
        "first_mismatch_block": None,
    }
    score_audit = None
    if scoring is not None:
        score_audit = {
            "parameters": asdict(scoring),
            "audited_paired_blocks": 0,
            "recomputed_score_equal_blocks": 0,
            "old_declared_score_inconsistencies": 0,
            "new_declared_score_inconsistencies": 0,
            "new_recomputed_score_improvements": 0,
            "old_recomputed_score_improvements": 0,
            "first_declared_score_inconsistency_block": None,
        }

    old_blocks = read_maf(old_path, old_stats)
    new_blocks = read_maf(new_path, new_stats)
    for block_number, pair in enumerate(
            zip_longest(old_blocks, new_blocks), start=1):
        old_block, new_block = pair
        if old_block is None:
            comparison["missing_old_blocks"] += 1
            continue
        if new_block is None:
            comparison["missing_new_blocks"] += 1
            continue

        comparison["paired_blocks"] += 1
        coordinates_equal = (
            tuple(row.coordinates for row in old_block.rows) ==
            tuple(row.coordinates for row in new_block.rows)
        )
        ungapped_equal = (
            tuple(row.ungapped for row in old_block.rows) ==
            tuple(row.ungapped for row in new_block.rows)
        )
        text_equal = (
            tuple(row.text for row in old_block.rows) ==
            tuple(row.text for row in new_block.rows)
        )
        score_equal = old_block.score == new_block.score

        if score_audit is not None:
            old_recomputed = score_pairwise_block(old_block, scoring)
            new_recomputed = score_pairwise_block(new_block, scoring)
            if old_recomputed is not None and new_recomputed is not None:
                score_audit["audited_paired_blocks"] += 1
                if old_recomputed == new_recomputed:
                    score_audit["recomputed_score_equal_blocks"] += 1
                elif new_recomputed > old_recomputed:
                    score_audit["new_recomputed_score_improvements"] += 1
                else:
                    score_audit["old_recomputed_score_improvements"] += 1
                old_inconsistent = old_block.score != str(old_recomputed)
                new_inconsistent = new_block.score != str(new_recomputed)
                score_audit["old_declared_score_inconsistencies"] += int(
                    old_inconsistent)
                score_audit["new_declared_score_inconsistencies"] += int(
                    new_inconsistent)
                if ((old_inconsistent or new_inconsistent) and
                        score_audit[
                            "first_declared_score_inconsistency_block"] is None):
                    score_audit[
                        "first_declared_score_inconsistency_block"] = block_number

        if coordinates_equal and ungapped_equal:
            comparison["coordinate_equivalent_blocks"] += 1
        if coordinates_equal and ungapped_equal and text_equal and score_equal:
            comparison["exact_blocks"] += 1
        if score_equal:
            comparison["score_equal_blocks"] += 1
        if not coordinates_equal:
            comparison["coordinate_mismatches"] += 1
        if not ungapped_equal:
            comparison["ungapped_sequence_mismatches"] += 1
        if not text_equal:
            comparison["alignment_text_mismatches"] += 1
        if not score_equal:
            comparison["score_mismatches"] += 1
        if ((not coordinates_equal or not ungapped_equal or
             not text_equal or not score_equal) and
                comparison["first_mismatch_block"] is None):
            comparison["first_mismatch_block"] = block_number

    invalid = old_stats.invalid_blocks or new_stats.invalid_blocks
    missing = (comparison["missing_old_blocks"] or
               comparison["missing_new_blocks"])
    biological_difference = (
        comparison["coordinate_mismatches"] or
        comparison["ungapped_sequence_mismatches"] or
        comparison["score_mismatches"]
    )
    if invalid:
        verdict = "INVALID_MAF"
    elif missing or biological_difference:
        verdict = "DIFFERENT"
    elif comparison["alignment_text_mismatches"]:
        verdict = "COORDINATE_EQUIVALENT"
    else:
        verdict = "EXACT"

    return {
        "old": {"path": str(old_path), **asdict(old_stats)},
        "new": {"path": str(new_path), **asdict(new_stats)},
        "comparison": comparison,
        "score_audit": score_audit,
        "verdict": verdict,
    }


def parse_score_penalties(value: str) -> ScoringParameters:
    try:
        values = [int(field) for field in value.split(",")]
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "score penalties must be integers") from error
    if len(values) != 5 or any(value < 0 for value in values):
        raise argparse.ArgumentTypeError(
            "expected five non-negative values: mismatch,O1,E1,O2,E2")
    return ScoringParameters(*values)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare old/new AnchorWave MAFs in streaming mode")
    parser.add_argument("old_maf", type=Path)
    parser.add_argument("new_maf", type=Path)
    parser.add_argument("--json", type=Path, dest="json_path")
    parser.add_argument(
        "--report-only", action="store_true",
        help="return success for any valid comparison, including DIFFERENT")
    parser.add_argument(
        "--score-penalties", type=parse_score_penalties,
        metavar="MISMATCH,O1,E1,O2,E2",
        help="recompute every pairwise block score with 2-piece affine costs")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        report = compare_mafs(args.old_maf, args.new_maf,
                              scoring=args.score_penalties)
    except (OSError, UnicodeError, ValueError) as error:
        print(json.dumps({"verdict": "ERROR", "error": str(error)}, indent=2))
        return 2

    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    if args.json_path:
        args.json_path.parent.mkdir(parents=True, exist_ok=True)
        args.json_path.write_text(rendered + "\n", encoding="utf-8")
    if args.report_only:
        return 0
    return 0 if report["verdict"] in {"EXACT", "COORDINATE_EQUIVALENT"} else 2


if __name__ == "__main__":
    raise SystemExit(main())
