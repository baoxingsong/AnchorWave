#!/usr/bin/env python3
"""Summarize paired exact-candidate measurements for selector calibration.

The report keeps censored attempts separate, computes time/memory residual
quantiles, and evaluates the standalone predicted-fastest ranking only on
tasks for which all requested engines completed. It does not rewrite selector
constants automatically; biological scoring and safety policy remain explicit
code-review decisions.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import sys
from collections import Counter, defaultdict
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


def finite_float(value: object) -> Optional[float]:
    try:
        result = float(str(value))
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def quantile(values: Iterable[float], probability: float) -> Optional[float]:
    ordered = sorted(values)
    if not ordered:
        return None
    if len(ordered) == 1:
        return ordered[0]
    position = probability * (len(ordered) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def rounded(value: Optional[float]) -> Optional[float]:
    return None if value is None else round(value, 6)


def load_results(path: Path) -> List[Dict[str, str]]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"empty results: {path}")
    return rows


def completed_measurements(rows: Iterable[Mapping[str, str]]) -> List[Mapping[str, str]]:
    return [row for row in rows
            if row.get("status") == "completed" and
            row.get("reconstruction_ok") == "yes" and
            finite_float(row.get("internal_seconds")) is not None]


def ratio(numerator: object, denominator: object, scale: float = 1.0) -> Optional[float]:
    first = finite_float(numerator)
    second = finite_float(denominator)
    if first is None or second is None or second <= 0:
        return None
    return first / (second * scale)


def algorithm_summary(algorithm: str, rows: Sequence[Mapping[str, str]]) -> Dict[str, object]:
    completed = completed_measurements(rows)
    actual = [float(row["internal_seconds"]) for row in completed]
    p50_residual = [value for row in completed
                    if (value := ratio(row["internal_seconds"],
                                       row.get("predicted_minutes_p50"), 60.0)) is not None]
    p90_residual = [value for row in completed
                    if (value := ratio(row["internal_seconds"],
                                       row.get("predicted_minutes_p90"), 60.0)) is not None]
    p90_covered = [value <= 1.0 for value in p90_residual]
    memory_residual: List[float] = []
    for row in completed:
        actual_memory = (row.get("wfa_memory_peak_bytes")
                         if algorithm != "ksw2-full"
                         else str((finite_float(row.get("max_rss_kib")) or 0) * 1024.0))
        value = ratio(actual_memory, row.get("predicted_candidate_memory_bytes"))
        if value is not None:
            memory_residual.append(value)
    return {
        "attempts": len(rows),
        "status_counts": dict(sorted(Counter(row.get("status", "") for row in rows).items())),
        "completed_reconstructed": len(completed),
        "sum_internal_seconds": rounded(sum(actual)),
        "actual_seconds_p50": rounded(quantile(actual, 0.50)),
        "actual_seconds_p90": rounded(quantile(actual, 0.90)),
        "actual_over_predicted_p50_p50": rounded(quantile(p50_residual, 0.50)),
        "actual_over_predicted_p50_p90": rounded(quantile(p50_residual, 0.90)),
        "actual_over_predicted_p50_p99": rounded(quantile(p50_residual, 0.99)),
        "actual_over_predicted_p90_p90": rounded(quantile(p90_residual, 0.90)),
        "predicted_p90_empirical_coverage": rounded(
            sum(p90_covered) / len(p90_covered) if p90_covered else None),
        "actual_over_predicted_memory_p50": rounded(quantile(memory_residual, 0.50)),
        "actual_over_predicted_memory_p99": rounded(quantile(memory_residual, 0.99)),
        "actual_over_predicted_memory_max": rounded(max(memory_residual) if memory_residual else None),
    }


def paired_ranking(rows: Sequence[Mapping[str, str]], algorithms: Sequence[str]) -> Dict[str, object]:
    by_task: Dict[str, Dict[str, Mapping[str, str]]] = defaultdict(dict)
    for row in completed_measurements(rows):
        by_task[row["task_id"]][row["algorithm"]] = row
    complete_tasks = {
        task: values for task, values in by_task.items()
        if all(algorithm in values for algorithm in algorithms)
    }
    actual_fastest = Counter()
    predicted_fastest = Counter()
    correct = 0
    regrets: List[float] = []
    for values in complete_tasks.values():
        actual_choice = min(
            algorithms, key=lambda name: float(values[name]["internal_seconds"]))
        predicted_choice = min(
            algorithms,
            key=lambda name: finite_float(
                values[name].get("predicted_minutes_p50")) or math.inf)
        actual_fastest[actual_choice] += 1
        predicted_fastest[predicted_choice] += 1
        correct += actual_choice == predicted_choice
        fastest_time = float(values[actual_choice]["internal_seconds"])
        predicted_time = float(values[predicted_choice]["internal_seconds"])
        if fastest_time > 0:
            regrets.append(predicted_time / fastest_time)
    count = len(complete_tasks)
    return {
        "tasks_with_all_algorithms_completed": count,
        "actual_fastest_counts": dict(sorted(actual_fastest.items())),
        "predicted_fastest_counts": dict(sorted(predicted_fastest.items())),
        "predicted_fastest_accuracy": rounded(correct / count if count else None),
        "standalone_runtime_regret_p50": rounded(quantile(regrets, 0.50)),
        "standalone_runtime_regret_p90": rounded(quantile(regrets, 0.90)),
        "standalone_runtime_regret_max": rounded(max(regrets) if regrets else None),
        "limitation": (
            "This is isolated-engine ranking. Production makespan also includes "
            "memory wait, concurrent opportunity cost, and failure probability."
        ),
    }


def write_strata(path: Path, rows: Sequence[Mapping[str, str]]) -> None:
    groups: Dict[Tuple[str, str, str], List[Mapping[str, str]]] = defaultdict(list)
    for row in rows:
        groups[(row["algorithm"], row["length_bin"], row["ratio_bin"])].append(row)
    fields = (
        "algorithm", "length_bin", "ratio_bin", "attempts", "completed",
        "actual_seconds_p50", "actual_seconds_p90",
        "actual_over_predicted_p50_p50", "actual_over_predicted_p50_p90",
    )
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for (algorithm, length_bin, ratio_bin), values in sorted(groups.items()):
            completed = completed_measurements(values)
            actual = [float(row["internal_seconds"]) for row in completed]
            residual = [value for row in completed
                        if (value := ratio(row["internal_seconds"],
                                           row.get("predicted_minutes_p50"), 60.0)) is not None]
            writer.writerow({
                "algorithm": algorithm,
                "length_bin": length_bin,
                "ratio_bin": ratio_bin,
                "attempts": len(values),
                "completed": len(completed),
                "actual_seconds_p50": rounded(quantile(actual, 0.50)),
                "actual_seconds_p90": rounded(quantile(actual, 0.90)),
                "actual_over_predicted_p50_p50": rounded(quantile(residual, 0.50)),
                "actual_over_predicted_p50_p90": rounded(quantile(residual, 0.90)),
            })


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True,
                        help="JSON output; a sibling .strata.tsv is also written")
    args = parser.parse_args(argv)
    try:
        rows = load_results(args.results.resolve())
        algorithms = sorted({row["algorithm"] for row in rows})
        by_algorithm = {
            algorithm: algorithm_summary(
                algorithm, [row for row in rows if row["algorithm"] == algorithm])
            for algorithm in algorithms
        }
        report = {
            "format_version": 1,
            "source_results": str(args.results.resolve()),
            "algorithms": by_algorithm,
            "paired_ranking": paired_ranking(rows, algorithms),
            "censoring_rule": (
                "timeout, memory_limit, and analytic infeasibility are excluded "
                "from runtime residuals and retained as censored status counts"
            ),
        }
        output = args.output.resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                          encoding="utf-8")
        write_strata(output.with_suffix(".strata.tsv"), rows)
        print(json.dumps(report, sort_keys=True))
        return 0
    except (FileNotFoundError, OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
