#!/usr/bin/env python3
"""Generate a deterministic random Newick fixture.

The generated topology is deliberately labelled as a software stress-test.  It
must not be interpreted as a biological phylogeny.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import random
import re
import tempfile
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple


ALGORITHM_VERSION = "random-binary-newick-v1"
BIOLOGICAL_WARNING = (
    "RANDOM TREE: software parser/scheduler/stress-test fixture only; do not use "
    "it for biological ancestral sequence, copy-number, rearrangement, branch-event, "
    "or ancestral-karyotype inference."
)
_SAFE_LABEL = re.compile(r"^[A-Za-z_][A-Za-z0-9_.-]*$")
_SAFE_SAMPLE = re.compile(r"^[A-Za-z0-9_.-]+$")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def quote_newick_label(label: str) -> str:
    if not label:
        raise ValueError("Newick labels must not be empty")
    if _SAFE_LABEL.fullmatch(label):
        return label
    return "'" + label.replace("'", "''") + "'"


def validate_samples(samples: Iterable[str]) -> List[str]:
    values = list(samples)
    if len(values) < 2:
        raise ValueError("at least two distinct sample IDs are required")
    if len(values) != len(set(values)):
        raise ValueError("sample IDs must be unique")
    for sample in values:
        if not _SAFE_SAMPLE.fullmatch(sample):
            raise ValueError(
                f"unsafe sample ID {sample!r}; allowed characters are letters, digits, ._-"
            )
    return sorted(values)


def _random_length(generator: random.Random, minimum: float, maximum: float) -> str:
    return f"{generator.uniform(minimum, maximum):.6f}"


def generate_random_newick(
    samples: Iterable[str],
    seed: int,
    forced_sisters: Optional[Tuple[str, str]] = None,
    ancestor_name: str = "INGROUP_ANCESTOR",
    root_name: str = "ROOT",
    branch_length_min: float = 0.01,
    branch_length_max: float = 0.2,
) -> str:
    """Return a reproducible rooted, bifurcating random tree.

    ``forced_sisters`` is useful for trio tests: the two ingroups are first
    joined under the explicitly named ancestral node, then that clade is mixed
    with the outgroup tips while constructing the remaining random topology.
    """

    ordered = validate_samples(samples)
    if not (0.0 < branch_length_min <= branch_length_max):
        raise ValueError("branch lengths require 0 < minimum <= maximum")
    if not _SAFE_LABEL.fullmatch(ancestor_name):
        raise ValueError("ancestor_name must be a Newick-safe identifier")
    if not _SAFE_LABEL.fullmatch(root_name):
        raise ValueError("root_name must be a Newick-safe identifier")

    generator = random.Random(seed)
    nodes: List[str] = [quote_newick_label(sample) for sample in ordered]

    if forced_sisters is not None:
        left, right = forced_sisters
        if left == right or left not in ordered or right not in ordered:
            raise ValueError("forced_sisters must name two distinct selected samples")
        nodes.remove(quote_newick_label(left))
        nodes.remove(quote_newick_label(right))
        sister_node = (
            f"({quote_newick_label(left)}:{_random_length(generator, branch_length_min, branch_length_max)},"
            f"{quote_newick_label(right)}:{_random_length(generator, branch_length_min, branch_length_max)})"
            f"{quote_newick_label(ancestor_name)}"
        )
        nodes.append(sister_node)

    merge_count = 0
    while len(nodes) > 1:
        first_index, second_index = sorted(generator.sample(range(len(nodes)), 2))
        right = nodes.pop(second_index)
        left = nodes.pop(first_index)
        if generator.getrandbits(1):
            left, right = right, left
        parent = (
            f"({left}:{_random_length(generator, branch_length_min, branch_length_max)},"
            f"{right}:{_random_length(generator, branch_length_min, branch_length_max)})"
        )
        nodes.append(parent)
        merge_count += 1

    tree = nodes[0]
    # The last constructed node is the root and does not have a branch length.
    # With exactly two forced sisters, that named ancestral node is also root.
    root_label = "" if forced_sisters is not None and merge_count == 0 else quote_newick_label(root_name)
    return f"{tree}{root_label};\n"


def load_catalog_samples(
    manifest: Path, requested_samples: Sequence[str], cohorts: Sequence[str]
) -> List[str]:
    with manifest.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or "sample_id" not in rows[0] or "cohort" not in rows[0]:
        raise ValueError(f"{manifest}: expected sample_id and cohort columns")
    available = {row["sample_id"]: row for row in rows}
    missing = sorted(set(requested_samples) - set(available))
    if missing:
        raise ValueError(f"unknown sample ID(s): {', '.join(missing)}")
    selected = list(requested_samples)
    cohort_set = set(cohorts)
    selected.extend(row["sample_id"] for row in rows if row["cohort"] in cohort_set)
    # Stable first occurrence, then validate/sort in generate_random_newick.
    return list(dict.fromkeys(selected))


def build_parser() -> argparse.ArgumentParser:
    default_manifest = Path(__file__).resolve().parent / "data" / "maizegdb_assets.tsv"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=default_manifest)
    parser.add_argument("--sample", action="append", default=[], help="sample ID; repeatable")
    parser.add_argument("--cohort", action="append", default=[], help="catalog cohort; repeatable")
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--force-sisters", nargs=2, metavar=("INGROUP1", "INGROUP2"))
    parser.add_argument("--ancestor-name", default="INGROUP_ANCESTOR")
    parser.add_argument("--root-name", default="ROOT")
    parser.add_argument("--branch-length-min", type=float, default=0.01)
    parser.add_argument("--branch-length-max", type=float, default=0.2)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--metadata",
        type=Path,
        help="metadata JSON (default: OUTPUT.metadata.json)",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.sample and not args.cohort:
        raise SystemExit("select samples with at least one --sample or --cohort")
    samples = load_catalog_samples(args.manifest, args.sample, args.cohort)
    tree = generate_random_newick(
        samples,
        args.seed,
        tuple(args.force_sisters) if args.force_sisters else None,
        args.ancestor_name,
        args.root_name,
        args.branch_length_min,
        args.branch_length_max,
    )
    metadata_path = args.metadata or Path(str(args.output) + ".metadata.json")
    atomic_write_text(args.output, tree)
    metadata = {
        "schema_version": 1,
        "algorithm": ALGORITHM_VERSION,
        "seed": args.seed,
        "samples": sorted(samples),
        "forced_sisters": list(args.force_sisters) if args.force_sisters else None,
        "ancestor_name": args.ancestor_name,
        "root_name": args.root_name,
        "branch_length_range": [args.branch_length_min, args.branch_length_max],
        "manifest": str(args.manifest.resolve()),
        "manifest_sha256": sha256_file(args.manifest),
        "newick_sha256": hashlib.sha256(tree.encode("utf-8")).hexdigest(),
        "biological_warning": BIOLOGICAL_WARNING,
    }
    atomic_write_text(metadata_path, json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(f"wrote deterministic stress-test tree: {args.output}")
    print(BIOLOGICAL_WARNING)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
