#!/usr/bin/env python3
"""Convert AnchorWave anchor rows into method-BED-compatible gap records.

AnchorWave anchor coordinates are one-based and inclusive.  The compact
B73/Mo17 trace workflow consumes zero-based half-open coordinates, so this
converter makes the coordinate convention explicit and validates every row.
It is intended for building a benchmark corpus from an anchor file before the
corresponding sequence-alignment stage has completed.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Optional, Sequence


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--anchors", type=Path, required=True)
    result.add_argument("--output", type=Path, required=True)
    result.add_argument("--reference-chromosome", default="1")
    result.add_argument("--query-chromosome", default="1")
    result.add_argument("--reference-start", type=int, default=0)
    result.add_argument("--reference-end", type=int, required=True)
    result.add_argument("--query-start", type=int, default=0)
    result.add_argument("--query-end", type=int, required=True)
    result.add_argument(
        "--include-anchor-sequences",
        action="store_true",
        help="also emit local/transcript anchor rows; default emits interanchor gaps only",
    )
    return result


def convert(args: argparse.Namespace) -> int:
    if args.reference_start < 0 or args.query_start < 0:
        raise ValueError("region starts must be non-negative")
    if args.reference_end <= args.reference_start:
        raise ValueError("reference end must be greater than start")
    if args.query_end <= args.query_start:
        raise ValueError("query end must be greater than start")
    source = args.anchors.resolve()
    output = args.output.resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    output.parent.mkdir(parents=True, exist_ok=True)

    emitted = 0
    excluded = 0
    with source.open("rt", encoding="utf-8") as reader, output.open(
        "xt", encoding="utf-8"
    ) as writer:
        writer.write(
            "# refChr\trefStart0\trefEnd0\tmethod\tscore\tstrand\t"
            "queryChr\tqueryStart0\tqueryEnd0\n"
        )
        for line_number, line in enumerate(reader, 1):
            if not line.strip() or line.startswith("#") or line.startswith("refChr\t"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                raise ValueError(f"malformed anchor row at line {line_number}")
            ref_chr, query_chr = fields[0], fields[3]
            if ref_chr != args.reference_chromosome or query_chr != args.query_chromosome:
                continue
            strand = fields[6]
            if strand not in ("+", "-"):
                raise ValueError(f"invalid strand at line {line_number}: {strand!r}")
            label = fields[7]
            if label != "interanchor" and not args.include_anchor_sequences:
                continue
            try:
                ref_begin_1 = int(fields[1])
                ref_end_1 = int(fields[2])
                query_begin_1 = int(fields[4])
                query_end_1 = int(fields[5])
            except ValueError as error:
                raise ValueError(f"invalid coordinate at line {line_number}") from error

            # Empty inter-anchor intervals are represented as begin=end+1.
            ref_start = ref_begin_1 - 1
            ref_end = max(ref_start, ref_end_1)
            query_start = query_begin_1 - 1
            query_end = max(query_start, query_end_1)
            if min(ref_start, ref_end, query_start, query_end) < 0:
                raise ValueError(f"negative converted coordinate at line {line_number}")
            if ref_end == ref_start or query_end == query_start:
                excluded += 1
                continue
            if not (
                args.reference_start <= ref_start
                and ref_end <= args.reference_end
                and args.query_start <= query_start
                and query_end <= args.query_end
            ):
                continue
            method = "ANCHOR_INTERANCHOR" if label == "interanchor" else "ANCHOR_SEQUENCE"
            writer.write(
                f"{ref_chr}\t{ref_start}\t{ref_end}\t{method}\t{fields[9]}\t"
                f"{strand}\t{query_chr}\t{query_start}\t{query_end}\n"
            )
            emitted += 1

    if emitted == 0:
        output.unlink()
        raise ValueError("no anchor intervals matched the requested regions")
    print(f"emitted={emitted}\texcluded_empty={excluded}\toutput={output}")
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    try:
        return convert(parser().parse_args(argv))
    except (FileExistsError, FileNotFoundError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
