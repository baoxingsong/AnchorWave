#!/usr/bin/env python3
"""Run a real three-pair proali -> trioGraphAli integration smoke test.

The fixture is generated deterministically and is intentionally small.  This
test exercises the real AnchorWave and minimap2 executables; it does not mock
pairwise MAF production.  It is kept out of the default CTest suite because a
system minimap2 installation is an optional dependency.
"""

from __future__ import annotations

import argparse
import random
import subprocess
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


Taxon = str


def mutate(sequence: str, seed: int, rate: float) -> str:
    generator = random.Random(seed)
    alphabet = "ACGT"
    result: List[str] = []
    for base in sequence:
        if generator.random() < rate:
            result.append(generator.choice(alphabet.replace(base, "")))
        else:
            result.append(base)
    return "".join(result)


def wrap_fasta(name: str, sequence: str, width: int = 80) -> str:
    lines = [f">{name}"]
    lines.extend(sequence[index : index + width] for index in range(0, len(sequence), width))
    return "\n".join(lines) + "\n"


def annotation(intervals: Iterable[Tuple[int, int]]) -> str:
    lines = ["##gff-version 3"]
    for index, (start, end) in enumerate(intervals, start=1):
        gene = f"gene{index}"
        transcript = f"tx{index}"
        lines.extend(
            [
                f"chr1\tsynthetic\tgene\t{start}\t{end}\t.\t+\t.\tID={gene};",
                f"chr1\tsynthetic\tmRNA\t{start}\t{end}\t.\t+\t.\tID={transcript};Parent={gene};",
                f"chr1\tsynthetic\tCDS\t{start}\t{end}\t.\t+\t0\tID=cds{index};Parent={transcript};",
            ]
        )
    return "\n".join(lines) + "\n"


def run(argv: Sequence[str], stdout: Optional[Path] = None) -> None:
    if stdout is None:
        subprocess.run(list(argv), check=True)
        return
    stdout.parent.mkdir(parents=True, exist_ok=True)
    with stdout.open("wb") as handle:
        subprocess.run(list(argv), check=True, stdout=handle)


def prepare_fixture(root: Path) -> Dict[Taxon, Path]:
    root.mkdir(parents=True, exist_ok=True)
    generator = random.Random(20260806)
    reference = "".join(generator.choice("ACGT") for _ in range(6000))
    sequences = {
        "I1": reference,
        "I2": mutate(reference, 20260807, 0.012),
        "O1": mutate(reference, 20260808, 0.035),
    }
    fasta_paths: Dict[Taxon, Path] = {}
    for taxon, sequence in sequences.items():
        path = root / f"{taxon}.fa"
        path.write_text(wrap_fasta("chr1", sequence), encoding="ascii")
        fasta_paths[taxon] = path
    intervals = [(301 + index * 900, 420 + index * 900) for index in range(6)]
    for reference_taxon in ("I1", "I2"):
        (root / f"{reference_taxon}.gff3").write_text(
            annotation(intervals), encoding="ascii"
        )
    return fasta_paths


def make_pair(
    anchorwave: str,
    root: Path,
    reference: Taxon,
    query: Taxon,
    fasta: Dict[Taxon, Path],
    cds: Dict[Taxon, Path],
    sam: Dict[Tuple[Taxon, Taxon], Path],
) -> Tuple[Path, Path]:
    pair_dir = root / f"{reference}__{query}"
    pair_dir.mkdir(parents=True, exist_ok=True)
    anchors = pair_dir / "anchors.tsv"
    maf = pair_dir / "alignment.maf"
    fragments = pair_dir / "fragments.maf"
    run(
        [
            anchorwave,
            "proali",
            "-i",
            str(root / f"{reference}.gff3"),
            "-as",
            str(cds[reference]),
            "-r",
            str(fasta[reference]),
            "-a",
            str(sam[(reference, query)]),
            "-ar",
            str(sam[(reference, reference)]),
            "-s",
            str(fasta[query]),
            "-n",
            str(anchors),
            "-R",
            "1",
            "-Q",
            "1",
            "-o",
            str(maf),
            "-f",
            str(fragments),
            "-t",
            "1",
            "-m",
            "20",
            "-e",
            "1",
            "-ns",
        ]
    )
    if not maf.is_file() or maf.stat().st_size == 0:
        raise RuntimeError(f"proali produced no MAF for {reference}/{query}")
    return maf, anchors


def execute(anchorwave: str, minimap2: str, root: Path) -> None:
    fasta = prepare_fixture(root)
    cds: Dict[Taxon, Path] = {}
    sam: Dict[Tuple[Taxon, Taxon], Path] = {}
    for reference in ("I1", "I2"):
        cds[reference] = root / f"{reference}.cds.fa"
        run(
            [
                anchorwave,
                "gff2seq",
                "-i",
                str(root / f"{reference}.gff3"),
                "-r",
                str(fasta[reference]),
                "-o",
                str(cds[reference]),
                "-m",
                "20",
            ]
        )
        targets = ("I1", "I2", "O1") if reference == "I1" else ("I2", "O1")
        for target in targets:
            sam[(reference, target)] = root / f"{reference}.cds.to.{target}.sam"
            run(
                [
                    minimap2,
                    "-x",
                    "splice",
                    "-t",
                    "1",
                    "-k",
                    "12",
                    "-a",
                    "-p",
                    "0.4",
                    "-N",
                    "20",
                    str(fasta[target]),
                    str(cds[reference]),
                ],
                stdout=sam[(reference, target)],
            )

    pairs = [
        ("I1", "I2"),
        ("I1", "O1"),
        ("I2", "O1"),
    ]
    pair_outputs = {
        pair: make_pair(anchorwave, root, pair[0], pair[1], fasta, cds, sam)
        for pair in pairs
    }
    taxa = root / "taxa.tsv"
    taxa.write_text(
        "taxon_id\trole\tfasta\tgff\tanchor_sam\tanchor_fasta\tcallability_bed\tquality_weight\n"
        + "\n".join(
            [
                f"I1\tingroup\t{fasta['I1']}\t{root / 'I1.gff3'}\t{sam[('I1', 'I1')]}\t{cds['I1']}\t.\t1",
                f"I2\tingroup\t{fasta['I2']}\t{root / 'I2.gff3'}\t{sam[('I2', 'I2')]}\t{cds['I2']}\t.\t1",
                f"O1\tprimary_outgroup\t{fasta['O1']}\t.\t.\t.\t.\t1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    pair_manifest = root / "pairwise.tsv"
    pair_manifest.write_text(
        "taxon_a\ttaxon_b\tmaf\tanchor_map\tscore_profile\tweight\n"
        + "\n".join(
            f"{left}\t{right}\t{pair_outputs[(left, right)][0]}\t"
            f"{pair_outputs[(left, right)][1]}\tproali-smoke\t1"
            for left, right in pairs
        )
        + "\n",
        encoding="utf-8",
    )
    tree = root / "species.nwk"
    tree.write_text("(O1,(I1,I2)A)ROOT;\n", encoding="ascii")
    output_prefix = root / "trio" / "run"
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    run(
        [
            anchorwave,
            "trioGraphAli",
            "--taxa",
            str(taxa),
            "--pairwise-manifest",
            str(pair_manifest),
            "--species-tree",
            str(tree),
            "--ancestor-node",
            "A",
            "--output-prefix",
            str(output_prefix),
            "--pairwise-scope",
            "triangles",
            "--validate-input-paths",
        ]
    )
    required = [
        ".graph.gfa",
        ".extant.maf",
        ".ancestor.blocks.fa",
        ".ancestor.children.maf",
        ".qc.tsv",
    ]
    for suffix in required:
        path = Path(str(output_prefix) + suffix)
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"missing integration output: {path}")
    print(f"real proali trio integration passed: {output_prefix}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--anchorwave", required=True)
    parser.add_argument("--minimap2", required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    execute(args.anchorwave, args.minimap2, args.work_dir.resolve())


if __name__ == "__main__":
    main()
