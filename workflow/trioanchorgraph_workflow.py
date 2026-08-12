#!/usr/bin/env python3
"""Plan or execute the reproducible MaizeGDB TrioAnchorGraph input workflow.

This program never downloads data.  Use ``download_maizegdb.py --download``
explicitly first.  The execution boundary is also explicit: ``--plan`` writes
an auditable plan, ``--dry-run`` changes nothing, and only ``--execute`` runs
AnchorWave/minimap2.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from random_newick import (
    ALGORITHM_VERSION as TREE_ALGORITHM_VERSION,
    BIOLOGICAL_WARNING,
    generate_random_newick,
)


WORKFLOW_VERSION = "trioanchorgraph-maizegdb-v2"
CATALOG_COLUMNS = (
    "sample_id",
    "cohort",
    "taxon",
    "assembly_url",
    "assembly_filename",
    "assembly_md5",
    "assembly_size_bytes",
    "annotation_url",
    "annotation_filename",
    "annotation_md5",
    "annotation_size_bytes",
    "annotation_status",
)
SAFE_ID = re.compile(r"^[A-Za-z0-9_.-]+$")
MD5_RE = re.compile(r"^[0-9a-f]{32}$")
FORBIDDEN_FIELD_CHARACTERS = "\x00\t\r\n"


@dataclass(frozen=True)
class CatalogRecord:
    sample_id: str
    cohort: str
    taxon: str
    assembly_url: str
    assembly_filename: str
    assembly_md5: str
    assembly_size_bytes: int
    annotation_url: str
    annotation_filename: str
    annotation_md5: str
    annotation_size_bytes: Optional[int]
    annotation_status: str
    model: str
    accession: str
    notice_url: str


@dataclass
class Step:
    step_id: str
    kind: str
    inputs: List[str]
    outputs: List[str]
    argv: List[str] = field(default_factory=list)
    stdout_output: Optional[str] = None
    source_expected_size: Optional[int] = None
    source_expected_md5: Optional[str] = None
    description: str = ""

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class PairJob:
    reference: str
    query: str
    maf: Path
    anchors: Path
    fragments: Path


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()


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


def hash_file(path: Path, algorithm: str = "sha256") -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_safe_component(value: str, field_name: str) -> None:
    if value in (".", "..") or not SAFE_ID.fullmatch(value):
        raise ValueError(
            f"unsafe {field_name} {value!r}; expected a non-traversing letters/digits/._- component"
        )


def require_no_field_controls(value: str, field_name: str) -> None:
    if any(character in value for character in FORBIDDEN_FIELD_CHARACTERS):
        raise ValueError(
            f"{field_name} contains NUL, tab, carriage-return, or newline characters"
        )


def parse_optional_size(value: str, field_name: str, line: int) -> Optional[int]:
    if value == "":
        return None
    try:
        parsed = int(value)
    except ValueError as error:
        raise ValueError(f"catalog line {line}: {field_name} must be an integer") from error
    if parsed < 0:
        raise ValueError(f"catalog line {line}: {field_name} must be non-negative")
    return parsed


def load_catalog(path: Path) -> Dict[str, CatalogRecord]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing_columns = sorted(set(CATALOG_COLUMNS) - set(reader.fieldnames or []))
        if missing_columns:
            raise ValueError(f"{path}: missing column(s): {', '.join(missing_columns)}")
        records: Dict[str, CatalogRecord] = {}
        for line, row in enumerate(reader, 2):
            sample = row["sample_id"]
            require_safe_component(sample, f"catalog line {line} sample_id")
            require_safe_component(row["cohort"], f"catalog line {line} cohort")
            if sample in records:
                raise ValueError(f"catalog line {line}: duplicate sample_id {sample!r}")
            assembly_size = parse_optional_size(
                row["assembly_size_bytes"], "assembly_size_bytes", line
            )
            if assembly_size is None:
                raise ValueError(f"catalog line {line}: assembly_size_bytes is required")
            if Path(row["assembly_filename"]).name != row["assembly_filename"]:
                raise ValueError(f"catalog line {line}: assembly_filename must be a basename")
            require_safe_component(
                row["assembly_filename"], f"catalog line {line} assembly_filename"
            )
            annotation_size = parse_optional_size(
                row["annotation_size_bytes"], "annotation_size_bytes", line
            )
            if row["annotation_filename"] and (
                Path(row["annotation_filename"]).name != row["annotation_filename"]
            ):
                raise ValueError(f"catalog line {line}: annotation_filename must be a basename")
            if row["annotation_filename"]:
                require_safe_component(
                    row["annotation_filename"],
                    f"catalog line {line} annotation_filename",
                )
            records[sample] = CatalogRecord(
                sample_id=sample,
                cohort=row["cohort"],
                taxon=row["taxon"],
                assembly_url=row["assembly_url"],
                assembly_filename=row["assembly_filename"],
                assembly_md5=row["assembly_md5"].lower(),
                assembly_size_bytes=assembly_size,
                annotation_url=row["annotation_url"],
                annotation_filename=row["annotation_filename"],
                annotation_md5=row["annotation_md5"].lower(),
                annotation_size_bytes=annotation_size,
                annotation_status=row["annotation_status"],
                model=row.get("model", ""),
                accession=row.get("accession", ""),
                notice_url=row.get("notice_url", ""),
            )
    if not records:
        raise ValueError(f"{path}: empty catalog")
    return records


def require_selected(
    catalog: Mapping[str, CatalogRecord], ingroup1: str, ingroup2: str, outgroups: Sequence[str]
) -> List[CatalogRecord]:
    samples = [ingroup1, ingroup2] + list(outgroups)
    if len(outgroups) < 1:
        raise ValueError("at least one --outgroup is required")
    if len(samples) != len(set(samples)):
        raise ValueError("ingroup and outgroup sample IDs must be distinct")
    missing = [sample for sample in samples if sample not in catalog]
    if missing:
        raise ValueError(f"sample ID(s) absent from catalog: {', '.join(missing)}")
    for sample in (ingroup1, ingroup2):
        record = catalog[sample]
        if not record.annotation_filename or record.annotation_status == "assembly_only":
            raise ValueError(
                f"ingroup reference {sample} has no official annotation and cannot drive proali"
            )
    return [catalog[sample] for sample in samples]


def downloaded_path(data_root: Path, record: CatalogRecord, filename: str) -> Path:
    return (data_root / record.cohort / record.sample_id / filename).resolve()


def uncompressed_name(filename: str) -> str:
    return filename[:-3] if filename.endswith(".gz") else filename


class WorkflowPlan:
    def __init__(self, args: argparse.Namespace, catalog: Mapping[str, CatalogRecord]):
        self.args = args
        self.catalog = catalog
        self.records = require_selected(
            catalog, args.ingroup1, args.ingroup2, args.outgroup
        )
        self.selected = [record.sample_id for record in self.records]
        self.run_dir = args.run_dir.resolve()
        self.prepared_fasta: Dict[str, Path] = {}
        self.prepared_gff: Dict[str, Path] = {}
        self.cds: Dict[str, Path] = {}
        self.ref_sam: Dict[str, Path] = {}
        self.liftover_sam: Dict[Tuple[str, str], Path] = {}
        self.pairs: List[PairJob] = []
        self.steps: List[Step] = []
        self.tree_text: Optional[str] = None
        self.tree_path = self.run_dir / "manifests" / "species_tree.nwk"
        self.taxa_manifest = self.run_dir / "manifests" / "taxa.tsv"
        self.pairwise_manifest = self.run_dir / "manifests" / "pairwise.tsv"
        self.trio_command_file = self.run_dir / "manifests" / "run_trioGraphAli.sh"
        self.plan_path = self.run_dir / "workflow.plan.json"
        self._build()

    def _prepare_asset(
        self,
        record: CatalogRecord,
        asset_kind: str,
        filename: str,
        expected_size: int,
        expected_md5: str,
    ) -> Path:
        source = downloaded_path(self.args.data_root.resolve(), record, filename)
        destination = (
            self.run_dir
            / "prepared"
            / record.sample_id
            / asset_kind
            / uncompressed_name(filename)
        ).resolve()
        kind = "gunzip" if filename.endswith(".gz") else "copy"
        self.steps.append(
            Step(
                step_id=f"prepare.{record.sample_id}.{asset_kind}",
                kind=kind,
                inputs=[str(source)],
                outputs=[str(destination)],
                source_expected_size=expected_size,
                source_expected_md5=expected_md5,
                description=f"prepare {asset_kind} for {record.sample_id}",
            )
        )
        return destination

    def _add_minimap(self, reference: str, target: str, output: Path) -> None:
        argv = [
            self.args.minimap2,
            "-x",
            "splice",
            "-t",
            str(self.args.minimap2_threads),
            "-k",
            "12",
            "-a",
            "-p",
            "0.4",
            "-N",
            "20",
            str(self.prepared_fasta[target]),
            str(self.cds[reference]),
        ]
        self.steps.append(
            Step(
                step_id=f"minimap2.{reference}.to.{target}",
                kind="subprocess",
                inputs=[str(self.prepared_fasta[target]), str(self.cds[reference])],
                outputs=[str(output)],
                argv=argv,
                stdout_output=str(output),
                description=f"map {reference} anchor CDS to {target}",
            )
        )

    def _add_pair(self, reference: str, query: str) -> None:
        pair_dir = (self.run_dir / "pairwise" / f"{reference}__{query}").resolve()
        anchors = pair_dir / "anchors.tsv"
        maf = pair_dir / "alignment.maf"
        fragments = pair_dir / "fragments.maf"
        query_sam = self.liftover_sam[(reference, query)]
        argv = [
            self.args.anchorwave,
            "proali",
            "-i",
            str(self.prepared_gff[reference]),
            "-as",
            str(self.cds[reference]),
            "-r",
            str(self.prepared_fasta[reference]),
            "-a",
            str(query_sam),
            "-ar",
            str(self.ref_sam[reference]),
            "-s",
            str(self.prepared_fasta[query]),
            "-n",
            str(anchors),
            "-R",
            str(self.args.reference_depth),
            "-Q",
            str(self.args.query_depth),
            "-o",
            str(maf),
            "-f",
            str(fragments),
            "-t",
            str(self.args.anchorwave_threads),
            "-m",
            str(self.args.min_exon),
            "-e",
            str(self.args.expected_gene_copies),
        ]
        if self.args.window_width is not None:
            argv.extend(["-w", str(self.args.window_width)])
        if self.args.fa3 is not None:
            argv.extend(["-fa3", str(self.args.fa3)])
        if self.args.no_new_anchors:
            argv.append("-ns")
        if self.args.exon_model:
            argv.append("-x")
        self.steps.append(
            Step(
                step_id=f"proali.{reference}.to.{query}",
                kind="subprocess",
                inputs=[
                    str(self.prepared_gff[reference]),
                    str(self.cds[reference]),
                    str(self.prepared_fasta[reference]),
                    str(query_sam),
                    str(self.ref_sam[reference]),
                    str(self.prepared_fasta[query]),
                ],
                outputs=[str(anchors), str(maf), str(fragments)],
                argv=argv,
                description=f"AnchorWave proali {reference} (reference) vs {query}",
            )
        )
        self.pairs.append(PairJob(reference, query, maf, anchors, fragments))

    def _build(self) -> None:
        ingroups = [self.args.ingroup1, self.args.ingroup2]
        for record in self.records:
            self.prepared_fasta[record.sample_id] = self._prepare_asset(
                record,
                "assembly",
                record.assembly_filename,
                record.assembly_size_bytes,
                record.assembly_md5,
            )
        for reference in ingroups:
            record = self.catalog[reference]
            assert record.annotation_size_bytes is not None
            self.prepared_gff[reference] = self._prepare_asset(
                record,
                "annotation",
                record.annotation_filename,
                record.annotation_size_bytes,
                record.annotation_md5,
            )

        for reference in ingroups:
            reference_dir = (self.run_dir / "references" / reference).resolve()
            cds = reference_dir / "anchors.fa"
            self.cds[reference] = cds
            argv = [
                self.args.anchorwave,
                "gff2seq",
                "-i",
                str(self.prepared_gff[reference]),
                "-r",
                str(self.prepared_fasta[reference]),
                "-o",
                str(cds),
                "-m",
                str(self.args.min_exon),
            ]
            if self.args.exon_model:
                argv.append("-x")
            self.steps.append(
                Step(
                    step_id=f"gff2seq.{reference}",
                    kind="subprocess",
                    inputs=[str(self.prepared_gff[reference]), str(self.prepared_fasta[reference])],
                    outputs=[str(cds)],
                    argv=argv,
                    description=f"extract anchor sequences from {reference}",
                )
            )

            targets = [reference]
            if reference == self.args.ingroup1:
                targets.append(self.args.ingroup2)
            targets.extend(self.args.outgroup)
            for target in targets:
                output = reference_dir / f"anchors.to.{target}.sam"
                self.liftover_sam[(reference, target)] = output
                self._add_minimap(reference, target, output)
            self.ref_sam[reference] = self.liftover_sam[(reference, reference)]

        # Triangle scope: one ingroup edge and two edges per outgroup.
        self._add_pair(self.args.ingroup1, self.args.ingroup2)
        for outgroup in self.args.outgroup:
            self._add_pair(self.args.ingroup1, outgroup)
            self._add_pair(self.args.ingroup2, outgroup)

        if self.args.random_tree_seed is not None:
            self.tree_text = generate_random_newick(
                self.selected,
                self.args.random_tree_seed,
                forced_sisters=(self.args.ingroup1, self.args.ingroup2),
                ancestor_name=self.args.ancestor_node,
            )
        elif self.args.species_tree is None:
            raise ValueError("select --species-tree or explicitly set --random-tree-seed")

    def taxon_manifest_text(self) -> str:
        lines = [
            "taxon_id\trole\tfasta\tgff\tanchor_sam\tanchor_fasta\tcallability_bed\tquality_weight"
        ]
        for index, sample in enumerate(self.selected):
            if sample == self.args.ingroup1:
                role = "ingroup_reference"
            elif sample == self.args.ingroup2:
                role = "ingroup"
            elif index == 2:
                role = "primary_outgroup"
            else:
                role = "outgroup"
            gff = str(self.prepared_gff[sample]) if sample in self.prepared_gff else "."
            anchor_sam = str(self.ref_sam[sample]) if sample in self.ref_sam else "."
            anchor_fasta = str(self.cds[sample]) if sample in self.cds else "."
            lines.append(
                "\t".join(
                    [
                        sample,
                        role,
                        str(self.prepared_fasta[sample]),
                        gff,
                        anchor_sam,
                        anchor_fasta,
                        ".",
                        "1",
                    ]
                )
            )
        return "\n".join(lines) + "\n"

    def pairwise_manifest_text(self) -> str:
        profile = (
            f"proali:R{self.args.reference_depth}:Q{self.args.query_depth}:"
            f"m{self.args.min_exon}:e{self.args.expected_gene_copies}"
        )
        lines = ["taxon_a\ttaxon_b\tmaf\tanchor_map\tscore_profile\tweight"]
        for pair in self.pairs:
            lines.append(
                "\t".join(
                    [
                        pair.reference,
                        pair.query,
                        str(pair.maf),
                        str(pair.anchors),
                        profile,
                        "1",
                    ]
                )
            )
        return "\n".join(lines) + "\n"

    def trio_argv(self) -> List[str]:
        argv = [
            self.args.anchorwave,
            "trioGraphAli",
            "--taxa",
            str(self.taxa_manifest),
            "--pairwise-manifest",
            str(self.pairwise_manifest),
            "--species-tree",
            str(self.tree_path),
            "--ancestor-node",
            self.args.ancestor_node,
            "--output-prefix",
            str((self.run_dir / "trio" / "trioanchorgraph").resolve()),
            "--pairwise-scope",
            "triangles",
            "--validate-input-paths",
        ]
        if self.args.copy_relations is not None:
            argv.extend(["--copy-relations", str(self.args.copy_relations.resolve())])
        if self.args.copy_mode != "constrained":
            argv.extend(["--copy-mode", self.args.copy_mode])
        for specification in self.args.copy_number:
            argv.extend(["--copy-number", specification])
        if self.args.block_projections is not None:
            argv.extend(["--block-projections", str(self.args.block_projections.resolve())])
        # FASTA evidence validation is intentionally enabled through the
        # trioGraphAli default. Never emit --skip-fasta-validation here.
        return argv

    def to_dict(self) -> Dict[str, object]:
        source_assets = []
        for record in self.records:
            source_assets.append(
                {
                    "sample_id": record.sample_id,
                    "cohort": record.cohort,
                    "taxon": record.taxon,
                    "assembly": {
                        "path": str(
                            downloaded_path(
                                self.args.data_root.resolve(), record, record.assembly_filename
                            )
                        ),
                        "url": record.assembly_url,
                        "size": record.assembly_size_bytes,
                        "md5": record.assembly_md5,
                    },
                    "annotation": (
                        {
                            "path": str(
                                downloaded_path(
                                    self.args.data_root.resolve(),
                                    record,
                                    record.annotation_filename,
                                )
                            ),
                            "url": record.annotation_url,
                            "size": record.annotation_size_bytes,
                            "md5": record.annotation_md5,
                            "status": record.annotation_status,
                        }
                        if record.annotation_filename
                        else None
                    ),
                    "accession": record.accession,
                    "notice_url": record.notice_url,
                }
            )
        return {
            "schema_version": 1,
            "workflow_version": WORKFLOW_VERSION,
            "catalog": str(self.args.manifest.resolve()),
            "catalog_sha256": hash_file(self.args.manifest.resolve()),
            "run_dir": str(self.run_dir),
            "selection": {
                "ingroup_reference": self.args.ingroup1,
                "ingroup": self.args.ingroup2,
                "primary_outgroup": self.args.outgroup[0],
                "additional_outgroups": self.args.outgroup[1:],
            },
            "random_tree": (
                {
                    "algorithm": TREE_ALGORITHM_VERSION,
                    "seed": self.args.random_tree_seed,
                    "newick": self.tree_text,
                    "warning": BIOLOGICAL_WARNING,
                }
                if self.tree_text is not None
                else None
            ),
            "source_assets": source_assets,
            "steps": [step.to_dict() for step in self.steps],
            "taxa_manifest": str(self.taxa_manifest),
            "pairwise_manifest": str(self.pairwise_manifest),
            "pair_count": len(self.pairs),
            "execution_policy": {
                "step_concurrency": 1,
                "pairwise_concurrency": 1,
                "minimap2_threads": self.args.minimap2_threads,
                "anchorwave_threads": self.args.anchorwave_threads,
                "source_size_validation": not self.args.skip_source_integrity,
                "source_md5_validation": not (
                    self.args.skip_source_integrity or self.args.skip_source_md5
                ),
            },
            "trioGraphAli_validation": {
                "fasta_evidence": "enabled_by_cli_default",
                "manifest_input_paths": True,
            },
            "trioGraphAli_argv": self.trio_argv(),
        }


class Executor:
    def __init__(
        self,
        plan: WorkflowPlan,
        resume: bool,
        verify_size: bool,
        verify_md5: bool,
        resolved_tools: Mapping[str, str],
    ):
        self.plan = plan
        self.resume = resume
        self.verify_size = verify_size
        self.verify_md5 = verify_md5
        self.resolved_tools = dict(resolved_tools)
        self.tool_provenance: Dict[str, Dict[str, object]] = {}
        for requested, resolved in self.resolved_tools.items():
            tool_path = Path(resolved)
            self.tool_provenance[requested] = {
                "requested": requested,
                "resolved": resolved,
                "size": tool_path.stat().st_size,
                "mtime_ns": tool_path.stat().st_mtime_ns,
                "sha256": hash_file(tool_path),
            }
        self.receipt_dir = plan.run_dir / "provenance" / "receipts"
        self.log_dir = plan.run_dir / "provenance" / "logs"
        self.event_path = plan.run_dir / "provenance" / "events.jsonl"

    def event(self, payload: Mapping[str, object]) -> None:
        self.event_path.parent.mkdir(parents=True, exist_ok=True)
        record = dict(payload)
        record["time_utc"] = utc_now()
        with self.event_path.open("a", encoding="utf-8", newline="\n") as handle:
            handle.write(json.dumps(record, sort_keys=True) + "\n")
            handle.flush()
            os.fsync(handle.fileno())

    @staticmethod
    def fingerprint(path: Path) -> Dict[str, object]:
        status = path.stat()
        return {
            "path": str(path),
            "size": status.st_size,
            "mtime_ns": status.st_mtime_ns,
        }

    def signature(self, step: Step) -> str:
        inputs = []
        for input_name in step.inputs:
            path = Path(input_name)
            if not path.is_file():
                raise FileNotFoundError(f"{step.step_id}: missing input {path}")
            inputs.append(self.fingerprint(path))
        payload = {
            "workflow_version": WORKFLOW_VERSION,
            "step": step.to_dict(),
            "inputs": inputs,
            "tool": self.tool_provenance.get(step.argv[0]) if step.argv else None,
        }
        return hashlib.sha256(
            json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest()

    def receipt_path(self, step: Step) -> Path:
        return self.receipt_dir / f"{step.step_id}.json"

    def reusable(self, step: Step, signature: str) -> bool:
        if not self.resume:
            return False
        receipt_path = self.receipt_path(step)
        if not receipt_path.is_file():
            return False
        try:
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            if receipt.get("signature") != signature:
                return False
            observed = [self.fingerprint(Path(path)) for path in step.outputs]
            return observed == receipt.get("outputs")
        except (OSError, ValueError, TypeError):
            return False

    def verify_source(self, step: Step) -> None:
        source = Path(step.inputs[0])
        if not source.is_file():
            raise FileNotFoundError(
                f"{source} is missing; this workflow never downloads data. "
                "Run workflow/download_maizegdb.py with explicit --download first."
            )
        if not self.verify_size:
            return
        observed_size = source.stat().st_size
        if step.source_expected_size is not None and observed_size != step.source_expected_size:
            raise ValueError(
                f"{source}: expected {step.source_expected_size} bytes, observed {observed_size}"
            )
        expected_md5 = step.source_expected_md5 or ""
        if self.verify_md5 and MD5_RE.fullmatch(expected_md5):
            observed_md5 = hash_file(source, "md5")
            if observed_md5 != expected_md5:
                raise ValueError(
                    f"{source}: expected MD5 {expected_md5}, observed {observed_md5}"
                )
        elif self.verify_md5 and expected_md5 and expected_md5 != "missing_in_source_manifest":
            raise ValueError(f"{source}: invalid catalog MD5 {expected_md5!r}")

    def execute_prepare(self, step: Step) -> None:
        self.verify_source(step)
        source = Path(step.inputs[0])
        destination = Path(step.outputs[0])
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = destination.parent / f".{destination.name}.{os.getpid()}.tmp"
        if temporary.exists():
            temporary.unlink()
        opener = gzip.open if step.kind == "gunzip" else open
        try:
            with opener(source, "rb") as input_handle, temporary.open("wb") as output_handle:
                shutil.copyfileobj(input_handle, output_handle, length=4 * 1024 * 1024)
                output_handle.flush()
                os.fsync(output_handle.fileno())
            if temporary.stat().st_size == 0:
                raise ValueError(f"{step.step_id}: prepared output is empty")
            os.replace(temporary, destination)
        finally:
            try:
                temporary.unlink()
            except FileNotFoundError:
                pass

    def execute_subprocess(self, step: Step) -> None:
        for input_name in step.inputs:
            if not Path(input_name).is_file():
                raise FileNotFoundError(f"{step.step_id}: missing input {input_name}")
        output_map: Dict[str, Path] = {}
        for output_name in step.outputs:
            output = Path(output_name)
            output.parent.mkdir(parents=True, exist_ok=True)
            temporary = output.parent / f".{output.name}.{step.step_id}.{os.getpid()}.tmp"
            try:
                temporary.unlink()
            except FileNotFoundError:
                pass
            output_map[output_name] = temporary
        argv = [str(output_map.get(token, Path(token))) for token in step.argv]
        if argv:
            argv[0] = self.resolved_tools.get(step.argv[0], argv[0])
        stderr_path = self.log_dir / f"{step.step_id}.stderr.log"
        stdout_log_path = self.log_dir / f"{step.step_id}.stdout.log"
        self.log_dir.mkdir(parents=True, exist_ok=True)
        try:
            with stderr_path.open("wb") as stderr_handle:
                if step.stdout_output is not None:
                    with output_map[step.stdout_output].open("wb") as stdout_handle:
                        subprocess.run(argv, check=True, stdout=stdout_handle, stderr=stderr_handle)
                        stdout_handle.flush()
                        os.fsync(stdout_handle.fileno())
                elif len(argv) > 1 and argv[1] == "proali":
                    # proali writes a terminal progress display containing very
                    # large numbers of backspace characters to stdout. On a
                    # whole-genome run this can inflate into multi-gigabyte log
                    # files in minutes, although all durable results are written
                    # through -n, -o, and -f. Keep a provenance note and discard
                    # only that progress stream; stderr remains fully captured.
                    with stdout_log_path.open("wb") as stdout_handle:
                        stdout_handle.write(
                            b"AnchorWave proali terminal progress stream discarded; "
                            b"see stderr log and declared result files.\n"
                        )
                        stdout_handle.flush()
                        os.fsync(stdout_handle.fileno())
                    subprocess.run(
                        argv, check=True, stdout=subprocess.DEVNULL, stderr=stderr_handle
                    )
                else:
                    with stdout_log_path.open("wb") as stdout_handle:
                        subprocess.run(argv, check=True, stdout=stdout_handle, stderr=stderr_handle)
            for final_name, temporary in output_map.items():
                if not temporary.is_file() or temporary.stat().st_size == 0:
                    raise ValueError(f"{step.step_id}: missing or empty output {temporary}")
            for final_name, temporary in output_map.items():
                os.replace(temporary, Path(final_name))
        finally:
            for temporary in output_map.values():
                try:
                    temporary.unlink()
                except FileNotFoundError:
                    pass

    def run_step(self, step: Step) -> None:
        signature = self.signature(step)
        if self.reusable(step, signature):
            print(f"[reuse] {step.step_id}")
            self.event({"event": "reuse", "step": step.step_id, "signature": signature})
            return
        existing_outputs = [path for path in step.outputs if Path(path).exists()]
        if existing_outputs and not self.resume:
            raise FileExistsError(
                f"{step.step_id}: output exists and --resume was not supplied: "
                + ", ".join(existing_outputs)
            )
        print(f"[run] {step.step_id}")
        self.event({"event": "start", "step": step.step_id, "signature": signature})
        try:
            if step.kind in ("gunzip", "copy"):
                self.execute_prepare(step)
            elif step.kind == "subprocess":
                self.execute_subprocess(step)
            else:
                raise ValueError(f"{step.step_id}: unsupported step kind {step.kind!r}")
            outputs = [self.fingerprint(Path(path)) for path in step.outputs]
            receipt = {
                "schema_version": 1,
                "step": step.step_id,
                "signature": signature,
                "outputs": outputs,
                "completed_utc": utc_now(),
            }
            atomic_write_text(
                self.receipt_path(step), json.dumps(receipt, indent=2, sort_keys=True) + "\n"
            )
            self.event({"event": "complete", "step": step.step_id, "signature": signature})
        except BaseException as error:
            self.event(
                {
                    "event": "failure",
                    "step": step.step_id,
                    "signature": signature,
                    "error": f"{type(error).__name__}: {error}",
                }
            )
            raise

    def run(self) -> None:
        atomic_write_text(
            self.plan.run_dir / "provenance" / "execution.json",
            json.dumps(
                {
                    "schema_version": 1,
                    "workflow_version": WORKFLOW_VERSION,
                    "started_utc": utc_now(),
                    "python_executable": sys.executable,
                    "python_version": sys.version,
                    "tools": self.tool_provenance,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
        )
        for step in self.plan.steps:
            self.run_step(step)


def resolve_executable(value: str) -> str:
    if os.sep in value or (os.altsep and os.altsep in value):
        path = Path(value).expanduser().resolve()
        if not path.is_file() or not os.access(path, os.X_OK):
            raise FileNotFoundError(f"executable is missing or not executable: {path}")
        return str(path)
    resolved = shutil.which(value)
    if resolved is None:
        raise FileNotFoundError(f"executable not found on PATH: {value}")
    return str(Path(resolved).resolve())


def materialize_plan(plan: WorkflowPlan, allow_replace: bool) -> None:
    plan.run_dir.mkdir(parents=True, exist_ok=True)
    artifacts = {
        plan.taxa_manifest: plan.taxon_manifest_text(),
        plan.pairwise_manifest: plan.pairwise_manifest_text(),
        plan.trio_command_file: "#!/bin/sh\nset -eu\n" + shlex.join(plan.trio_argv()) + "\n",
        plan.plan_path: json.dumps(plan.to_dict(), indent=2, sort_keys=True) + "\n",
    }
    if plan.tree_text is not None:
        artifacts[plan.tree_path] = plan.tree_text
        artifacts[Path(str(plan.tree_path) + ".metadata.json")] = json.dumps(
            {
                "schema_version": 1,
                "algorithm": TREE_ALGORITHM_VERSION,
                "seed": plan.args.random_tree_seed,
                "samples": sorted(plan.selected),
                "forced_sisters": [plan.args.ingroup1, plan.args.ingroup2],
                "ancestor_name": plan.args.ancestor_node,
                "biological_warning": BIOLOGICAL_WARNING,
                "newick_sha256": hashlib.sha256(plan.tree_text.encode("utf-8")).hexdigest(),
            },
            indent=2,
            sort_keys=True,
        ) + "\n"
    else:
        assert plan.args.species_tree is not None
        source_tree = plan.args.species_tree.resolve()
        if not source_tree.is_file():
            raise FileNotFoundError(f"species tree does not exist: {source_tree}")
        artifacts[plan.tree_path] = source_tree.read_text(encoding="utf-8")

    for path, content in artifacts.items():
        if path.exists():
            existing = path.read_text(encoding="utf-8")
            if existing == content:
                continue
            if not allow_replace:
                raise FileExistsError(
                    f"refusing to replace differing plan artifact {path}; use --replace-plan"
                )
        atomic_write_text(path, content)
    mode = plan.trio_command_file.stat().st_mode
    plan.trio_command_file.chmod(mode | 0o100)
    (plan.run_dir / "trio").mkdir(parents=True, exist_ok=True)


def build_parser() -> argparse.ArgumentParser:
    workflow_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--plan", action="store_true", help="write manifests and plan; run nothing")
    action.add_argument("--dry-run", action="store_true", help="print plan; change nothing")
    action.add_argument("--execute", action="store_true", help="materialize and execute the plan")
    parser.add_argument(
        "--manifest", type=Path, default=workflow_dir / "data" / "maizegdb_assets.tsv"
    )
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--ingroup1", required=True, help="preferred ingroup/proali reference")
    parser.add_argument("--ingroup2", required=True, help="second ingroup and second reference")
    parser.add_argument("--outgroup", action="append", required=True, help="repeatable; first is primary")
    tree_group = parser.add_mutually_exclusive_group(required=True)
    tree_group.add_argument("--species-tree", type=Path)
    tree_group.add_argument("--random-tree-seed", type=int)
    parser.add_argument("--ancestor-node", default="INGROUP_ANCESTOR")
    parser.add_argument("--anchorwave", default="anchorwave")
    parser.add_argument("--minimap2", default="minimap2")
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="legacy fallback used by both tools unless a tool-specific value is set",
    )
    parser.add_argument(
        "--minimap2-threads",
        type=int,
        help="threads for each minimap2 job (defaults to --threads)",
    )
    parser.add_argument(
        "--anchorwave-threads",
        type=int,
        help="threads for each AnchorWave proali job (defaults to --threads)",
    )
    parser.add_argument("--min-exon", type=int, default=20)
    parser.add_argument("--reference-depth", type=int, default=1, metavar="R")
    parser.add_argument("--query-depth", type=int, default=1, metavar="Q")
    parser.add_argument("--expected-gene-copies", type=int, default=1, metavar="E")
    parser.add_argument("--window-width", type=int)
    parser.add_argument("--fa3", type=int)
    parser.add_argument("--no-new-anchors", action="store_true")
    parser.add_argument("--exon-model", action="store_true")
    parser.add_argument("--copy-relations", type=Path)
    parser.add_argument("--copy-mode", choices=("constrained", "strict"), default="constrained")
    parser.add_argument("--copy-number", action="append", default=[])
    parser.add_argument(
        "--block-projections",
        type=Path,
        help="optional extant macro-block projection TSV passed to trioGraphAli",
    )
    parser.add_argument("--resume", action="store_true", help="reuse only matching receipts/outputs")
    parser.add_argument(
        "--skip-source-md5",
        action="store_true",
        help="check pinned byte counts but skip source MD5 (not recommended)",
    )
    parser.add_argument(
        "--skip-source-integrity",
        action="store_true",
        help="check only that source files exist; skip pinned byte counts and MD5",
    )
    parser.add_argument(
        "--replace-plan", action="store_true", help="replace differing plan/manifest artifacts"
    )
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if args.minimap2_threads is None:
        args.minimap2_threads = args.threads
    if args.anchorwave_threads is None:
        args.anchorwave_threads = args.threads
    for name in ("threads", "reference_depth", "query_depth", "expected_gene_copies"):
        if getattr(args, name) < 1:
            raise ValueError(f"--{name.replace('_', '-')} must be at least 1")
    for name in ("minimap2_threads", "anchorwave_threads"):
        if getattr(args, name) < 1:
            raise ValueError(f"--{name.replace('_', '-')} must be at least 1")
    if args.min_exon < 0:
        raise ValueError("--min-exon must be non-negative")
    if args.window_width is not None and args.window_width < 1:
        raise ValueError("--window-width must be positive")
    if args.fa3 is not None and args.fa3 < 1:
        raise ValueError("--fa3 must be positive")
    for sample in [args.ingroup1, args.ingroup2] + list(args.outgroup):
        require_safe_component(sample, "selected sample ID")
    if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_.-]*", args.ancestor_node):
        raise ValueError("--ancestor-node must be a Newick-safe identifier")
    for attribute in (
        "manifest",
        "data_root",
        "run_dir",
        "species_tree",
        "copy_relations",
        "block_projections",
    ):
        value = getattr(args, attribute)
        if value is not None:
            require_no_field_controls(str(value), f"--{attribute.replace('_', '-')}")
    require_no_field_controls(args.anchorwave, "--anchorwave")
    require_no_field_controls(args.minimap2, "--minimap2")
    for specification in args.copy_number:
        require_no_field_controls(specification, "--copy-number")
    if args.resume and not args.execute:
        raise ValueError("--resume is meaningful only with --execute")
    if args.skip_source_md5 and not args.execute:
        raise ValueError("--skip-source-md5 is meaningful only with --execute")
    if args.skip_source_integrity and not args.execute:
        raise ValueError("--skip-source-integrity is meaningful only with --execute")


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        validate_args(args)
        args.manifest = args.manifest.resolve()
        catalog = load_catalog(args.manifest)
        for executable_attribute in ("anchorwave", "minimap2"):
            value = getattr(args, executable_attribute)
            if os.sep in value or (os.altsep and os.altsep in value):
                setattr(args, executable_attribute, str(Path(value).expanduser().resolve()))
        resolved_tools: Dict[str, str] = {}
        if args.execute:
            # Validate now, but keep a bare command name in the plan so a
            # preceding --plan and later --execute are byte-for-byte stable.
            resolved_tools[args.anchorwave] = resolve_executable(args.anchorwave)
            resolved_tools[args.minimap2] = resolve_executable(args.minimap2)
        plan = WorkflowPlan(args, catalog)
        if args.dry_run:
            print(json.dumps(plan.to_dict(), indent=2, sort_keys=True))
            print("\n# trioGraphAli (emitted, not executed)")
            print(shlex.join(plan.trio_argv()))
            if plan.tree_text is not None:
                print(f"\n# {BIOLOGICAL_WARNING}")
            return 0

        materialize_plan(plan, allow_replace=args.replace_plan or args.resume)
        print(f"wrote workflow plan: {plan.plan_path}")
        print(f"pairwise jobs: {len(plan.pairs)} (= 1 + 2 * {len(args.outgroup)})")
        if plan.tree_text is not None:
            print(BIOLOGICAL_WARNING)
        if args.plan:
            print("plan only: no data preparation or alignment command was executed")
            return 0

        Executor(
            plan,
            resume=args.resume,
            verify_size=not args.skip_source_integrity,
            verify_md5=not args.skip_source_md5 and not args.skip_source_integrity,
            resolved_tools=resolved_tools,
        ).run()
        print(f"pairwise evidence ready: {plan.pairwise_manifest}")
        print(f"next command: {shlex.join(plan.trio_argv())}")
        return 0
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f"workflow error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
