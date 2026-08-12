#!/usr/bin/env python3
from __future__ import annotations

import contextlib
import csv
import gzip
import hashlib
import io
import json
import os
import shlex
import stat
import sys
import tempfile
import unittest
from pathlib import Path


WORKFLOW_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WORKFLOW_DIR))

import random_newick  # noqa: E402
import trioanchorgraph_workflow as workflow  # noqa: E402


CATALOG_FIELDS = [
    "sample_id",
    "cohort",
    "taxon",
    "assembly_dir",
    "assembly_url",
    "assembly_filename",
    "assembly_md5",
    "assembly_size_bytes",
    "annotation_url",
    "annotation_filename",
    "annotation_md5",
    "annotation_size_bytes",
    "annotation_status",
    "model",
    "accession",
    "notice_url",
]


class RandomNewickTests(unittest.TestCase):
    def test_tree_is_seeded_and_forced_ingroup_node_is_named(self) -> None:
        samples = ["B73", "B97", "Zv-TIL01", "Zx-TIL25"]
        first = random_newick.generate_random_newick(
            samples, 20260806, ("B73", "B97"), "ANCESTOR"
        )
        second = random_newick.generate_random_newick(
            reversed(samples), 20260806, ("B73", "B97"), "ANCESTOR"
        )
        different = random_newick.generate_random_newick(
            samples, 20260807, ("B73", "B97"), "ANCESTOR"
        )
        self.assertEqual(first, second)
        self.assertNotEqual(first, different)
        self.assertIn(")ANCESTOR", first)
        self.assertTrue(first.endswith("ROOT;\n"))
        for sample in samples:
            self.assertEqual(first.count(sample), 1)

    def test_two_forced_sisters_have_one_valid_internal_label(self) -> None:
        tree = random_newick.generate_random_newick(
            ["I1", "I2"], 1, ("I1", "I2"), "ANCESTOR"
        )
        self.assertTrue(tree.endswith(")ANCESTOR;\n"))
        self.assertNotIn("ANCESTORROOT", tree)


class WorkflowTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.data_root = self.root / "downloads"
        self.run_dir = self.root / "run"
        self.catalog_path = self.root / "assets.tsv"
        self.tool_log = self.root / "fake-tools.log"
        self._records = []

    def tearDown(self) -> None:
        self.temporary.cleanup()

    @staticmethod
    def _md5(data: bytes) -> str:
        return hashlib.md5(data).hexdigest()

    def add_sample(self, sample: str, cohort: str, annotated: bool) -> None:
        assembly_filename = f"{sample}.fa.gz"
        assembly_raw = f">chr1\nACGTACGT{sample}\n".encode()
        assembly_path = self.data_root / cohort / sample / assembly_filename
        assembly_path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.GzipFile(filename=str(assembly_path), mode="wb", mtime=0) as handle:
            handle.write(assembly_raw)
        assembly_data = assembly_path.read_bytes()

        annotation_filename = f"{sample}.gff3.gz" if annotated else ""
        annotation_md5 = ""
        annotation_size = ""
        if annotated:
            annotation_path = self.data_root / cohort / sample / annotation_filename
            with gzip.GzipFile(filename=str(annotation_path), mode="wb", mtime=0) as handle:
                handle.write(
                    b"##gff-version 3\nchr1\ttest\tgene\t1\t4\t.\t+\t.\tID=g1\n"
                )
            annotation_data = annotation_path.read_bytes()
            annotation_md5 = self._md5(annotation_data)
            annotation_size = str(len(annotation_data))

        self._records.append(
            {
                "sample_id": sample,
                "cohort": cohort,
                "taxon": "Zea test",
                "assembly_dir": sample,
                "assembly_url": f"https://invalid.example/{assembly_filename}",
                "assembly_filename": assembly_filename,
                "assembly_md5": self._md5(assembly_data),
                "assembly_size_bytes": str(len(assembly_data)),
                "annotation_url": (
                    f"https://invalid.example/{annotation_filename}" if annotated else ""
                ),
                "annotation_filename": annotation_filename,
                "annotation_md5": annotation_md5,
                "annotation_size_bytes": annotation_size,
                "annotation_status": "official" if annotated else "assembly_only",
                "model": "test.1" if annotated else "",
                "accession": f"TEST_{sample}",
                "notice_url": "https://invalid.example/notice",
            }
        )

    def write_catalog(self) -> None:
        with self.catalog_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=CATALOG_FIELDS, delimiter="\t")
            writer.writeheader()
            writer.writerows(self._records)

    def make_fake_tools(self) -> tuple[Path, Path]:
        anchorwave = self.root / "fake-anchorwave"
        minimap2 = self.root / "fake-minimap2"
        anchorwave.write_text(
            """#!/usr/bin/env python3
import os, pathlib, sys
with open(os.environ['FAKE_TOOL_LOG'], 'a') as h: h.write('anchorwave ' + ' '.join(sys.argv[1:]) + '\\n')
args = sys.argv[1:]
def value(flag): return pathlib.Path(args[args.index(flag) + 1])
if args[0] == 'gff2seq':
    value('-o').write_text('>anchor\\nACGT\\n')
elif args[0] == 'proali':
    sys.stdout.write('progress\\b' * 1000)
    value('-n').write_text('# anchors\\nanchor1\\n')
    value('-o').write_text('##maf version=1\\n\\na score=1\\ns ref.chr1 0 4 + 4 ACGT\\ns qry.chr1 0 4 + 4 ACGT\\n')
    value('-f').write_text('##maf version=1\\n')
else:
    raise SystemExit(3)
""",
            encoding="utf-8",
        )
        minimap2.write_text(
            """#!/usr/bin/env python3
import os, sys
with open(os.environ['FAKE_TOOL_LOG'], 'a') as h: h.write('minimap2 ' + ' '.join(sys.argv[1:]) + '\\n')
print('@HD\\tVN:1.6')
print('anchor\\t4\\t*\\t0\\t0\\t*\\t*\\t0\\t0\\tACGT\\t*')
""",
            encoding="utf-8",
        )
        anchorwave.chmod(anchorwave.stat().st_mode | stat.S_IXUSR)
        minimap2.chmod(minimap2.stat().st_mode | stat.S_IXUSR)
        return anchorwave, minimap2

    def base_arguments(self) -> list[str]:
        return [
            "--manifest",
            str(self.catalog_path),
            "--data-root",
            str(self.data_root),
            "--run-dir",
            str(self.run_dir),
            "--ingroup1",
            "I1",
            "--ingroup2",
            "I2",
            "--outgroup",
            "O1",
            "--random-tree-seed",
            "17",
            "--ancestor-node",
            "ANCESTOR",
        ]

    def test_multiple_outgroups_create_one_plus_two_k_pairs(self) -> None:
        self.add_sample("I1", "NAM", True)
        self.add_sample("I2", "NAM", True)
        self.add_sample("O1", "PanAnd", False)
        self.add_sample("O2", "PanAnd", False)
        self.write_catalog()
        parser = workflow.build_parser()
        args = parser.parse_args(
            ["--dry-run"] + self.base_arguments()[2:] + ["--outgroup", "O2"]
        )
        # base_arguments()[2:] intentionally drops duplicate --manifest, so add it.
        args.manifest = self.catalog_path
        workflow.validate_args(args)
        plan = workflow.WorkflowPlan(args, workflow.load_catalog(self.catalog_path))
        self.assertEqual(len(plan.pairs), 5)
        self.assertEqual(
            [(pair.reference, pair.query) for pair in plan.pairs],
            [("I1", "I2"), ("I1", "O1"), ("I2", "O1"), ("I1", "O2"), ("I2", "O2")],
        )
        self.assertFalse(self.run_dir.exists(), "constructing/dry-running a plan must not write")

    def test_current_trio_cli_defaults_block_passthrough_and_shell_quoting(self) -> None:
        self.add_sample("I1", "NAM", True)
        self.add_sample("I2", "NAM", True)
        self.add_sample("O1", "PanAnd", False)
        self.write_catalog()
        self.run_dir = self.root / "run with spaces and 'apostrophe'"
        block_projections = self.root / "macro blocks 'v1'.tsv"
        block_projections.write_text("#anchorwave-block-projections-v1\n", encoding="utf-8")
        parser = workflow.build_parser()
        args = parser.parse_args(
            ["--plan"]
            + self.base_arguments()
            + ["--block-projections", str(block_projections)]
        )
        workflow.validate_args(args)
        plan = workflow.WorkflowPlan(args, workflow.load_catalog(self.catalog_path))
        argv = plan.trio_argv()

        self.assertIn("--validate-input-paths", argv)
        self.assertNotIn("--skip-fasta-validation", argv)
        block_index = argv.index("--block-projections")
        self.assertEqual(argv[block_index + 1], str(block_projections.resolve()))
        self.assertEqual(
            plan.to_dict()["trioGraphAli_validation"],
            {"fasta_evidence": "enabled_by_cli_default", "manifest_input_paths": True},
        )

        workflow.materialize_plan(plan, allow_replace=False)
        script_lines = plan.trio_command_file.read_text(encoding="utf-8").splitlines()
        self.assertEqual(script_lines[:2], ["#!/bin/sh", "set -eu"])
        # shlex splitting is the shell-tokenization round trip: neither spaces
        # nor apostrophes in real paths can add or alter an argument.
        self.assertEqual(shlex.split(script_lines[2], posix=True), argv)
        self.assertIn("'\"'\"'", script_lines[2])

    def test_manifest_breaking_control_characters_are_rejected(self) -> None:
        self.add_sample("I1", "NAM", True)
        self.add_sample("I2", "NAM", True)
        self.add_sample("O1", "PanAnd", False)
        self.write_catalog()
        parser = workflow.build_parser()
        arguments = self.base_arguments()
        arguments[arguments.index("--run-dir") + 1] = str(self.root / "bad\tpath")
        args = parser.parse_args(["--dry-run"] + arguments)
        with self.assertRaisesRegex(ValueError, "contains NUL, tab"):
            workflow.validate_args(args)

    def test_catalog_path_components_reject_directory_traversal(self) -> None:
        self.add_sample("I1", "NAM", True)
        self._records[0]["cohort"] = ".."
        self.write_catalog()
        with self.assertRaisesRegex(ValueError, "non-traversing"):
            workflow.load_catalog(self.catalog_path)

    def test_execute_and_receipt_validated_resume_with_fake_tools(self) -> None:
        self.add_sample("I1", "NAM", True)
        self.add_sample("I2", "NAM", True)
        self.add_sample("O1", "PanAnd", False)
        self.write_catalog()
        anchorwave, minimap2 = self.make_fake_tools()
        old_log = os.environ.get("FAKE_TOOL_LOG")
        os.environ["FAKE_TOOL_LOG"] = str(self.tool_log)
        try:
            arguments = (
                ["--execute"]
                + self.base_arguments()
                + [
                    "--anchorwave",
                    str(anchorwave),
                    "--minimap2",
                    str(minimap2),
                    "--minimap2-threads",
                    "12",
                    "--anchorwave-threads",
                    "6",
                    "--skip-source-integrity",
                ]
            )
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(workflow.main(arguments), 0)
            first_log = self.tool_log.read_text(encoding="utf-8")
            self.assertEqual(len(first_log.splitlines()), 10)
            self.assertTrue((self.run_dir / "pairwise/I1__I2/alignment.maf").is_file())
            progress_log = self.run_dir / "provenance/logs/proali.I1.to.I2.stdout.log"
            self.assertLess(progress_log.stat().st_size, 256)
            self.assertIn("progress stream discarded", progress_log.read_text())
            pairwise = (self.run_dir / "manifests/pairwise.tsv").read_text(encoding="utf-8")
            self.assertEqual(len(pairwise.strip().splitlines()), 4)
            plan_json = json.loads((self.run_dir / "workflow.plan.json").read_text())
            self.assertEqual(plan_json["pair_count"], 3)
            self.assertEqual(
                plan_json["execution_policy"],
                {
                    "step_concurrency": 1,
                    "pairwise_concurrency": 1,
                    "minimap2_threads": 12,
                    "anchorwave_threads": 6,
                    "source_size_validation": False,
                    "source_md5_validation": False,
                },
            )
            for line in first_log.splitlines():
                if line.startswith("minimap2 "):
                    self.assertIn("-t 12", line)
                if line.startswith("anchorwave proali "):
                    self.assertIn("-t 6", line)
            self.assertIn("R1:Q1:m20:e1", pairwise)
            self.assertIn("RANDOM TREE", (self.run_dir / "manifests/species_tree.nwk.metadata.json").read_text())

            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(workflow.main(arguments + ["--resume"]), 0)
            self.assertEqual(self.tool_log.read_text(encoding="utf-8"), first_log)
            receipts = list((self.run_dir / "provenance/receipts").glob("*.json"))
            self.assertEqual(len(receipts), 15)  # 5 prepare + 2 gff2seq + 5 minimap2 + 3 proali

            # A changed executable invalidates its five mapping receipts; the
            # resulting SAM mtimes then invalidate the three dependent proali
            # receipts. Preparation and gff2seq remain reusable.
            minimap2.write_text(
                minimap2.read_text(encoding="utf-8") + "\n# changed tool build\n",
                encoding="utf-8",
            )
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(workflow.main(arguments + ["--resume"]), 0)
            self.assertEqual(len(self.tool_log.read_text(encoding="utf-8").splitlines()), 18)
        finally:
            if old_log is None:
                os.environ.pop("FAKE_TOOL_LOG", None)
            else:
                os.environ["FAKE_TOOL_LOG"] = old_log


if __name__ == "__main__":
    unittest.main()
