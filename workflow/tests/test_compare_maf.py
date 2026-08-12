import io
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

from workflow.compare_maf import ScoringParameters, compare_mafs, main


class CompareMafTests(unittest.TestCase):
    def write_maf(self, directory: Path, name: str, body: str) -> Path:
        path = directory / name
        path.write_text(body, encoding="ascii")
        return path

    def test_exact_mafs(self):
        maf = """##maf version=1
a score=-10
s ref 0 4 + 4 AC-GT
s qry 0 4 + 4 A-CGT

"""
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            old = self.write_maf(directory, "old.maf", maf)
            new = self.write_maf(directory, "new.maf", maf)
            report = compare_mafs(old, new)

        self.assertEqual("EXACT", report["verdict"])
        self.assertEqual(1, report["comparison"]["exact_blocks"])

    def test_accepts_tab_delimited_anchorwave_rows(self):
        maf = """a\tscore=0
s\tref\t0\t4\t+\t4\tACGT
s\tqry\t0\t4\t+\t4\tACGT

"""
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            old = self.write_maf(directory, "old.maf", maf)
            new = self.write_maf(directory, "new.maf", maf)
            report = compare_mafs(old, new)

        self.assertEqual("EXACT", report["verdict"])
        self.assertEqual(2, report["old"]["rows"])

    def test_tie_resolved_cigars_are_coordinate_equivalent(self):
        old_maf = """a score=-4
s ref 0 4 + 4 A-AAA
s qry 0 5 + 5 AAAAA

"""
        new_maf = """a score=-4
s ref 0 4 + 4 AAA-A
s qry 0 5 + 5 AAAAA

"""
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            old = self.write_maf(directory, "old.maf", old_maf)
            new = self.write_maf(directory, "new.maf", new_maf)
            report = compare_mafs(old, new)

        self.assertEqual("COORDINATE_EQUIVALENT", report["verdict"])
        self.assertEqual(1, report["comparison"]["alignment_text_mismatches"])
        self.assertEqual(0, report["comparison"]["coordinate_mismatches"])

    def test_coordinate_change_is_different(self):
        old_maf = """a score=0
s ref 0 4 + 10 ACGT
s qry 0 4 + 10 ACGT

"""
        new_maf = """a score=0
s ref 1 4 + 10 ACGT
s qry 0 4 + 10 ACGT

"""
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            old = self.write_maf(directory, "old.maf", old_maf)
            new = self.write_maf(directory, "new.maf", new_maf)
            report = compare_mafs(old, new)

        self.assertEqual("DIFFERENT", report["verdict"])
        self.assertEqual(1, report["comparison"]["coordinate_mismatches"])

    def test_invalid_maf_is_reported(self):
        invalid = """a score=0
s ref 0 3 + 4 ACGT
s qry 0 4 + 4 ACGT

"""
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            old = self.write_maf(directory, "old.maf", invalid)
            new = self.write_maf(directory, "new.maf", invalid)
            report = compare_mafs(old, new)

        self.assertEqual("INVALID_MAF", report["verdict"])

    def test_report_only_keeps_valid_differences_nonfatal(self):
        old_maf = """a score=0
s ref 0 4 + 10 ACGT
s qry 0 4 + 10 ACGT

"""
        new_maf = """a score=0
s ref 1 4 + 10 ACGT
s qry 0 4 + 10 ACGT

"""
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            old = self.write_maf(directory, "old.maf", old_maf)
            new = self.write_maf(directory, "new.maf", new_maf)
            report_path = directory / "report.json"
            with redirect_stdout(io.StringIO()):
                strict_status = main([str(old), str(new)])
                report_status = main([
                    str(old), str(new), "--json", str(report_path),
                    "--report-only"])

            self.assertEqual(2, strict_status)
            self.assertEqual(0, report_status)
            self.assertIn('"verdict": "DIFFERENT"',
                          report_path.read_text(encoding="utf-8"))

    def test_score_audit_detects_an_incorrect_declared_fallback_score(self):
        old_maf = """a score=-1
s ref 0 4 + 4 ACGT
s qry 0 4 + 4 TGCA

"""
        new_maf = """a score=-24
s ref 0 4 + 4 ACGT
s qry 0 4 + 4 TGCA

"""
        scoring = ScoringParameters(6, 8, 2, 75, 1)
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            old = self.write_maf(directory, "old.maf", old_maf)
            new = self.write_maf(directory, "new.maf", new_maf)
            report = compare_mafs(old, new, scoring=scoring)

        audit = report["score_audit"]
        self.assertEqual(1, audit["audited_paired_blocks"])
        self.assertEqual(1, audit["recomputed_score_equal_blocks"])
        self.assertEqual(1, audit["old_declared_score_inconsistencies"])
        self.assertEqual(0, audit["new_declared_score_inconsistencies"])


if __name__ == "__main__":
    unittest.main()
