#!/usr/bin/env python3
from __future__ import annotations

import contextlib
import io
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


WORKFLOW_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WORKFLOW_DIR))

import download_maizegdb as downloader  # noqa: E402


class MaizeGdbManifestTests(unittest.TestCase):
    def test_pinned_catalog_has_expected_cohorts_and_asset_count(self) -> None:
        records = downloader.load_manifest(downloader.DEFAULT_MANIFEST)
        self.assertEqual(35, len(records))
        self.assertEqual(26, sum(record.cohort == "NAM" for record in records))
        self.assertEqual(9, sum(record.cohort == "PanAnd" for record in records))
        assets = downloader.records_to_assets(records, "both")
        self.assertEqual(69, len(assets))
        self.assertEqual(
            1,
            sum(asset.md5 == downloader.MISSING_CHECKSUM for asset in assets),
        )

    def test_partial_bytes_reduce_disk_preflight_requirement(self) -> None:
        asset = self._asset("I1", 100)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            final = downloader.asset_path(root, asset)
            final.parent.mkdir(parents=True)
            Path(str(final) + ".part").write_bytes(b"x" * 40)
            self.assertEqual(60, downloader.remaining_bytes(root, [asset]))

    @staticmethod
    def _asset(sample: str, size: int = 10) -> downloader.Asset:
        return downloader.Asset(
            sample_id=sample,
            cohort="NAM",
            kind="assembly",
            url=f"https://download.maizegdb.org/Genomes/NAM_Founders/{sample}.fa.gz",
            filename=f"{sample}.fa.gz",
            md5="0" * 32,
            size=size,
            status="official",
        )


class MaizeGdbDownloadSchedulingTests(unittest.TestCase):
    def test_bounded_parallel_scheduler_visits_each_asset_once(self) -> None:
        assets = [MaizeGdbManifestTests._asset(f"I{index}") for index in range(5)]
        with tempfile.TemporaryDirectory() as temporary:
            with mock.patch.object(downloader, "download_asset") as worker:
                downloader.download_assets(
                    assets, Path(temporary), retries=3, timeout=5.0, jobs=3
                )
        self.assertEqual(5, worker.call_count)
        self.assertEqual({asset.sample_id for asset in assets}, {
            call.args[0].sample_id for call in worker.call_args_list
        })

    def test_parallel_failures_are_aggregated(self) -> None:
        assets = [MaizeGdbManifestTests._asset("I1")]
        with tempfile.TemporaryDirectory() as temporary:
            with mock.patch.object(
                downloader, "download_asset", side_effect=downloader.DownloadError("boom")
            ):
                with contextlib.redirect_stderr(io.StringIO()):
                    with self.assertRaisesRegex(
                        downloader.DownloadError, "1 asset.*I1/assembly.*boom"
                    ):
                        downloader.download_assets(
                            assets, Path(temporary), retries=1, timeout=1.0, jobs=2
                        )

    def test_jobs_are_capped_for_server_politeness(self) -> None:
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                downloader.build_parser().parse_args(["--list", "--jobs", "17"])


if __name__ == "__main__":
    unittest.main()
