#!/usr/bin/env python3
"""Reproducible, checksum-aware downloader for the pinned MaizeGDB manifest.

The program has no third-party dependencies and requires an explicit action.
In particular, importing it or invoking it without ``--download`` never starts
a data transfer.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import http.client
import json
import os
import re
import shutil
import socket
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


DEFAULT_MANIFEST = Path(__file__).resolve().parent / "data" / "maizegdb_assets.tsv"
ALLOWED_HOST = "download.maizegdb.org"
MISSING_CHECKSUM = "MISSING_IN_SOURCE_MANIFEST"
MD5_RE = re.compile(r"^[0-9a-f]{32}$")
SAFE_COMPONENT_RE = re.compile(r"^[A-Za-z0-9._-]+$")
CHUNK_SIZE = 4 * 1024 * 1024
MIN_FREE_MARGIN = 256 * 1024 * 1024
USER_AGENT = "AnchorWave-MaizeGDB-workflow/1.0"

REQUIRED_COLUMNS = (
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
)


class ManifestError(ValueError):
    """The pinned TSV is malformed or violates a safety invariant."""


class DownloadError(RuntimeError):
    """A transfer failed or produced an incomplete file."""


class IntegrityError(DownloadError):
    """A file failed its expected size or checksum."""


@dataclass(frozen=True)
class Record:
    sample_id: str
    cohort: str
    taxon: str
    assembly_dir: str
    assembly_url: str
    assembly_filename: str
    assembly_md5: str
    assembly_size: int
    annotation_url: str
    annotation_filename: str
    annotation_md5: str
    annotation_size: Optional[int]
    annotation_status: str
    model: str
    accession: str
    notice_url: str


@dataclass(frozen=True)
class Asset:
    sample_id: str
    cohort: str
    kind: str
    url: str
    filename: str
    md5: str
    size: int
    status: str


def parse_positive_int(text: str, field: str, row_number: int) -> Optional[int]:
    if not text:
        return None
    try:
        value = int(text)
    except ValueError as exc:
        raise ManifestError(
            f"row {row_number}: {field} must be an integer, got {text!r}"
        ) from exc
    if value <= 0:
        raise ManifestError(f"row {row_number}: {field} must be positive")
    return value


def validate_component(value: str, field: str, row_number: int) -> None:
    if not SAFE_COMPONENT_RE.fullmatch(value):
        raise ManifestError(
            f"row {row_number}: unsafe {field} {value!r}; expected letters, digits, '.', '_' or '-'"
        )


def validate_url(url: str, filename: str, field: str, row_number: int) -> None:
    parsed = urllib.parse.urlsplit(url)
    if parsed.scheme != "https" or parsed.hostname != ALLOWED_HOST:
        raise ManifestError(
            f"row {row_number}: {field} must be HTTPS on {ALLOWED_HOST}: {url!r}"
        )
    if parsed.username or parsed.password or parsed.query or parsed.fragment:
        raise ManifestError(f"row {row_number}: unsupported URL components in {field}")
    remote_name = urllib.parse.unquote(Path(parsed.path).name)
    if remote_name != filename:
        raise ManifestError(
            f"row {row_number}: {field} basename {remote_name!r} does not match {filename!r}"
        )


def validate_md5(value: str, field: str, row_number: int, allow_missing: bool) -> None:
    if allow_missing and value == MISSING_CHECKSUM:
        return
    if not MD5_RE.fullmatch(value):
        raise ManifestError(f"row {row_number}: invalid {field} value {value!r}")


def load_manifest(path: Path) -> List[Record]:
    records: List[Record] = []
    seen = set()
    try:
        handle = path.open("r", encoding="utf-8", newline="")
    except OSError as exc:
        raise ManifestError(f"cannot open manifest {path}: {exc}") from exc

    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ManifestError(f"manifest {path} has no header")
        missing = [name for name in REQUIRED_COLUMNS if name not in reader.fieldnames]
        if missing:
            raise ManifestError(f"manifest is missing columns: {', '.join(missing)}")

        for row_number, raw in enumerate(reader, start=2):
            row = {key: (raw.get(key) or "").strip() for key in REQUIRED_COLUMNS}
            sample_id = row["sample_id"]
            cohort = row["cohort"]
            validate_component(sample_id, "sample_id", row_number)
            validate_component(cohort, "cohort", row_number)
            validate_component(row["assembly_dir"], "assembly_dir", row_number)
            validate_component(row["assembly_filename"], "assembly_filename", row_number)
            if sample_id in seen:
                raise ManifestError(f"row {row_number}: duplicate sample_id {sample_id!r}")
            seen.add(sample_id)
            if cohort not in {"NAM", "PanAnd"}:
                raise ManifestError(f"row {row_number}: unsupported cohort {cohort!r}")
            if not row["taxon"] or not row["accession"]:
                raise ManifestError(f"row {row_number}: taxon and accession are required")

            assembly_size = parse_positive_int(
                row["assembly_size_bytes"], "assembly_size_bytes", row_number
            )
            if assembly_size is None:
                raise ManifestError(f"row {row_number}: assembly_size_bytes is required")
            validate_url(
                row["assembly_url"], row["assembly_filename"], "assembly_url", row_number
            )
            validate_md5(row["assembly_md5"], "assembly_md5", row_number, False)

            status = row["annotation_status"]
            if status == "assembly_only":
                annotation_fields = (
                    "annotation_url",
                    "annotation_filename",
                    "annotation_md5",
                    "annotation_size_bytes",
                    "model",
                )
                populated = [field for field in annotation_fields if row[field]]
                if populated:
                    raise ManifestError(
                        f"row {row_number}: assembly_only row has annotation fields: {populated}"
                    )
                annotation_size = None
            elif status in {"official", "official_checksum_missing"}:
                if not all(
                    row[field]
                    for field in (
                        "annotation_url",
                        "annotation_filename",
                        "annotation_md5",
                        "annotation_size_bytes",
                        "model",
                    )
                ):
                    raise ManifestError(
                        f"row {row_number}: {status} annotation is missing required fields"
                    )
                validate_component(
                    row["annotation_filename"], "annotation_filename", row_number
                )
                validate_url(
                    row["annotation_url"],
                    row["annotation_filename"],
                    "annotation_url",
                    row_number,
                )
                validate_md5(
                    row["annotation_md5"], "annotation_md5", row_number, True
                )
                if status == "official_checksum_missing":
                    if row["annotation_md5"] != MISSING_CHECKSUM:
                        raise ManifestError(
                            f"row {row_number}: checksum-missing status requires {MISSING_CHECKSUM}"
                        )
                elif row["annotation_md5"] == MISSING_CHECKSUM:
                    raise ManifestError(
                        f"row {row_number}: missing checksum requires official_checksum_missing status"
                    )
                annotation_size = parse_positive_int(
                    row["annotation_size_bytes"],
                    "annotation_size_bytes",
                    row_number,
                )
                if annotation_size is None:
                    raise ManifestError(
                        f"row {row_number}: annotation_size_bytes is required"
                    )
            else:
                raise ManifestError(
                    f"row {row_number}: unsupported annotation_status {status!r}"
                )

            records.append(
                Record(
                    sample_id=sample_id,
                    cohort=cohort,
                    taxon=row["taxon"],
                    assembly_dir=row["assembly_dir"],
                    assembly_url=row["assembly_url"],
                    assembly_filename=row["assembly_filename"],
                    assembly_md5=row["assembly_md5"],
                    assembly_size=assembly_size,
                    annotation_url=row["annotation_url"],
                    annotation_filename=row["annotation_filename"],
                    annotation_md5=row["annotation_md5"],
                    annotation_size=annotation_size,
                    annotation_status=status,
                    model=row["model"],
                    accession=row["accession"],
                    notice_url=row["notice_url"],
                )
            )

    if not records:
        raise ManifestError(f"manifest {path} contains no data rows")
    return records


def select_records(
    records: Sequence[Record], cohorts: Sequence[str], samples: Sequence[str]
) -> List[Record]:
    known_samples = {record.sample_id for record in records}
    unknown = sorted(set(samples) - known_samples)
    if unknown:
        raise ManifestError(f"unknown sample_id(s): {', '.join(unknown)}")
    cohort_set = set(cohorts)
    sample_set = set(samples)
    return [
        record
        for record in records
        if (not cohort_set or record.cohort in cohort_set)
        and (not sample_set or record.sample_id in sample_set)
    ]


def records_to_assets(records: Iterable[Record], asset_choice: str) -> List[Asset]:
    assets: List[Asset] = []
    include_assembly = asset_choice in {"assembly", "both"}
    include_annotation = asset_choice in {"annotation", "both"}
    for record in records:
        if include_assembly:
            assets.append(
                Asset(
                    sample_id=record.sample_id,
                    cohort=record.cohort,
                    kind="assembly",
                    url=record.assembly_url,
                    filename=record.assembly_filename,
                    md5=record.assembly_md5,
                    size=record.assembly_size,
                    status="official",
                )
            )
        if include_annotation and record.annotation_status != "assembly_only":
            assert record.annotation_size is not None
            assets.append(
                Asset(
                    sample_id=record.sample_id,
                    cohort=record.cohort,
                    kind="annotation",
                    url=record.annotation_url,
                    filename=record.annotation_filename,
                    md5=record.annotation_md5,
                    size=record.annotation_size,
                    status=record.annotation_status,
                )
            )
    return assets


def asset_path(output_dir: Path, asset: Asset) -> Path:
    return output_dir / asset.cohort / asset.sample_id / asset.filename


def human_bytes(value: int) -> str:
    units = ("B", "KiB", "MiB", "GiB", "TiB")
    amount = float(value)
    for unit in units:
        if amount < 1024.0 or unit == units[-1]:
            return f"{amount:.2f} {unit}"
        amount /= 1024.0
    raise AssertionError("unreachable")


def print_assets(assets: Sequence[Asset], output_dir: Path) -> None:
    print(
        "sample_id\tcohort\tasset\tsize_bytes\tmd5\tstatus\turl\tdestination"
    )
    for asset in assets:
        print(
            "\t".join(
                (
                    asset.sample_id,
                    asset.cohort,
                    asset.kind,
                    str(asset.size),
                    asset.md5,
                    asset.status,
                    asset.url,
                    str(asset_path(output_dir, asset)),
                )
            )
        )


def file_md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_file(path: Path, asset: Asset) -> Tuple[bool, str]:
    if not path.is_file():
        return False, "missing"
    actual_size = path.stat().st_size
    if actual_size != asset.size:
        return False, f"size mismatch: expected {asset.size}, found {actual_size}"
    if asset.md5 == MISSING_CHECKSUM:
        return True, "size verified; source manifest has no MD5"
    actual_md5 = file_md5(path)
    if actual_md5 != asset.md5:
        return False, f"MD5 mismatch: expected {asset.md5}, found {actual_md5}"
    return True, "size and MD5 verified"


def nearest_existing_parent(path: Path) -> Path:
    candidate = path.absolute()
    while not candidate.exists():
        parent = candidate.parent
        if parent == candidate:
            raise DownloadError(f"cannot find an existing parent for {path}")
        candidate = parent
    return candidate


def remaining_bytes(output_dir: Path, assets: Sequence[Asset]) -> int:
    total = 0
    for asset in assets:
        final = asset_path(output_dir, asset)
        if final.is_file() and final.stat().st_size == asset.size:
            continue
        part = Path(str(final) + ".part")
        partial_size = part.stat().st_size if part.is_file() else 0
        if partial_size > asset.size:
            raise IntegrityError(
                f"partial file is larger than expected: {part} ({partial_size} > {asset.size})"
            )
        total += asset.size - partial_size
    return total


def disk_preflight(output_dir: Path, assets: Sequence[Asset], create: bool) -> Tuple[int, int, int]:
    if create:
        output_dir.mkdir(parents=True, exist_ok=True)
        anchor = output_dir
    else:
        anchor = nearest_existing_parent(output_dir)
    required = remaining_bytes(output_dir, assets)
    margin = max(MIN_FREE_MARGIN, int(required * 0.02)) if required else 0
    free = shutil.disk_usage(anchor).free
    if free < required + margin:
        raise DownloadError(
            "insufficient disk space: "
            f"need {human_bytes(required)} plus {human_bytes(margin)} margin; "
            f"only {human_bytes(free)} free on {anchor}"
        )
    return required, margin, free


def meta_paths(final: Path) -> Tuple[Path, Path]:
    part = Path(str(final) + ".part")
    meta = Path(str(final) + ".part.meta.json")
    return part, meta


def load_resume_meta(meta_path: Path) -> Dict[str, object]:
    try:
        with meta_path.open("r", encoding="utf-8") as handle:
            value = json.load(handle)
    except (OSError, ValueError) as exc:
        raise DownloadError(f"cannot read resume metadata {meta_path}: {exc}") from exc
    if not isinstance(value, dict):
        raise DownloadError(f"resume metadata is not an object: {meta_path}")
    return value


def write_resume_meta(meta_path: Path, value: Dict[str, object]) -> None:
    temp_path = Path(str(meta_path) + ".tmp")
    with temp_path.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temp_path, meta_path)


def expected_meta(asset: Asset) -> Dict[str, object]:
    return {"url": asset.url, "expected_md5": asset.md5, "expected_size": asset.size}


def validate_resume_meta(meta: Dict[str, object], asset: Asset, meta_path: Path) -> None:
    expected = expected_meta(asset)
    mismatches = [key for key, value in expected.items() if meta.get(key) != value]
    if mismatches:
        raise DownloadError(
            f"resume metadata {meta_path} does not match the selected asset ({', '.join(mismatches)})"
        )


def response_status(response: object) -> int:
    status = getattr(response, "status", None)
    if status is None:
        status = response.getcode()  # type: ignore[attr-defined]
    return int(status)


def parse_content_range(value: str) -> Tuple[int, int, int]:
    match = re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)", value.strip())
    if not match:
        raise DownloadError(f"invalid Content-Range header: {value!r}")
    return tuple(int(part) for part in match.groups())  # type: ignore[return-value]


def transfer_once(asset: Asset, final: Path, timeout: float) -> None:
    part, meta_path = meta_paths(final)
    final.parent.mkdir(parents=True, exist_ok=True)

    if part.exists() and not part.is_file():
        raise DownloadError(f"partial path is not a regular file: {part}")
    if meta_path.exists() and not meta_path.is_file():
        raise DownloadError(f"resume metadata path is not a regular file: {meta_path}")
    if part.is_file() and not meta_path.is_file():
        raise DownloadError(
            f"refusing to resume untracked partial file {part}; remove it or restore {meta_path.name}"
        )
    if meta_path.is_file() and not part.is_file():
        meta_path.unlink()

    offset = part.stat().st_size if part.is_file() else 0
    if offset > asset.size:
        raise IntegrityError(f"partial file is larger than expected: {part}")
    if offset == asset.size:
        ok, detail = verify_file(part, asset)
        if not ok:
            raise IntegrityError(f"complete partial file failed verification: {detail}")
        if final.exists():
            raise DownloadError(f"destination appeared while resuming: {final}")
        os.replace(part, final)
        if meta_path.exists():
            meta_path.unlink()
        return

    meta: Dict[str, object] = {}
    headers = {"User-Agent": USER_AGENT, "Accept-Encoding": "identity"}
    if offset:
        meta = load_resume_meta(meta_path)
        validate_resume_meta(meta, asset, meta_path)
        headers["Range"] = f"bytes={offset}-"
        if meta.get("etag"):
            headers["If-Range"] = str(meta["etag"])
        elif meta.get("last_modified"):
            headers["If-Range"] = str(meta["last_modified"])

    request = urllib.request.Request(asset.url, headers=headers, method="GET")
    with urllib.request.urlopen(request, timeout=timeout) as response:
        status = response_status(response)
        mode = "wb"
        if offset:
            if status == 206:
                content_range = response.headers.get("Content-Range", "")
                start, _end, total = parse_content_range(content_range)
                if start != offset or total != asset.size:
                    raise DownloadError(
                        f"unsafe range response for {asset.url}: start={start}, total={total}"
                    )
                remote_etag = response.headers.get("ETag", "")
                remote_modified = response.headers.get("Last-Modified", "")
                if meta.get("etag") and remote_etag and meta["etag"] != remote_etag:
                    raise DownloadError("ETag changed during resumed transfer")
                if (
                    not meta.get("etag")
                    and meta.get("last_modified")
                    and remote_modified
                    and meta["last_modified"] != remote_modified
                ):
                    raise DownloadError("Last-Modified changed during resumed transfer")
                mode = "ab"
            elif status == 200:
                remote_length = response.headers.get("Content-Length")
                if remote_length and int(remote_length) != asset.size:
                    raise DownloadError(
                        f"source size changed: manifest={asset.size}, remote={remote_length}"
                    )
                print(
                    f"WARN\t{asset.sample_id}\t{asset.kind}\tserver did not honor resume; restarting",
                    file=sys.stderr,
                )
                offset = 0
                mode = "wb"
            else:
                raise DownloadError(f"unexpected HTTP status {status} for resume")
        elif status != 200:
            raise DownloadError(f"unexpected HTTP status {status} for fresh download")

        if status == 200:
            remote_length = response.headers.get("Content-Length")
            if remote_length and int(remote_length) != asset.size:
                raise DownloadError(
                    f"source size changed: manifest={asset.size}, remote={remote_length}"
                )

        remote_meta = expected_meta(asset)
        remote_meta.update(
            {
                "etag": response.headers.get("ETag", ""),
                "last_modified": response.headers.get("Last-Modified", ""),
            }
        )
        write_resume_meta(meta_path, remote_meta)

        with part.open(mode) as handle:
            while True:
                block = response.read(CHUNK_SIZE)
                if not block:
                    break
                handle.write(block)
            handle.flush()
            os.fsync(handle.fileno())

    actual_size = part.stat().st_size
    if actual_size != asset.size:
        raise DownloadError(
            f"incomplete transfer for {asset.filename}: expected {asset.size}, found {actual_size}"
        )

    ok, detail = verify_file(part, asset)
    if not ok:
        raise IntegrityError(detail)
    if final.exists():
        raise DownloadError(f"destination appeared during transfer: {final}")
    os.replace(part, final)
    if meta_path.exists():
        meta_path.unlink()


def retryable_http(error: urllib.error.HTTPError) -> bool:
    return error.code in {408, 425, 429} or 500 <= error.code <= 599


def download_asset(asset: Asset, output_dir: Path, retries: int, timeout: float) -> None:
    final = asset_path(output_dir, asset)
    if final.exists():
        ok, detail = verify_file(final, asset)
        if ok:
            print(f"SKIP\t{asset.sample_id}\t{asset.kind}\t{detail}")
            return
        raise IntegrityError(
            f"refusing to overwrite existing destination {final}: {detail}"
        )

    for attempt in range(1, retries + 1):
        try:
            transfer_once(asset, final, timeout)
            ok, detail = verify_file(final, asset)
            if not ok:
                raise IntegrityError(detail)
            print(f"OK\t{asset.sample_id}\t{asset.kind}\t{detail}\t{final}")
            return
        except urllib.error.HTTPError as exc:
            if not retryable_http(exc) or attempt == retries:
                raise DownloadError(
                    f"HTTP {exc.code} for {asset.url}: {exc.reason}"
                ) from exc
            error: BaseException = exc
        except (
            urllib.error.URLError,
            http.client.IncompleteRead,
            http.client.RemoteDisconnected,
            socket.timeout,
            TimeoutError,
            ConnectionError,
            DownloadError,
        ) as exc:
            if isinstance(exc, IntegrityError):
                part, meta_path = meta_paths(final)
                if part.exists():
                    part.unlink()
                if meta_path.exists():
                    meta_path.unlink()
            if attempt == retries:
                raise DownloadError(
                    f"failed after {retries} attempt(s): {asset.url}: {exc}"
                ) from exc
            error = exc

        delay = min(2 ** (attempt - 1), 30)
        print(
            f"RETRY\t{asset.sample_id}\t{asset.kind}\tattempt={attempt}/{retries}"
            f"\tdelay={delay}s\t{error}",
            file=sys.stderr,
        )
        time.sleep(delay)


def download_assets(
    assets: Sequence[Asset],
    output_dir: Path,
    retries: int,
    timeout: float,
    jobs: int,
) -> None:
    """Download independent assets with bounded, file-level concurrency.

    Every worker owns a distinct final/partial/metadata path.  A failed asset
    therefore cannot corrupt another asset, and all already-running transfers
    are allowed to reach a resumable boundary before an aggregate error is
    reported.
    """
    if jobs == 1:
        for asset in assets:
            download_asset(asset, output_dir, retries, timeout)
        return

    failures: List[Tuple[Asset, BaseException]] = []
    with ThreadPoolExecutor(max_workers=jobs, thread_name_prefix="maizegdb") as pool:
        future_to_asset = {
            pool.submit(download_asset, asset, output_dir, retries, timeout): asset
            for asset in assets
        }
        for future in as_completed(future_to_asset):
            asset = future_to_asset[future]
            try:
                future.result()
            except Exception as exc:
                failures.append((asset, exc))
                print(
                    f"FAIL\t{asset.sample_id}\t{asset.kind}\t{exc}",
                    file=sys.stderr,
                )

    if failures:
        first_asset, first_error = failures[0]
        raise DownloadError(
            f"{len(failures)} asset(s) failed; first failure was "
            f"{first_asset.sample_id}/{first_asset.kind}: {first_error}"
        )


def positive_int(text: str) -> int:
    value = int(text)
    if value <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return value


def bounded_jobs(text: str) -> int:
    value = positive_int(text)
    if value > 16:
        raise argparse.ArgumentTypeError("must not exceed 16")
    return value


def positive_float(text: str) -> float:
    value = float(text)
    if value <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return value


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="List, verify, or explicitly download pinned MaizeGDB assets."
    )
    parser.add_argument(
        "--manifest", type=Path, default=DEFAULT_MANIFEST, help="pinned TSV manifest"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("maizegdb-downloads"),
        help="destination root (files are grouped by cohort/sample_id)",
    )
    parser.add_argument(
        "--cohort",
        action="append",
        choices=("NAM", "PanAnd"),
        default=[],
        help="select a cohort; repeat to select more than one",
    )
    parser.add_argument(
        "--sample",
        action="append",
        default=[],
        metavar="SAMPLE_ID",
        help="select an exact sample_id; repeat to select more than one",
    )
    parser.add_argument(
        "--asset",
        choices=("assembly", "annotation", "both"),
        default="both",
        help="asset type to process (default: both)",
    )
    parser.add_argument(
        "--retries",
        type=positive_int,
        default=4,
        help="maximum transfer attempts per file (default: 4)",
    )
    parser.add_argument(
        "--timeout",
        type=positive_float,
        default=60.0,
        help="per-request timeout in seconds (default: 60)",
    )
    parser.add_argument(
        "--jobs",
        type=bounded_jobs,
        default=1,
        help="concurrent file downloads, bounded to 1..16 (default: 1)",
    )
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--list", dest="action", action="store_const", const="list")
    action.add_argument(
        "--dry-run", dest="action", action="store_const", const="dry-run"
    )
    action.add_argument(
        "--verify-only", dest="action", action="store_const", const="verify"
    )
    action.add_argument(
        "--download",
        dest="action",
        action="store_const",
        const="download",
        help="explicitly authorize network downloads",
    )
    return parser


def run(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        records = load_manifest(args.manifest)
        selected_records = select_records(records, args.cohort, args.sample)
        assets = records_to_assets(selected_records, args.asset)
        if not assets:
            raise ManifestError("selection contains no downloadable assets")

        if args.action == "list":
            print_assets(assets, args.output_dir)
            return 0

        if args.action == "dry-run":
            print_assets(assets, args.output_dir)
            required, margin, free = disk_preflight(args.output_dir, assets, False)
            print(
                f"DRY-RUN: {len(assets)} asset(s); remaining={human_bytes(required)}; "
                f"margin={human_bytes(margin)}; free={human_bytes(free)}; no files downloaded",
                file=sys.stderr,
            )
            return 0

        if args.action == "verify":
            failures = 0
            for asset in assets:
                final = asset_path(args.output_dir, asset)
                ok, detail = verify_file(final, asset)
                label = "OK" if ok else "FAIL"
                print(f"{label}\t{asset.sample_id}\t{asset.kind}\t{detail}\t{final}")
                failures += 0 if ok else 1
            return 1 if failures else 0

        required, margin, free = disk_preflight(args.output_dir, assets, True)
        print(
            f"PREFLIGHT\tassets={len(assets)}\tremaining={required}\tmargin={margin}\tfree={free}",
            file=sys.stderr,
        )
        download_assets(
            assets, args.output_dir, args.retries, args.timeout, args.jobs
        )
        return 0
    except (ManifestError, DownloadError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


def main() -> None:
    raise SystemExit(run())


if __name__ == "__main__":
    main()
