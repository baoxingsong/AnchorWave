# Pinned MaizeGDB assets for AnchorWave trio testing

`maizegdb_assets.tsv` is a machine-readable, version-pinned inventory of the
maize and teosinte data selected for the AnchorWave trio workflow. It contains
exactly 35 rows:

- 26 NAM genomes: B73 v5 and the 25 NAM founder inbreds;
- all nine `Zea` assemblies currently listed in the MaizeGDB PanAnd collection.

The unrelated B73-Ab10 assembly
`Zm-B73_AB10-REFERENCE-NAM-1.0` is deliberately excluded. The PanAnd project
also contains non-`Zea` Andropogoneae genomes; those are outside this manifest.

## Provenance and pinning

The inventory was audited on 2026-08-06 against these official resources:

- [MaizeGDB NAM founder download index](https://download.maizegdb.org/Genomes/NAM_Founders/)
- [MaizeGDB PanAnd download index](https://download.maizegdb.org/Genomes/PanAnd/)
- [NCBI BioProject PRJEB31061](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB31061), containing the 25 non-B73 NAM assemblies
- [NCBI BioProject PRJEB32225](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB32225), containing B73 v5

MD5 values come from the `ASSEMBLY_DIR.md5` files in the corresponding
MaizeGDB directories. Exact byte sizes were checked from the official HTTP
object metadata. Assembly and annotation URLs intentionally point into the
same directory so that coordinate-compatible releases are not silently mixed.

The selected official annotations are the current `*ab.1` gene-model sets, not
the older `*aa.1` sets and not Helixer predictions. B73 is the one packaging
exception: its selected `Zm00001eb.1.gff3` is uncompressed. The other selected
GFF3 files are gzip-compressed.

Two exceptions are explicit in the TSV:

- `Zx-TIL25` has an official `Zx00003ab.1` GFF3, but the current MaizeGDB MD5
  manifest only lists the older `aa.1` annotation. Its `annotation_md5` is
  `MISSING_IN_SOURCE_MANIFEST`, its exact byte size is retained, and its status
  is `official_checksum_missing`. The downloader can verify its size but must
  warn that source-level checksum verification is unavailable.
- `Zl-RIL003` is an assembly-only `Zea luxurians` release. Its MaizeGDB README
  explicitly states that no official annotation exists. The available Helixer
  prediction is therefore not presented as an official gene annotation.

NCBI currently labels `GCA_902373975.1` as
`Zm-CML52-REFERENCE-NAM-2.1`, whereas the current MaizeGDB package and its
coordinate-matched annotation remain under directory
`Zm-CML52-REFERENCE-NAM-1.0`. The workflow pins the MaizeGDB package name and
records the stable GCA accession; do not derive a replacement URL from NCBI's
display name.

Pinned payload totals are approximately:

| Cohort | Assemblies | Selected annotations | Combined |
|---|---:|---:|---:|
| NAM | 17,267,418,815 bytes (16.082 GiB) | 455,454,401 bytes (0.424 GiB) | 16.506 GiB |
| PanAnd `Zea` | 7,106,420,544 bytes (6.618 GiB) | 116,849,149 bytes (0.109 GiB) | 6.727 GiB |

These are compressed download sizes except for the B73 GFF3 noted above.
Substantially more working storage is required after decompression and for
pairwise AnchorWave output.

## TSV fields

Each row contains:

- `sample_id`, `cohort`, `taxon`: stable workflow identifiers and biological grouping;
- `assembly_dir`: the exact MaizeGDB release directory;
- `assembly_url`, `assembly_filename`, `assembly_md5`, `assembly_size_bytes`;
- `annotation_url`, `annotation_filename`, `annotation_md5`, `annotation_size_bytes`;
- `annotation_status`: `official`, `official_checksum_missing`, or `assembly_only`;
- `model`: official gene-model release identifier when present;
- `accession`: NCBI GCA/GCF for NAM, otherwise the versioned MaizeGDB assembly identifier;
- `notice_url`: the relevant catalog or PanAnd data-use notice.

Empty annotation fields are intentional only for `assembly_only` rows.

## Downloader

`../download_maizegdb.py` uses only the Python standard library. It requires an
explicit mutually exclusive action, so a typo or import cannot initiate these
large downloads.

List a selection without network or filesystem changes:

```bash
python3 workflow/download_maizegdb.py \
  --cohort NAM --sample B73 --asset both --list
```

Check the plan and disk-space requirement without downloading:

```bash
python3 workflow/download_maizegdb.py \
  --cohort PanAnd --asset assembly \
  --output-dir /data/anchorwave/maizegdb \
  --dry-run
```

Download only after explicitly supplying `--download`:

```bash
python3 workflow/download_maizegdb.py \
  --cohort NAM --sample B73 --sample B97 --asset both \
  --output-dir /data/anchorwave/maizegdb \
  --download
```

Verify already downloaded files without network access:

```bash
python3 workflow/download_maizegdb.py \
  --cohort NAM --sample B73 --asset both \
  --output-dir /data/anchorwave/maizegdb \
  --verify-only
```

Files are organized as `OUTPUT_DIR/COHORT/SAMPLE_ID/FILENAME`. Interrupted
transfers remain as `FILENAME.part` plus a small resume-metadata file. Resume
uses HTTP Range and If-Range, validates the returned Content-Range, then checks
the pinned byte count and MD5 before an atomic rename. If a server version
changes, the pinned size/checksum causes a hard failure instead of silently
accepting new content. Existing invalid final files are never overwritten.

## Data-use terms

The [USDA data.gov entry for MaizeGDB](https://catalog.data.gov/dataset/maizegdb)
marks the database as public domain. That catalog-level designation must not be
used to erase producer-specific conditions attached to individual files.

PanAnd packages are distributed under the
[Toronto Agreement notice shipped by MaizeGDB](https://download.maizegdb.org/Genomes/PanAnd/Zd-Momo-REFERENCE-PanAnd-1.0/TORONTO_AGREEMENT),
which permits research use and originally imposed a prepublication moratorium.
The consortium study has since been published as Stitzer et al.,
“[Extensive genome evolution distinguishes maize within a stable tribe of grasses](https://pmc.ncbi.nlm.nih.gov/articles/PMC11785232/)”.
Because MaizeGDB still serves the notice for these releases, retain it with the
data, cite the producers, and confirm current publication expectations before
publishing benchmark results.

For the NAM assemblies and annotations, cite Hufford et al.,
“[De novo assembly, annotation, and comparative analysis of 26 diverse maize genomes](https://pmc.ncbi.nlm.nih.gov/articles/PMC8733867/)”,
along with MaizeGDB and the exact assembly/gene-model identifiers in the TSV.

## Phylogenetic-tree warning

A randomly generated 26-tip tree is acceptable only as a deterministic parser,
scheduler, and stress-test fixture. It is not a biological estimate and must
not be used to assess ancestral sequence, ancestral copy number, branch-specific
events, or ancestral karyotype accuracy. Biological validation requires a
taxonomically and phylogenetically defensible tree with documented branch
lengths. Store the random seed and generated Newick whenever a random fixture
is used so test failures are reproducible.
