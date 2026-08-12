# Reproducible MaizeGDB TrioAnchorGraph workflow

This directory separates three expensive or scientifically distinct actions:

1. `download_maizegdb.py` downloads pinned MaizeGDB assets, but only with an
   explicit `--download` action;
2. `trioanchorgraph_workflow.py` prepares local inputs and creates pairwise
   AnchorWave evidence, but only with `--execute`;
3. the generated `run_trioGraphAli.sh` command consumes the pairwise MAFs and
   builds the authoritative multi-genome graph. The workflow emits this command
   but does not run it automatically.

No workflow command performs an implicit network request.

## Pairwise schedule

For ingroups `I1`, `I2` and outgroups `O1..Ok`, the runner creates exactly
`1 + 2k` pairwise MAFs:

```text
I1 (reference) -- I2
I1 (reference) -- Oj    for every outgroup j
I2 (reference) -- Oj    for every outgroup j
```

Both ingroups therefore require coordinate-matched annotations. Outgroups are
queries and may be assembly-only; this supports the PanAnd `Zl-RIL003` package,
which has no official annotation. For each ingroup reference the workflow:

1. decompresses/copies its assembly and GFF3 with an atomic rename;
2. runs `anchorwave gff2seq` once;
3. maps those anchor CDSs to the reference and every required query with the
   documented minimap2 splice settings
   `-x splice -k 12 -a -p 0.4 -N 20`;
4. reuses the reference SAM and CDS to run each `anchorwave proali` edge.

`--min-exon` is passed identically to `gff2seq` and `proali`. The proali
coverage parameters default to `-R 1 -Q 1`; change them only when the expected
copy-depth relationship supports a different value. These alignment-depth
limits are independent from copy constraints later supplied to `trioGraphAli`.

## Download the pinned data

The exact 26 NAM and nine PanAnd `Zea` releases, checksums, sizes, provenance,
and data-use notes are documented in [data/README.md](data/README.md). First
inspect the 23+ GiB compressed plan:

```bash
python3 workflow/download_maizegdb.py \
  --asset both --output-dir /data/anchorwave/maizegdb --dry-run
```

Then explicitly download all pinned rows:

```bash
python3 workflow/download_maizegdb.py \
  --asset both --output-dir /data/anchorwave/maizegdb --download
```

The downloader validates pinned sizes and available MD5 values. Allow ample
additional space for decompressed FASTA, SAM, MAF, graph, and temporary files.

## Plan before execution

`--dry-run` prints the complete JSON plan and final `trioGraphAli` argv without
creating the run directory or checking for downloaded payloads:

```bash
python3 workflow/trioanchorgraph_workflow.py \
  --dry-run \
  --data-root /data/anchorwave/maizegdb \
  --run-dir /work/anchorwave/B73_B97_Zv \
  --ingroup1 B73 --ingroup2 B97 \
  --outgroup Zv-TIL01 \
  --random-tree-seed 20260806 \
  --ancestor-node NAM_ANCESTOR \
  --minimap2-threads 16 \
  --anchorwave-threads 6
```

`--plan` writes the same plan, manifests, tree, and executable final command,
but performs no decompression or alignment:

```bash
python3 workflow/trioanchorgraph_workflow.py \
  --plan \
  --data-root /data/anchorwave/maizegdb \
  --run-dir /work/anchorwave/B73_B97_Zv \
  --ingroup1 B73 --ingroup2 B97 \
  --outgroup Zv-TIL01 --outgroup Zl-RIL003 \
  --random-tree-seed 20260806 \
  --ancestor-node NAM_ANCESTOR \
  --minimap2-threads 16 \
  --anchorwave-threads 6
```

The first outgroup is written as `primary_outgroup`; additional outgroups are
written as `outgroup`. To use a biological phylogeny, replace
`--random-tree-seed` with `--species-tree FILE` and name the ingroup ancestral
node using `--ancestor-node`.

## Execute and resume safely

Only `--execute` runs tools:

```bash
python3 workflow/trioanchorgraph_workflow.py \
  --execute \
  --data-root /data/anchorwave/maizegdb \
  --run-dir /work/anchorwave/B73_B97_Zv \
  --ingroup1 B73 --ingroup2 B97 \
  --outgroup Zv-TIL01 \
  --species-tree /work/trees/zea_rooted.nwk \
  --ancestor-node NAM_ANCESTOR \
  --anchorwave /path/to/anchorwave \
  --minimap2 /path/to/minimap2 \
  --minimap2-threads 16 \
  --anchorwave-threads 6
```

The execution uses subprocess argv arrays, never shell interpolation. Every
large output is first written to a same-directory temporary file and renamed
only after successful completion. Per-step stdout/stderr logs, events, input
signatures, and output fingerprints are stored under `provenance/`.

After interruption, add `--resume`. A result is reused only if its receipt,
step definition, input sizes/timestamps, and output fingerprints still match.
An existing output without a matching receipt is recomputed atomically. Without
`--resume`, pre-existing step outputs cause a hard error.

Source byte counts and pinned MD5 values are checked before preparation. The
exceptional PanAnd checksum status remains explicit in the catalog. The
`--skip-source-md5` escape hatch retains byte-count checking but is discouraged.
For intentionally modified local input files, `--skip-source-integrity` checks
only that each selected source file exists and skips both catalog byte counts
and MD5. This weaker mode is recorded in `workflow.plan.json`.

Execution is deliberately serial: one workflow step and therefore one
pairwise `proali` job runs at a time. Thread-level parallelism remains available
inside each tool through independent `--minimap2-threads` and
`--anchorwave-threads` values. The legacy `--threads` value is only a fallback
when a tool-specific value is omitted.

AnchorWave `proali` uses stdout for a terminal progress display containing many
backspace characters. On whole-genome jobs that stream can otherwise create
multi-gigabyte logs in minutes. The workflow discards only this progress stream,
records a note in the per-step stdout log, captures stderr normally, and validates
the three declared `-n`, `-o`, and `-f` result files before marking the step done.

After all pairwise jobs finish, inspect and run:

```bash
/work/anchorwave/B73_B97_Zv/manifests/run_trioGraphAli.sh
```

The generated taxon and pairwise manifests use absolute paths and include the
`proali` score/depth profile. Pass graph-level copy information with
`--copy-relations FILE` and/or repeatable `--copy-number TAXON=N` arguments to
the workflow; they are forwarded to the emitted command. Add
`--copy-mode strict` when every retained relationship must come from explicit
hard copy constraints rather than the default constrained inference mode.
An existing versioned macro-synteny projection TSV can be forwarded independently
with `--block-projections FILE`; it is used for ancestral adjacency/karyotype
inference and does not alter the sequence-alignment graph.

The emitted `trioGraphAli` argv always contains `--validate-input-paths`.
FASTA-to-MAF residue and source-length validation remains enabled through the
current `trioGraphAli` default; the workflow deliberately has no option that
emits `--skip-fasta-validation`. Thus an invalid FASTA/MAF combination fails
before graph construction rather than silently entering the evidence graph.

## Random-tree fixture warning

`random_newick.py` is deterministic for a fixed sample list, seed, and branch
range. It records the seed, input-manifest hash, Newick hash, algorithm version,
and a prominent warning in JSON metadata. A random tree is only a parser,
scheduler, and stress-test fixture. It is invalid evidence for ancestral
sequence, copy number, rearrangement, branch-specific events, or ancestral
karyotype.

Generate the 26-tip NAM software fixture explicitly:

```bash
python3 workflow/random_newick.py \
  --cohort NAM --seed 20260806 \
  --output /work/trees/NAM26.random.stress-test.nwk
```

## Tests

The unit suite uses tiny gzipped assemblies plus fake `anchorwave` and
`minimap2` executables. It exercises all three primary-trio pairwise jobs and a
receipt-validated resume without downloading data:

```bash
python3 -m unittest discover -s workflow/tests -v
```
