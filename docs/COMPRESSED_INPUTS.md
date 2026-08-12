# Compressed genome and annotation inputs

AnchorWave accepts plain-text, gzip, and BGZF genome FASTA and GFF/GTF input
files. Compression is detected from the leading gzip magic bytes (`1f 8b`),
not from `.gz`, `.bgz`, or any other filename suffix. No command-line flag is
needed.

This applies to the genome and annotation inputs of `gff2seq`, `genoAli`, and
`proali`. SAM, anchor, MAF, and output compression are outside this input layer.

## Random-access behavior

AnchorWave's sequence extraction code uses byte offsets and `pread` to fetch
FASTA intervals. Ordinary gzip streams do not provide compatible random
access. A compressed FASTA is therefore decompressed once to a mode-0600 file
under `$TMPDIR` (or `/tmp` when `$TMPDIR` is unset). The materialized path is
cached and reused throughout the process. Compressed annotations use the same
validated materialization layer before the existing parser runs.

Managed files are removed on normal process exit. A forcibly killed process
can leave a file named `anchorwave-input-*`; such files can be removed after
confirming that no AnchorWave process is using them. `$TMPDIR` must have enough
space for the uncompressed inputs. BGZF is accepted as a concatenated gzip
member stream, but is likewise materialized rather than accessed through a
`.gzi` index.

## Examples

```bash
anchorwave gff2seq \
  -i reference.gff3.gz \
  -r reference.fa.gz \
  -o cds.fa

anchorwave genoAli \
  -i reference.gff3.bgz \
  -r reference.fa.gz \
  -s query.fa.gz \
  -as cds.fa -a query.sam -ar reference.sam \
  -n output.anchors
```

A gzip file without a gzip-like suffix is still decompressed. Conversely, a
plain file whose name ends in `.gz` is read as plain text. Corrupt gzip input
fails with the original input path in the error message and removes any partial
temporary file.

## Validation

Unit tests cover signature-based detection, misleading suffixes, concatenated
gzip members, corrupt streams, compressed FASTA random interval reads, and
compressed GFF parsing. A complete B73/Mo17 `genoAli -t 10` anchor-only run
using gzip FASTA/GFF inputs produced the same canonical anchor file as plain
input. That historical run predates the v1.3.0 correction of the `-fa3^2`
overflow and produced 109,413 lines:

```text
24fa6b9fac565f5860a0379ea3b51f065c9f0388724b5835fd92dc6ad6866ade
```

The hash validates compression equivalence for that build; it is not the
post-fix v1.3.0 anchor golden hash.
