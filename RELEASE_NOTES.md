# AnchorWave v1.3.0 release notes

AnchorWave v1.3.0 introduces resource-aware pairwise genome alignment while
preserving the two-piece-affine scoring objective of the exact algorithms.

## Alignment quality and selection

The production selector has two quality tiers:

1. exact dynamic programming: WFA-Singletrack high, WFA high, WFA medium,
   WFA low, full KSW2, and rescue-only score-certified KSW2;
2. fallback alignment: KSW2 banded and `alignSlidingWindow`.

BiWFA remains available only to the developer benchmark API and is not used by
the production `genoAli`, `proali`, or `ali` selector. Within the exact tier,
the scheduler chooses the predicted earliest safe completion under the current
thread and memory load. Temporary memory pressure does not lower alignment
quality. Within the fallback tier, the selector predicts which method will
produce the larger alignment score and normally executes only that method.

`-bt 0` is the default and selects exact-first mode. A positive value selects
balanced mode and limits every time quantile used to admit an exact candidate,
currently P50 and P90, as well as its cumulative runtime on that interval.

## Resource controls

- `-w` retains its historical window role and defines a hard per-alignment
  algorithm budget of `w^2` bytes. The default 100,000 gives 10,000,000,000
  bytes (about 9.31 GiB). It is not a sequence-length or score threshold for
  WFA.
- `-t` is the maximum number of workers. Both anchor organization and sequence
  alignment use it in `genoAli`; `proali` also parallelizes quota DP, novel
  anchor searches, and sequence alignment.
- `-M` is the predictive process-wide memory envelope in GiB. It must be at
  least `w^2`, includes an internal safety reserve, and dynamically controls
  admission. It never increases the per-alignment `w^2` ceiling. Use a cgroup
  or batch-scheduler limit when an operating-system-enforced hard cap is
  required.

## User-visible behavior

- Plain, gzip, and BGZF FASTA and GFF/GTF inputs are detected automatically.
- The optional `-b` BED output reports every WFA execution mode as
  `WAVEFRONT`. Full and successfully score-certified KSW2 are `MINIMAP2`.
  Fallback records retain the algorithm actually executed:
  `BANDED_MINIMAP2` or `SLIDING_WINDOW`.
- `anchorwave --version` and `anchorwave -v` print `anchorwave v1.3.0`.
- The per-interval alignment trace remains an internal diagnostic interface
  and is intentionally absent from ordinary command help.
- The experimental multiple-genome graph command remains available for
  development and regression testing but is intentionally absent from the
  public command list and README usage guide.

## Correctness fixes

- Novel-anchor admission now evaluates `reference_length * query_length >
  fa3^2` with checked 64-bit arithmetic. This removes the historical signed
  32-bit overflow at the default `-fa3 100000`; the resulting change in novel
  anchor count is intentional.
- Empty quota-DP candidate rows are rejected before accessing
  `AlignmentMatch`, removing the libc++ container-overflow reported by ASan.
- The bundled minimap2 radix-sort wrapper now handles empty ranges without
  forming a pointer outside a null array, removing the corresponding UBSan
  report.

## Validation summary

- Release core tests pass for native AVX-512 and explicit AVX2, SSE4.1, and
  SSE2 builds.
- On the real 4-Mb B73/Mo17 input, AVX2, SSE4.1, and SSE2 produced
  byte-identical anchors, method BED, whole MAF, and fragment MAF for both
  `proali` and `genoAli`.
- A sequence-producing B73/Mo17 `genoAli -t 34 -M 102 -bt 0` selector run made
  before the final `-fa3` overflow fix completed in 1:32:18 with 102,139,780
  KiB peak RSS and kept 96,727 of 97,549 evaluated non-filling intervals in
  the exact tier. It remains a resource measurement, not a post-fix anchor
  golden result.
- Final post-fix, full-genome anchor-only validation completed successfully:

  | command | threads | wall time | peak RSS | anchor lines | canonical body SHA-256 |
  |---|---:|---:|---:|---:|---|
  | `proali -R 1 -Q 1` | 34 | 1:17.31 | 2,247,984 KiB | 66,744 | `48ed1fc723800644fa81b4667b0f686e98a4a7d275f30ee74a835cfe05c86523` |
  | `proali -R 1 -Q 1` | 1 | 10:03.39 | 527,236 KiB | 66,744 | same |
  | `genoAli` | 34 | 1:11.49 | 2,276,688 KiB | 67,842 | `a12dd94fba69babe3d8da539841c221bb0d6a88fddaf9e11f29f21256a102d14` |
  | `genoAli` | 1 | 9:35.10 | 527,812 KiB | 67,842 | same |

  For each command, the one-thread and 34-thread bodies are byte-identical;
  only the provenance command line in the first row differs.
