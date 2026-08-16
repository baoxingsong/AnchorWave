# AnchorWave v1.3.0 release notes

AnchorWave v1.3.0 introduces resource-aware pairwise genome alignment while
preserving the two-piece-affine scoring objective of the exact algorithms.

## Alignment quality and selection

The production selector has two quality tiers:

1. exact dynamic programming: WFA-Singletrack high, WFA high, WFA medium,
   WFA low, KSW2-Singletrack, full KSW2, and score-certified KSW2;
2. fallback alignment: KSW2 banded and `alignSlidingWindow`.

BiWFA remains available only to the developer benchmark API and is not used by
the production `genoAli`, `proali`, or `ali` selector. Within the exact tier,
the scheduler chooses the predicted earliest safe completion under the current
thread and memory load. Temporary memory pressure does not lower alignment
quality. Within the fallback tier, the selector predicts which method will
produce the larger alignment score and normally executes only that method. At
the default scoring profile and `-w 100000`, the versioned model was fitted to
paired B73/Mo17 counterfactuals with chromosome-grouped validation and a
one-sided selective decision rule. Other settings use the deterministic
analytic quality rule rather than extrapolating the fitted model.

`-bt 0` is the default and selects exact-first mode. A positive value selects
balanced mode and limits every time quantile used to admit an exact candidate,
currently P50 and P90, as well as its cumulative runtime on that interval. The
historical `query_length * reference_length <= w^2` quality contract is the
explicit exception: such intervals remain in Tier 1 even in balanced mode.

## Resource controls

- `-w` retains its historical window role and normally defines a per-alignment
  algorithm budget of `w^2` bytes. The default 100,000 gives 10,000,000,000
  bytes (about 9.31 GiB). For `query_length * reference_length <= w^2`, all
  Tier-1 engines may request their larger predicted footprints to preserve the
  exact-alignment contract; process-wide `-M` admission still applies.
- `-t` is the maximum number of workers. Both anchor organization and sequence
  alignment use it in `genoAli`; `proali` also parallelizes quota DP, novel
  anchor searches, and sequence alignment.
- `-M` is the optional predictive process-wide memory envelope in GiB. When
  specified, it must be at least `w^2`, includes an internal safety reserve,
  and dynamically controls admission; omitting it disables the process-wide
  envelope. It is not an operating-system-enforced RSS cap: allocator peaks,
  caches, and sampling lag can exceed the prediction. Use a cgroup or
  batch-scheduler limit when a hard cap is required.

## User-visible behavior

- Plain, gzip, and BGZF FASTA and GFF/GTF inputs are detected automatically.
- The optional `-b` BED output reports every WFA execution mode as
  `WAVEFRONT`. KSW2-full, KSW2-Singletrack, and successfully
  score-certified KSW2 are all reported as `KSW2`.
  Fallback records retain the algorithm actually executed:
  `BANDED_KSW2` or `SLIDING_WINDOW`.
- `anchorwave --version` and `anchorwave -v` print `anchorwave v1.3.0`.
- The per-interval alignment trace remains an internal diagnostic interface
  and is intentionally absent from ordinary command help.
- The experimental multiple-genome graph command remains available for
  development and regression testing but is intentionally absent from the
  public command list and README usage guide.
- The `genoAli` and `proali` default for `-fa3` is now 50,000 (formerly
  100,000). This enables novel-anchor search in more gaps and can change anchor
  counts and interval fragmentation. Specify `-fa3 100000` to retain the
  former threshold.

## Correctness fixes

- All exact KSW2 paths now use generic 5-symbol scoring when a sequence
  contains ambiguous bases, and emitted CIGARs are re-scored with the same
  two-piece-affine objective used by WFA. This fixes the v1.2.6 behavior in
  which KSW2 treated the fifth symbol as a zero-cost wildcard.
- KSW2-Singletrack now supports both unbanded and native diagonal-banded
  traceback. Its zero-based leading-gap boundary bug was repaired, and the
  AVX2/AVX512 kernels now retain Singletrack's two difference tracks while
  sharing KSW2's wide score recurrence. ARM64 and Apple Silicon builds use the
  SSE2 interface through SSE2NEON/NEON. After 220,000 independent full-width
  executions across SSE, AVX2, AVX512, and NEON completed without a failure,
  the temporary cross-engine result validator and forced KSW2-full rescue were
  removed; ordinary allocation/deadline failure handling remains in the exact
  tier.
- Novel-anchor admission now evaluates `reference_length * query_length >
  fa3^2` with checked 64-bit arithmetic. This removes the historical signed
  32-bit overflow at the former default `-fa3 100000`; the resulting change in
  novel anchor count is intentional.
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
- The fallback-quality model was evaluated on 305 paired B73/Mo17 intervals
  where both approximate methods fit the production memory model. Nested
  chromosome-held-out validation selected banded KSW2 four times with no
  harmful banded choice and reduced total score regret from 502,447 for an
  always-sliding policy to 421,600. A replay of the final C++ coefficients had
  zero decision mismatch with the Python model; these are B73/Mo17 results,
  not a cross-species calibration claim.
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

- Final post-fix B73/Mo17 sequence-producing runs completed successfully on an
  18-core/36-thread Xeon W-2295:

  | command | settings | wall time | CPU time | peak RSS | exit |
  |---|---|---:|---:|---:|---:|
  | `proali -R 1 -Q 1` | former default `-fa3 100000`, `-t 34 -M 102 -bt 0` | 4:44:06 | 27.40 core-h | 106.61 GiB | 0 |
  | `proali -R 1 -Q 1` | new default `-fa3 50000`, `-t 34 -M 102 -bt 0` | 1:43:00 | 21.05 core-h | 97.92 GiB | 0 |
  | `genoAli` | former default `-fa3 100000`, `-t 34 -M 102 -bt 0` | 2:33:12 | 29.20 core-h | 95.69 GiB | 0 |
  | `genoAli` | new default `-fa3 50000`, `-t 34 -M 102 -bt 0` | 1:40:14 | 21.69 core-h | 97.30 GiB | 0 |

  Relative to the former v1.3.0 default on the same hardware, `-fa3 50000`
  reduced wall time by 63.75% for `proali` and 34.57% for `genoAli`. It added
  13,056 and 13,175 local anchors, respectively. Exact-tier reference coverage
  changed by -0.47 and -0.46 percentage points. Every common-coordinate
  fragment retained the same score and ungapped sequences (55,270 `proali`
  fragments and 56,361 `genoAli` fragments); the changed anchor boundaries
  mean complete fragment files are not expected to be byte-identical.

- Coordinate-keyed comparison found no invalid MAF blocks and identical total
  ungapped sequence content between v1.2.6 and v1.3.0. For common-coordinate
  Tier-1 intervals, all score differences were confined to N-containing
  intervals (463 for `proali`, 470 for `genoAli`), consistent with the generic
  KSW2 scoring fix above. Transcript anchors were unchanged; the corrected
  `-fa3` rule retained fewer local anchors, so fragment boundaries are not
  expected to be byte-identical to v1.2.6.

KSW2-Singletrack was adapted from `LorienLV/singletrack` commit
`7185195cb4049fd5290875ab4fc503384c4891dd`. See
`minimap2/KSW2_SINGLETRACK_NOTICE` and López-Villellas *et al.*,
*Bioinformatics* 42(5), btag183 (2026), DOI
`10.1093/bioinformatics/btag183`.
