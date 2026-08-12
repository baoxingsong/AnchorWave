# Compact B73/Mo17 selector calibration

## Scope and safety

This workflow uses only the existing B73 v4 and Mo17 CAU data. It never
downloads data and does not start a whole-genome alignment automatically. Its
purpose is to collect repeatable selector evidence before another expensive
whole-genome run.

There are two complementary measurements:

1. `anchorwave ali` replays extracted inter-anchor sequence pairs. This is
   inexpensive and reports per-pair wall time, peak process RSS, output hash,
   and sequence-reconstruction validity.
2. bounded `genoAli` runs the first 4 Mb of chromosome 1 and exercises the real
   process-wide `-t/-M` scheduler. It reports whole-process wall time/RSS, the
   actual method distribution, and anchor/MAF hashes.

`ali` has no `-t` or `-M` option and does not install the production memory
scheduler. The replay runner accepts these options only to limit the number of
independent `ali` processes conservatively. Its peak RSS field is the maximum
of one `ali` process, not aggregate concurrent RSS. Any conclusion about
dynamic memory admission or tail draining must come from bounded `genoAli`,
not replay. `-M` is an in-process predictive envelope; an operating-system or
batch-scheduler limit is required for hard kill-at-cap enforcement.

Do not request a memory scenario larger than the physical RAM of the host.
In particular, `-M 500` and `-M 1000` must be run only on matching large-memory
nodes even though they are accepted by the scripts.

## Existing B73/Mo17 evidence

With `DATA_ROOT` set to
`../AnchorWave-wfa2-benchmark-data/installation` relative to the repository:

| purpose | path | status |
|---|---|---|
| full prepared B73/Mo17 inputs | `$DATA_ROOT/work/` | FASTA, GFF, CDS and both SAM files present |
| compact input | `$DATA_ROOT/bounded-chr1-4000000/work/` | first 4 Mb of chromosome 1, complete |
| compact source trace | `$DATA_ROOT/bounded-chr1-4000000/results/selector.methods.bed` | about 244 non-filling intervals |
| broader source trace | `$DATA_ROOT/biwfa-time-limit-20260809/full-genome-calibrated/t12.methods.bed` | 25,247 rows across chromosomes 1–10 |
| compact regression outputs | `$DATA_ROOT/bounded-chr1-4000000/results/` | anchors, MAF, fragment MAF, BED and GNU time |
| stopped Singletrack run | `$DATA_ROOT/full-b73-mo17-singletrack-bt30-20260809/` | partial only; exit 137; never use as ground truth |

The source method in a trace manifest is only a stratification label. It is not
a forced algorithm and is not necessarily the algorithm chosen by a later
binary during `ali` replay.

## Build once

Use a Release binary. Existing unit/sanitizer builds remain separate from
performance results.

```bash
cmake -S . -B build-selector-benchmark \
  -DCMAKE_BUILD_TYPE=Release \
  -DANCHORWAVE_BUILD_TESTS=OFF \
  -DANCHORWAVE_BUILD_BENCHMARKS=ON
cmake --build build-selector-benchmark \
  --target anchorwave anchorwave_exact_candidate_benchmark -j 8
```

On macOS the replay runner expects GNU time as `gtime`.

## 100-Mb paired-candidate corpus

The larger calibration uses every inter-anchor gap whose B73 and Mo17
coordinates both fall within the first 100 Mb of chromosome 1. It does **not**
ask an exact aligner to align two monolithic 100-Mb strings. Instead, it
replays the real 100-Mb genome-alignment workload as independent gaps, which
is the unit seen by the production selector. In the 2026-08-10 anchor set this
produces 2,465 pairs (about 65.3 Mb of B73 and 64.5 Mb of Mo17 sequence), with
individual lengths from 1 bp through 464,378 bp.

Convert the completed anchor organization and extract the immutable corpus:

```bash
DATA_ROOT=../AnchorWave-wfa2-benchmark-data/installation
CAL_ROOT="$DATA_ROOT/selector-100m-20260810"
WORK_100M="$DATA_ROOT/bounded-chr1-100000000/work"

python3 workflow/anchor_gaps_to_bed.py \
  --anchors "$DATA_ROOT/full-b73-mo17-dynamic-selector-20260810/full.anchors" \
  --output "$CAL_ROOT/corpus/chr1-100m.anchor-gaps.bed" \
  --reference-chromosome 1 --query-chromosome 1 \
  --reference-end 100000000 --query-end 100000000

python3 workflow/b73_mo17_trace.py prepare \
  --methods "$CAL_ROOT/corpus/chr1-100m.anchor-gaps.bed" \
  --reference "$WORK_100M/B73.chr1.100000000.fa" \
  --query "$WORK_100M/Mo17.chr1.100000000.fa" \
  --output "$CAL_ROOT/corpus/all-gaps" \
  --sample-size 0 --max-sequence-bases 500000 --seed 20260810
```

Run a 192-pair quantile-stratified counterfactual benchmark. Each invocation
contains exactly one pair and one algorithm; the runner never overlaps them.
High WFA and elastic KSW2-full may use up to 64 GiB in isolation, while the
succinct WFA modes retain the production `w^2 = 10 GB`-class budget. The
21-GiB difference from the 85-GiB process setting represents non-alignment
state and scheduler reserve.

```bash
python3 workflow/benchmark_exact_candidates.py \
  --trace "$CAL_ROOT/corpus/all-gaps" \
  --benchmark "$PWD/build-selector-benchmark/anchorwave_exact_candidate_benchmark" \
  --output "$CAL_ROOT/exact-paired-round1" \
  --sample-size 192 \
  --algorithms ksw2-score-certified,ksw2-full,singletrack,wfa-high,wfa-medium,wfa-low,biwfa \
  --memory-gib 64 \
  --algorithm-memory-gib wfa-medium=10 \
  --algorithm-memory-gib wfa-low=10 \
  --algorithm-memory-gib biwfa=10 \
  --timeout-minutes 10 \
  --algorithm-timeout-minutes wfa-medium=30 \
  --algorithm-timeout-minutes wfa-low=30 \
  --algorithm-timeout-minutes biwfa=30 \
  --wfa-threads 1 --biwfa-leaf-score 0

python3 workflow/summarize_exact_candidate_calibration.py \
  --results "$CAL_ROOT/exact-paired-round1/results.tsv" \
  --output "$CAL_ROOT/exact-paired-round1/calibration.json"
```

Timeouts and analytic KSW2 memory rejections are right-censored records; they
must not be treated as measured runtimes. `results.tsv` includes the current
selector's P50/P90 and memory predictions beside actual internal time, process
RSS, WFA allocator peak, exact score, reconstruction status, adaptive BiWFA
leaf statistics, and actual OpenMP team size. `score_inconsistencies.tsv` must
contain only its header before a parameter set is accepted. The receipt tree
makes a stopped run resumable with the identical command plus `--resume`.
For `ksw2-score-certified`, the result also records the global score, planned
initial/maximum band, final certified band and number of traceback attempts.
`alignment_path_differences.tsv` separately records equal-score tie-breaking
differences; these are not score errors, but they matter when reproducible MAF
paths are required across hardware configurations.
The calibration summary reports P50/P90 time residuals, empirical P90
coverage, allocator/process-memory residuals, length/ratio strata, and the
paired standalone fastest-engine regret. It intentionally does not rewrite C++
constants automatically.

Use a small hard-task subset for isolated WFA 1/2/4/8-thread scaling. Do not
enable inner threads in production merely because a microbenchmark scales:
the genome scheduler must first prove that enough outer workers are idle and
that the total active-thread cap remains respected.

## Round 1: exhaustive 4-Mb trace

Extract every eligible non-filling interval. The output is immutable: choose a
new directory if parameters or source files change.

```bash
DATA_ROOT=../AnchorWave-wfa2-benchmark-data/installation
TRACE_ROOT="$DATA_ROOT/selector-calibration/round1-4mb"

python3 workflow/b73_mo17_trace.py prepare \
  --methods "$DATA_ROOT/bounded-chr1-4000000/results/selector.methods.bed" \
  --reference "$DATA_ROOT/bounded-chr1-4000000/work/B73.chr1.4000000.fa" \
  --query "$DATA_ROOT/bounded-chr1-4000000/work/Mo17.chr1.4000000.fa" \
  --output "$TRACE_ROOT/trace" \
  --sample-size 0 --max-sequence-bases 500000 --seed 20260809
```

First run isolated replay to measure segment-level behavior without scheduler
contention:

```bash
workflow/benchmark_b73_mo17_trace.sh replay \
  --trace-dir "$TRACE_ROOT/trace" \
  --output "$TRACE_ROOT/replay-current-t1" \
  --anchorwave "$PWD/build-selector-benchmark/anchorwave" \
  -t 1 -M 0 -w 100000 -bt 30 --timeout-seconds 300
```

Then run bounded `genoAli` for the relevant hardware scenarios. This is the
authoritative dynamic-scheduler measurement. The first 4-Mb run normally takes
minutes, but the script does not start it unless `bounded` is invoked.

```bash
workflow/benchmark_b73_mo17_trace.sh bounded \
  --output "$TRACE_ROOT/bounded-t12-m80" \
  --anchorwave "$PWD/build-selector-benchmark/anchorwave" \
  -t 12 -M 80 -w 100000 -bt 30
```

For the six target configurations, use new output directories:

```text
t12-m80    -t 12  -M 80
t36-m120   -t 36  -M 120
t96-m250   -t 96  -M 250
t128-m250  -t 128 -M 250
t128-m500  -t 128 -M 500
t128-m1000 -t 128 -M 1000
```

Round 1 is used to correct gross model errors and choose a small number of
candidate settings. Inspect:

- `summary.json`: wall/RSS, scenario, hashes and method distribution;
- `tasks.tsv`: per-gap replay wall/RSS/status/hash;
- `run.stderr.txt`: selector and scheduler summaries from bounded `genoAli`;
- exact-tier interval and base fractions, memory failures, deferrals, waits,
  preferred reservations, and fallback counts.

Do not tune against only total wall time. Reject any setting that reduces
exact-tier coverage or exceeds `-M`.

## Round 2: stratified genome-wide trace

After freezing the Round-1 candidate settings, extract a compact deterministic
sample from the completed B73/Mo17 method BED. Sampling round-robins over old
method, length bin, length-ratio bin and strand, so rare long or imbalanced
tasks are not hidden by the many small WFA intervals.

```bash
TRACE_ROOT="$DATA_ROOT/selector-calibration/round2-genomewide"

python3 workflow/b73_mo17_trace.py prepare \
  --methods "$DATA_ROOT/biwfa-time-limit-20260809/full-genome-calibrated/t12.methods.bed" \
  --reference "$DATA_ROOT/work/Zea_mays.AGPv4.dna.toplevel.fa" \
  --query "$DATA_ROOT/work/Zm-Mo17-REFERENCE-CAU-1.0.fa" \
  --output "$TRACE_ROOT/trace" \
  --sample-size 192 --max-sequence-bases 500000 --seed 20260809
```

Run the same immutable trace with the baseline and candidate binaries. The
`-t/-M` replay scenario is conservative external process admission, not the
production scheduler; use it to detect sequence-dependent regressions and use
the Round-1 bounded runs to validate scheduling.

```bash
workflow/benchmark_b73_mo17_trace.sh replay \
  --trace-dir "$TRACE_ROOT/trace" \
  --output "$TRACE_ROOT/replay-baseline-t12-m80" \
  --anchorwave /path/to/baseline/anchorwave \
  -t 12 -M 80 -w 100000 -bt 30 --timeout-seconds 300

workflow/benchmark_b73_mo17_trace.sh replay \
  --trace-dir "$TRACE_ROOT/trace" \
  --output "$TRACE_ROOT/replay-candidate-t12-m80" \
  --anchorwave /path/to/candidate/anchorwave \
  -t 12 -M 80 -w 100000 -bt 30 --timeout-seconds 300

workflow/benchmark_b73_mo17_trace.sh compare \
  --left "$TRACE_ROOT/replay-baseline-t12-m80" \
  --right "$TRACE_ROOT/replay-candidate-t12-m80" \
  --output "$TRACE_ROOT/baseline-vs-candidate-t12-m80.json"
```

The comparison reports common-task output hashes, wall-time change and
single-process RSS change. For bounded pipeline directories it compares
semantic anchor hashes, MAF hashes, actual method distributions and
whole-process RSS.

Round 2 should change only one parameter family at a time:

1. score/profile confidence and candidate time coefficients;
2. memory coefficients and safety reserve;
3. wait/backfill/tail/aging policy;
4. `-bt` only after the first three are stable.

Keep a candidate only when it satisfies all of these conditions on the held
trace: no new invalid/timeout tasks, no loss of exact-tier bases, no `-M`
violation, no unexpected score/output reconstruction failure, and a repeatable
wall-time reduction. Repeat noisy cells; do not launch a whole-genome run until
both rounds are frozen.

## Completed calibration on 2026-08-09

The implemented selector was calibrated and then retested with the bounded
B73/Mo17 chromosome-1 input, `-t 12 -M 85 -w 100000 -bt 30 -IV`. No
whole-genome alignment was started.

| measurement | Round 1 | tuned Round 2 | change |
|---|---:|---:|---:|
| wall time | 116.28 s | 90.00 s | -22.60% (1.292x speedup) |
| sampled peak RSS | 45.97 GiB | 19.25 GiB | -58.12% |
| exact intervals | 244/244 | 244/244 | unchanged |
| exact runtime memory failures | 4 | 0 | -4 |
| banded/sliding intervals | 0/0 | 0/0 | unchanged |

Round 2 selected WFA-Singletrack high for 149 intervals, full KSW2 for 94,
and BiWFA for one. The final BiWFA interval took 82.41 s and was the critical
path, so most remaining wall time is tail latency rather than aggregate work.
Semantic anchor rows were identical between rounds. The final fragment MAF
contains 253 blocks, reconstructs all 8,000,000 ungapped input bases, and all
253 declared scores equal an independent CIGAR/alignment score recomputation.

Calibration exposed an objective mismatch for ambiguous bases: KSW2's
five-symbol fast matrix had treated its fifth symbol as a zero-cost wildcard,
whereas WFA charged unequal ambiguous symbols as mismatches. The production
KSW2 paths now enable the generic 5x5 scoring kernel only when an ambiguous
symbol is present and use the same strict substitution objective as WFA.
Ordinary A/C/G/T intervals retain the optimized KSW2 kernel. Five old N-rich
blocks changed score after this correction; their coordinates and ungapped
sequences are unchanged, and every corrected score is better under the common
strict objective.

The Release suite passed 8/8 tests. The complete AddressSanitizer and
UndefinedBehaviorSanitizer suite also passed 8/8. Unit resource-plan cases
cover all six requested `(-t, -M)` inputs, but only `(12, 85)` was physically
benchmarked on this host.

After Round 2, three edge-condition fixes were added: persistent RSS remains a
temporary admission state, the deferred `-bt` gate activates after both high
WFA attempts fail, and inter-anchor results use a bounded per-result rolling
window instead of a whole-batch barrier. These paths passed unit and sanitizer
tests. Round 2 itself reported zero admission deferrals and zero runtime
failures, so its measurements remain the parameter-calibration reference.

The post-audit source was then confirmed once on the same 4-Mb input without
starting a new tuning round:

| post-audit confirmation | value |
|---|---:|
| wall time | 86.46 s |
| GNU-time peak RSS | 25.58 GiB |
| exact / banded / sliding intervals | 244 / 0 / 0 |
| exact runtime failures | 0 |
| Singletrack / full KSW2 / BiWFA | 149 / 94 / 1 |
| maximum concurrent reservations | 12 |
| admission waits / deferrals | 0 / 0 |

Both MAF files are byte-identical to tuned Round 2, and the semantic anchor
hash is identical. Method records are also identical; only their explanatory
header gained the missing `WAVEFRONT_SINGLETRACK` definition. The 3.54-second
wall reduction is one repeat and should not be treated as a statistically
calibrated speedup. Peak RSS was 6.33 GiB above Round 2 but remained only 30.1%
of the requested 85-GiB envelope; the scheduler's retained-RSS accounting made
the projected peak conservative rather than misclassifying this rise as
structural infeasibility.

The immutable result directories are:

```text
$DATA_ROOT/selector-global-v1-20260809/round1
$DATA_ROOT/selector-global-v1-20260809/round2-strict
$DATA_ROOT/selector-global-v1-20260810/post-audit-strict
```

For a one-pair, one-engine counterfactual outside the batch runner, configure
with `-DANCHORWAVE_BUILD_BENCHMARKS=ON` and run, for example:

```bash
build-selector-benchmark/anchorwave_exact_candidate_benchmark \
  --reference /path/to/gap.ref.fa \
  --query /path/to/gap.query.fa \
  --algorithm singletrack \
  --memory-bytes 68719476736 \
  --wfa-threads 1 --biwfa-leaf-score 0
```

Valid candidate names are `ksw2-score-certified`, `ksw2-full`, `singletrack`,
`wfa-high`, `wfa-medium`, `wfa-low`, and `biwfa`. The historical MAF/block positional
form is retained for small diagnostics, but it runs all seven engines in one
process and therefore must not be used for large gaps. This target is not part
of the production command-line interface.

## Interpretation limits

This compact workflow can expose selector regressions and calibrate aggregate
policy using only B73/Mo17. Production `ali` still does not expose a
force-algorithm interface, so its replay output cannot reveal unchosen
counterfactuals. Use the optional non-production exact-candidate benchmark for
that purpose and do not infer engine coefficients solely by grouping replay
tasks by the old source-method label. KSW2 currently reports process-level
timing but not allocator high-water memory; scheduler safety must therefore
continue to use its conservative analytic reservation until finer telemetry is
implemented.
