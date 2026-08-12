# WFA2/BiWFA upgrade and validation

The current resource-aware preselector, including its algorithms, telemetry,
and complete B73/Mo17 4 Mb evaluation, is documented in
[`ALIGNMENT_ALGORITHM_SELECTOR.md`](ALIGNMENT_ALGORITHM_SELECTOR.md).

> Production status in v1.3.0: BiWFA has been removed from the `genoAli`,
> `proali`, and `ali` selector. It remains only in the low-level adapter and
> developer benchmark for counterfactual measurements. Historical sections
> below document the upgrade path, not the current candidate chain.

## Scope and provenance

AnchorWave vendors WFA2-lib directly. This upgrade synchronizes the complete
`WFA2-lib/` tree with the official `smarco/WFA2-lib` repository at:

```
main@2fc4a1afac0f624e7020ee5aadb7692b38157eaa
VERSION.txt: v2.3
```

The version string is still `v2.3`, but this upstream `main` is substantially
newer than the previously vendored `v2.3` snapshot. It adds the current public
`wavefront/wfa.h` interface, sequence and termination modules, extension
kernels, CIGAR utilities, and newer BiWFA implementation/fixes. AnchorWave does
not change WFA2's scoring objective. The adapter is implemented in
`src/myImportandFunction/WfaAlignment.cpp`; the production selection and
resource policy live in `AlignmentAlgorithmSelector.cpp`,
`AlignmentResourceScheduler.cpp`, and `alignSlidingWindow.cpp`.

Small, audited fixes are carried on top of that pinned revision. Four
portability fixes were exposed by Clang
AddressSanitizer/UndefinedBehaviorSanitizer on arm64:

- copy potentially unaligned 64-bit sequence words with `memcpy` in the match
  extension, packed-CIGAR, and backtrace paths;
- encode negative piggyback values through their unsigned representation before
  shifting, instead of left-shifting a negative signed integer;
- avoid constructing a pointer one element before the allocator request array
  when its last request is removed.

These changes remove C-language undefined behavior without changing the WFA or
BiWFA algorithm. They are deliberately narrow so that a future upstream sync is
easy to audit.

One additional integration fix restores the advertised memory-limit behavior
for BiWFA. Upstream's classic WFA score loop periodically calls its step/memory
limit probe, but the separate bidirectional score loop did not. AnchorWave calls
that same upstream probe after each forward/reverse score advance (the probe
itself remains interval-controlled). Without this fix, BiWFA inherited
`max_memory_abort` but never evaluated it during the main bidirectional loop.

The initial integration changed seven WFA2 files for the portability and BiWFA
probe fixes above. The current v1.3.0 tree also contains audited integration
changes for Singletrack/backtrace memory modes, cooperative runtime and memory
probes, and optional OpenMP telemetry. Therefore the present vendored-tree diff
is intentionally larger than that initial seven-file patch set.

## Current alignment policy

For every non-filling anchor or inter-anchor interval, AnchorWave now follows
this order:

1. Convert `-w` to the historical per-worker byte budget `w^2`. This is the
   only information derived from `-w` that is passed into the WFA adapter.
2. Profile the oriented sequences and retain the memory-feasible exact
   candidates: Singletrack high WFA, standard high WFA, medium/low WFA,
   score-certified KSW2, and full KSW2. Runtime ranking, rather than this
   deterministic tie-break order, chooses among them. BiWFA is not admitted.
3. `-bt 0` selects exact-first mode and keeps every memory-feasible exact
   candidate. A positive `-bt` selects balanced mode (`-bt` itself defaults to
   zero) and
   applies the same calibrated time gate to every exact candidate: every
   admission quantile, currently P50 and P90, must be within `-bt`. The same
   value is also a cooperative cumulative runtime deadline for exact attempts
   on one interval.
4. Cap every WFA/KSW2 attempt at `w^2`, including both high-WFA modes. With
   `-M`, compare predicted wait plus runtime among same-quality exact
   candidates and park/requeue the gap task when waiting wins.
5. In the workload tail, or after task aging, use a preferred FIFO reservation
   to drain memory for high WFA and prevent starvation. Reap a failed WFA
   working set before trying the next candidate.
6. Score-certified KSW2 accepts a banded traceback only when its score equals
   an unbanded score-only global optimum; otherwise it doubles the band within
   its guarded resource cap and then continues to another exact candidate.
7. Only after the admitted exact tier is exhausted or its explicit runtime
   deadline is reached, enter the shared fallback tier. Predict whether banded
   KSW2 or sliding-window alignment will return the larger score and normally
   execute only that method; runtime breaks close predictions.
8. Record the successful method in the optional BED output. All production WFA
   modes, including Singletrack and medium/low, are exported as `WAVEFRONT`;
   internal traces retain their detailed names. Other labels remain
   `MINIMAP2`, `BANDED_MINIMAP2`, and `SLIDING_WINDOW`.

The former rule that compared sequence length/difference/score directly with
`-w` remains removed: `-w` is not a WFA work threshold. The separate `-bt`
parameter is the only predicted-time cutoff and, in balanced mode, applies to
all exact engines before execution.
The standard WFA modes and benchmark-only BiWFA adapter use the same pinned
WFA2 library, so no duplicate legacy WFA2 source tree is required. Failed WFA
working sets are reaped before another candidate is attempted; successful
large retained sets are reaped after alignment. All-`N` filling regions retain
their existing behavior.

## Meaning of `-w`

The public meaning and default of `-w` are preserved for `genoAli` and
`proali`:

```
approximate per-thread alignment budget in bytes = w * w
```

Thus, the default `-w 100000` represents 10,000,000,000 bytes, or about
9.31 GiB, for one alignment engine. The benchmark-only BiWFA adapter divides
that aggregate target among its subsidiary WFA objects and invokes WFA2's
memory probe in its bidirectional score loop.

`-w` remains a per-attempt algorithm-memory ceiling and the matrix/window bound
for KSW2 and the sliding-window fallback. Inside WFA it is used only as memory:
Singletrack high, standard high, and succinct modes receive guarded,
prediction-sized budgets up to `w^2`. `-M` controls process-wide concurrency
but cannot expand an individual
attempt beyond `w^2`. Sequence length, length difference and
alignment score are not directly compared with `-w`. Values must be positive
and the `w*w` multiplication is overflow-safe.

BiWFA's lower-memory behavior remains useful for developer comparisons, but
its long-tail runtime is the reason it is excluded from the v1.3.0 production
selector. The same `w^2` ceiling applies when the benchmark API invokes it.

## Process-wide memory and dynamic thread scheduling

`genoAli` and `proali` accept two complementary resource limits:

```text
-t INT     maximum number of alignment workers
-M FLOAT   maximum process-wide memory reservation in GiB
```

When `-M` is present, it must represent at least `w^2` bytes. AnchorWave reads
its resident set after the reference and query genomes have been loaded and
deducts both that baseline and an in-limit safety reserve:

```text
safety reserve       = max(1 GiB, 5% of -M)
task capacity        = -M - baseline RSS - safety reserve
maximum workers      = -t
```

Each concrete alignment attempt reserves its predicted algorithm peak,
transient sequence/output storage, and a WFA allocator/probe guard. A
condition-variable scheduler compares the greater of sampled RSS and the
baseline-plus-reservations projection with the usable limit. As soon as a task
returns its token, waiting work is reconsidered. Fine-grained gaps can be
parked in a deferred queue so temporary memory pressure does not occupy a
worker. Tail/aged high-WFA requests enter a preferred FIFO queue that stops
younger reservations from consuming the memory they need. This permits many
small alignments to use all requested workers while preventing several large
tasks from independently assuming the same process memory is available.

The cap is an application-level admission limit rather than a portable OS
`RLIMIT_RSS`. Prediction margins, current-RSS sampling, and the safety reserve
are all accounted within `-M`; allocator growth therefore reduces subsequent
admission capacity rather than being treated as memory outside the setting.
Omitting `-M` preserves legacy `-t`-only scheduling.

## Result-equivalence rules

Standard WFA and BiWFA are exact and should return the same optimal score when
no heuristic is enabled. WFA2 documents one important distinction: BiWFA is
optimal but does not guarantee the same tie-breaking path as standard WFA.
Consequently, a byte-for-byte MAF difference can be valid when repetitive
sequence permits multiple equal-scoring gap placements.

Validation uses three levels:

- `EXACT`: block coordinates, score, and aligned text are identical.
- `COORDINATE_EQUIVALENT`: block coordinates, ungapped sequences, and score
  are identical, but an equal-scoring gapped path differs.
- `DIFFERENT`: coordinates, ungapped sequence, score, or block count differs.

`workflow/compare_maf.py` implements this comparison in streaming mode and
validates every MAF row, so multi-gigabyte outputs do not have to be loaded into
memory. Its strict default exits nonzero for `DIFFERENT`; the benchmark uses
`--report-only` so a valid difference is fully reported and does not prevent the
fragment MAF comparison. Parse or I/O errors remain fatal.

For the default AnchorWave penalties, the benchmark also passes
`--score-penalties 6,8,2,75,1`. This independently recomputes every pairwise
fragment score from the exported alignment text, reports inconsistencies
between a MAF `a score=` value and its rows, and counts which implementation has
the better recomputed score. This is important for distinguishing a BiWFA
regression from a pre-existing fallback score/reporting discrepancy.

## Tests

Configure and run the code/unit tests:

```bash
cmake -S . -B /private/tmp/anchorwave-wfa-new-build \
  -DCMAKE_BUILD_TYPE=Release
cmake --build /private/tmp/anchorwave-wfa-new-build -j 8
ctest --test-dir /private/tmp/anchorwave-wfa-new-build --output-on-failure
```

The WFA-adapter-specific tests cover:

- the exact `-w -> w^2` mapping, standard-WFA trial ceiling, BiWFA component
  budget, invalid values, and overflow;
- identical sequence alignment;
- two-piece affine score reconstruction from the returned alignment;
- optimal-score equality with standard exact WFA on 100 seeded sequence pairs
  containing substitutions, insertions, and deletions;
- CIGAR consumption and reconstruction of both ungapped inputs;
- a strongly length-imbalanced interval handled without a `-w` admission rule;
- bounded standard-WFA fast-path selection;
- independent standard-WFA and BiWFA thread-local cache reuse;
- simultaneous alignment in six worker threads, each with isolated thread-local
  state;
- a deliberately tiny memory budget that must return `MemoryLimit` with no
  partial alignment.

The MAF comparator tests cover exact, tie-equivalent, biologically different,
invalid, and nonfatal report-only cases. Existing TrioAnchorGraph C++, CLI, and
workflow tests remain in the same CTest run. The focused WFA adapter suite now
contains 15 tests. All six CTest targets pass in the release build on arm64
macOS, and all 15 focused tests pass under combined AddressSanitizer and
UndefinedBehaviorSanitizer with leak detection disabled on this platform.

An integration check runs `anchorwave ali -w 1` on the 20 kbp stress pair. It
forces BiWFA's memory-abort path and then exercises the existing fallback; the
two ungapped output sequences reconstruct both inputs exactly.

### Same-library WFA versus BiWFA microbenchmark

`benchmarks/WfaMemoryModeBenchmark.cpp` compares the two exact modes from the
same pinned WFA2 library. It generates deterministic sequence pairs, reuses one
aligner per mode, uses the same two-piece-affine penalties, checks the completed
count and score checksum, and reports WFA2's working-memory counter. Build and
run it with:

```bash
cmake -S . -B /private/tmp/anchorwave-wfa-bench-build \
  -DANCHORWAVE_BUILD_TESTS=OFF \
  -DANCHORWAVE_BUILD_BENCHMARKS=ON \
  -DCMAKE_BUILD_TYPE=Release
cmake --build /private/tmp/anchorwave-wfa-bench-build \
  --target anchorwave_wfa_mode_benchmark -j 8
/private/tmp/anchorwave-wfa-bench-build/anchorwave_wfa_mode_benchmark \
  10000 0.10 5 1
```

Representative single-thread arm64 results are:

| length | simulated error | alignments | standard WFA | BiWFA | WFA speedup | WFA / BiWFA peak working memory |
|---:|---:|---:|---:|---:|---:|---:|
| 1,000 | 5% | 1,000 | 0.271 s | 4.009 s | 14.8x | 10.1 MB / 6.9 MB |
| 5,000 | 20% | 10 | 1.057 s | 6.072 s | 5.7x | 978.5 MB / 27.0 MB |
| 10,000 | 10% | 5 | 0.562 s | 3.456 s | 6.1x | 1.28 GB / 39.0 MB |
| 10,000 | 20% | 1 | 0.553 s | 1.624 s | 2.9x | 3.49 GB / 56.0 MB |

Scores and completed-alignment counts were identical in every row. These
historical data confirmed the standard-WFA speed advantage but do not define
the v1.3.0 production candidate order or current memory thresholds.

### Focused memory stress check

A deterministic 20,000 bp pair with 1,000 substitutions (5%) was aligned with
`ali -w 100000` on the same arm64 macOS host. This is a mechanism check, not a
replacement for the whole-genome benchmark. GNU time reported:

| implementation | peak RSS (KiB) | wall time | output |
|---|---:|---:|---|
| previous standard WFA (`memory_high`) | 722,080 | 0.11 s | SHA-256 `73b172fb...53af` |
| upgraded BiWFA (`memory_ultralow`) | 52,672 | 0.41 s | SHA-256 `73b172fb...53af` |

In this single run, BiWFA used about 13.7-fold less peak RSS and was about
3.7-fold slower. The aligned FASTA outputs were byte-identical. The speed/memory
tradeoff is expected to change with interval length, divergence, gaps, and CPU;
the B73/Mo17 benchmark below is the decision-grade comparison.

## Historical BiWFA-first B73 v4 versus Mo17 CAU benchmark

The following measurements describe the superseded BiWFA-first policy and are
retained as the evidence that motivated the adaptive standard-WFA/BiWFA policy.
They are not performance measurements of the current coexistence selector.

`workflow/benchmark_wfa2_upgrade.sh` reproduces the data and commands in
`installation.md`, while running the old and new executables serially for a
fair resource comparison. Preparation uses all logical CPUs for minimap2;
`genoAli` defaults to one thread. Terminal progress stdout is discarded because
the authoritative products are the anchor and MAF files.

```bash
OLD_BIN=/private/tmp/anchorwave-wfa-old-build/anchorwave \
NEW_BIN=/private/tmp/anchorwave-wfa-new-build/anchorwave \
ANCHORWAVE_THREADS=1 \
WINDOW_WIDTH=100000 \
workflow/benchmark_wfa2_upgrade.sh all
```

The workflow is resumable by stage:

```bash
workflow/benchmark_wfa2_upgrade.sh download
workflow/benchmark_wfa2_upgrade.sh prepare
workflow/benchmark_wfa2_upgrade.sh provenance
workflow/benchmark_wfa2_upgrade.sh run-old
workflow/benchmark_wfa2_upgrade.sh run-new
workflow/benchmark_wfa2_upgrade.sh compare
```

On macOS the workflow uses `caffeinate` when available, scoped to the script
PID, so multi-hour stages are not interrupted by system sleep. This helper is
outside the GNU-time measurement.

Results include GNU-time peak RSS and wall time, SHA-256 checksums, exact anchor
file comparison, old/new alignment-method counts from `-b`, a streaming
comparison of final MAFs, and a separate comparison of fragment MAFs. A
`.complete` sentinel is created only after each `genoAli` run exits
successfully and all temporary products are moved into place.
`benchmark-provenance.txt` records both executable hashes, the repository and
WFA2 revisions, thread counts, `-w` budget, minimap2 path, and host information.

The full B73/Mo17 test is intentionally not run concurrently with the paused
26-genome workflow. It is executed old first and new second, each with one
AnchorWave thread.

The installation data were downloaded and preparation completed: 31,310 CDS
features were mapped to B73 and Mo17. A full old-version `genoAli -t 1` probe
spent 10 minutes in global anchor organization without producing alignment
outputs; extrapolation was approximately 18 hours. It was deliberately stopped
under the resource-conservation requirement, so no full-genome equivalence or
performance result is claimed. The inputs, logs, and resumable workflow remain
available for a later long-running benchmark.

For an executable real-data check with bounded cost,
`workflow/benchmark_wfa2_bounded.sh` applies the same installation workflow to
the first 4 Mb of chromosome 1 by default. This interval includes the initial
39,609 bp versus 167,149 bp pathological length-imbalance case as well as
hundreds of ordinary intervals. It reuses the official input data, runs old and
new versions serially, and produces the same time, checksum, method-count, and
streaming MAF-comparison reports. Set `WINDOW_BASES` to request another prefix.
This bounded comparison is a regression test, not a substitute for the
deferred whole-genome measurement.

An exploratory 20 Mb new-version run was also stopped after it reached 4.5 Mb.
It showed the intended tradeoff on a 88,758 bp versus 67,651 bp interval with
optimal score -111,742: the old full-CIGAR WFA exhausted its budget and used
minimap2, while BiWFA completed the exact interval with much lower working
memory but required several minutes of CPU time. The incomplete outputs and GNU
time report are retained with the suffix `20mb-exploratory`; they are not used
as comparison results.

### Completed 4 Mb bounded result

The bounded old/new runs completed successfully with `genoAli -t 1 -w 100000`
on the same arm64 macOS host:

| implementation | wall time | peak RSS (KiB) | peak RSS (GiB) |
|---|---:|---:|---:|
| previous standard WFA | 1 min 16.90 s | 14,746,736 | 14.06 |
| upgraded BiWFA policy | 3 min 11.96 s | 10,316,016 | 9.84 |

The new run reduced peak process RSS by 30.0% (1.43-fold lower) and took
2.50-fold longer. The peak is dominated by the bounded fallback rather than a
successful BiWFA interval, illustrating why `-w` is an approximate per-worker
target rather than a process-wide hard cap.

The semantic anchor rows are identical; only the first provenance/command line
differs because it contains the `old.*` or `new.*` output paths. Both MAFs are
structurally valid, contain 253 fragments and 506 sequence rows, preserve all
coordinates, and reconstruct exactly 8,000,000 ungapped bases. Of the 253
fragments, 198 are byte-exact and 55 use a different gapped path.

The method distribution changed from 198 `WAVEFRONT`, 46 `MINIMAP2`, and 9
`FILLING` intervals to 241 `BIWFA`, 3 `MINIMAP2`, and 9 `FILLING` intervals.
The three new-version fallbacks are precisely the intervals whose unavoidable
length-difference score exceeds 100,000; thus 95.3% of all intervals use exact
BiWFA while the three pathological intervals remain bounded.

Four fragments have different declared scores. Independent recomputation from
the exported alignments finds 249 equal scores, four cases where the BiWFA
alignment is better, and zero cases where the old alignment is better. The old
fragment MAF has five declared-score inconsistencies: one shared minimap2
fallback inconsistency plus the four cases that old WFA could not complete.
The new fragment MAF has only the shared fallback inconsistency; every exported
BiWFA block's declared score matches its alignment text. Consequently the MAF
comparison correctly reports `DIFFERENT`, but the difference consists of
equal-coordinate alternative paths and four verified score improvements, not
lost sequence or changed anchors.

The original fixed-reservation validation used `-t 10 -M 80 -w 100000` and is
retained as a historical baseline. The newer task-memory scheduler was tested
on the same 4-Mb B73/Mo17 input with `-t 12 -M 85 -w 100000`: it admitted 12
concurrent predicted-memory tasks, made 244 reservations without waiting, and
reported a 24,704,696,320-byte peak sampled RSS. Its safety reserve was
4,563,402,752 bytes and was deducted inside the 85-GiB limit.
