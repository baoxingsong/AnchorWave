# Fine-grained sequence-alignment parallelism

## Scope

The first implementation keeps anchor alignments and all output reduction in
their original collinear-block order.  Only expensive, two-sided inter-anchor
regions are split into finer tasks.  This preserves the existing alignment
boundaries, scoring, strand handling, MAF records, and fallback algorithms.

Both `genoAli` and `proali` use the implementation.

## Scheduler

The former implementation created one detached operating-system thread per
collinear block.  A chromosome with one large block could therefore use only
one core.  The new implementation uses the existing fixed-size
`AnchorTaskExecutor` for both block and inter-anchor work:

1. Convert consecutive same-strand anchors into lightweight gap descriptors.
2. Select two-sided gaps whose saturated `reference_length * query_length` is
   at least `100,000,000` (approximately two 10-kb sequences), build a compact
   selector plan, and replace geometric cost with predicted risk-adjusted
   runtime for scheduling priority.
3. Submit a worker-sized rolling window per block. Half protects the earliest
   ordered frontier; half prefetches the longest predicted gaps from the whole
   block.
4. Give each result its own task group and future, keyed by the following
   anchor index.
5. Consume each result at the original serial position and replenish both
   channels, so long work starts early while computation overlaps ordered
   anchor output.

Workers never write MAF or BED records.  Existing block code consumes the
aligned strings, score, and method in anchor order and performs the original
validation and output operations. Waiting for one required result does not
wait for unrelated gaps in the same window. Every pending task is drained
before the rolling scheduler is destroyed, including exception paths.

The bounded rolling window is important: submitting every maize gap at once
could keep many completed alignment strings resident even though only a few
WFA calls are active. A slow gap no longer holds all earlier completed results
behind a whole-batch barrier.

The executor nevertheless accounts for every descriptor outside the current
rolling window as planned work. Submission atomically transfers one cost from
the planned backlog to the runnable queue; destruction cancels any unsubmitted
costs. Tail detection and critical-work estimates therefore see the full
remaining block rather than only `-t` visible futures. This matters especially
on 96/128-core nodes, where a worker-sized window previously looked like the
end of the job even when thousands of later gaps remained.

Profiling is performed one gap at a time and stores only the compact selection
plan; extracted sequence strings are released immediately. The actual attempt
reuses that plan, so long-task priority does not duplicate the bounded sequence
profile and does not restore a whole-chromosome sequence cache.

## Memory model

`-w` always defines the maximum algorithm working memory of one alignment
attempt as `w^2`. With `-M`, the process does not assume that every active
thread consumes the maximum; it calculates:

```text
legacy_WFA_budget           = w^2 bytes
safety_reserve              = clamp(2.5% of -M + 64 MiB * -t,
                                    1 GiB, 32 GiB)
task_memory_capacity        = -M - baseline_RSS - safety_reserve
effective_workers           = -t
```

Before each concrete alignment candidate starts, it reserves its predicted
algorithm memory plus transient sequence/output storage. WFA candidates receive
a runtime `max_memory_abort` budget derived from that prediction; an additional
15% guard for WFA2 probe/slab granularity is part of the reservation. The
condition-variable scheduler tracks baseline RSS, active reservations and the
observed retained-RSS residual from output/caches. It admits a task only when
the guarded projection stays within `-M - safety_reserve`, and immediately
reconsiders work when a token is returned. A candidate too large to fit even
by itself is permanently rejected; high current RSS is temporary pressure,
not structural infeasibility.

For the validated command `-t 12 -w 100000 -M 85`:

```text
legacy w^2 WFA budget       = 10,000,000,000 bytes
process limit               = 91,268,055,040 bytes (85 GiB)
safety reserve              = 4,563,402,752 bytes (inside -M)
```

The WFA value is a ceiling, not eagerly allocated memory. Singletrack high and
medium/low WFA are capped at `w^2`; standard high retains its smaller
`min(w^2, 2 GiB)` speculative-trial cap. Full, score-certified and banded KSW2
also remain within `w^2`, while `alignSlidingWindow` keeps its `-w` sizing.
`-M` changes concurrent admission, never the per-engine ceiling. The dynamic
exact-tier ranking still decides which same-quality exact engine runs first.
Failed WFA working sets are
reaped before the next candidate starts. Current RSS is sampled at admission
and release; resident growth outside active reservations therefore reduces
later admission capacity instead of being treated as memory outside `-M`.

High-WFA admission distinguishes three states:

- acquired: run immediately;
- temporarily unavailable: estimate memory wait and compare
  `wait + high-WFA runtime` against the fastest exact fallback runtime;
- structurally infeasible: continue to the remaining exact candidates.

When waiting is predicted to finish first, an eligible gap task is placed in a
timed deferred queue and gives its worker back to the executor. Near the end of
the current executor workload, it instead enters a preferred FIFO reservation
queue and drains enough memory to run high WFA. The same preferred path is
activated after 30 seconds of task aging, preventing large intervals from
starving behind a stream of small reservations. `-bt 0` keeps every
memory-feasible exact engine (exact-first); a positive value applies the same
P50/P90 balanced-mode gate to every exact engine. Temporary pressure alone
never opens the banded or sliding tier.

KSW2 memory is predicted from its traceback representation rather than raw
matrix area. For SIMD traceback, the conservative model is approximately
`16 MiB + (q+r)*(retained_band+160)` bytes. This fixes the earlier `q*r`
estimate, which could under-reserve full KSW2 by about twofold for equal-length
sequences.

WFA2's inner OpenMP team is compiled as an optional execution capability.
Every ordinary genome task still requests one inner thread. The adapter and
benchmark can request 2/4/8 and report the actual team size, but production
must not do so until a global CPU-token scheduler can hold the corresponding
outer slots idle. On a real 14.8/13.7-kb B73/Mo17 gap, 1/2/4/8 threads retained
the same exact score and CIGAR but took 0.549/0.645/0.713/1.310 seconds, so
nested parallelism would be a regression at this scale. It remains a measured
tail-only option for substantially larger critical-path gaps, not a default.

On the 4-Mb B73/Mo17 regression, `-t 12 -w 100000 -M 85` admitted all 12
workers even though 12 fixed `w^2` reservations would not fit. The scheduler
made 244 candidate reservations, reached 12 concurrent reservations, and did
not wait or reject a candidate at admission. After the rolling-window and
retained-RSS fixes, its peak reservation was 30,545,767,316 bytes; peak
projected process memory was 36,440,655,025 bytes and peak sampled RSS was
27,465,449,472 bytes. This is the intended
benefit of task-specific admission: short and similar intervals do not reserve
the worst-case `w^2` amount.

## B73/Mo17 validation

The regression uses the real first 4 Mb of B73 AGPv4 chromosome 1 and Mo17 CAU
1.0 chromosome 1, with the corresponding GFF3, CDS FASTA, and SAM files from
the installation workflow.  It exercises one chromosome-scale block, for
which the previous block-only sequence scheduler reached only one active
thread.

Final command options were `genoAli -t 6 -w 100000 -M 80 -IV` with anchor,
whole MAF, fragment MAF, and method BED outputs.  The threshold selected 68
fine-grained inter-anchor tasks.

| Sequence scheduler | Requested threads | Peak active workers | Wall time | Mean CPU | Peak RSS |
|---|---:|---:|---:|---:|---:|
| previous block-only scheduler | 10 | 1 | 194.33 s | 99% | 10,319,696 KiB |
| fine-grained gap scheduler | 6 | 6 | 100.35 s | 207% | 11,182,032 KiB |

The fine-grained run is 1.94-fold faster on this one-block workload.  Peak RSS
increased by about 8.4% and remained far below the 80-GiB process limit.  The
moderate mean CPU value reflects serial anchor work and unequal gap costs; it
does not indicate that the scheduler created more than six workers.

The following comparisons are byte-identical:

- one-thread versus six-thread whole MAF;
- one-thread versus six-thread fragment MAF;
- one-thread versus six-thread method BED;
- six-thread outputs versus the previous BiWFA implementation;
- normalized anchor files (the command header is excluded because `-t` and
  output paths differ).

Benchmark outputs and timing logs are in:

```text
AnchorWave-wfa2-benchmark-data/installation/sequence-gap-v1-20260808/
```

The same 4-Mb B73/Mo17 inputs were also run through `proali -ns -t 6
-w 100000 -M 80 -R 1 -Q 1`.  That independent command path completed with
exit status zero, reached six active workers, executed 36 fine-grained
inter-anchor tasks, and produced complete whole MAF, fragment MAF, method BED,
and anchor outputs.  Its peak RSS was 14,612,208 KiB.

## Tests

`AlignmentGapPlanner_test.cpp` covers forward and reverse geometry, mixed
strands, overlaps, the high-cost threshold, single-worker bypass, and
overflow-saturating cost estimation.  The complete CTest suite is also run,
and real-data equality checks exercise the task execution and ordered
block-level consumption paths.

`WfaAlignment_test.cpp` additionally exercises in-limit reserve accounting,
dynamic token replenishment, resident-memory drift, process-wide scoped
installation, transient-string accounting, nonblocking temporary/permanent
admission, wait estimation, and preferred FIFO draining.
`AnchorTaskExecutor_test.cpp` verifies that a deferred gap returns its worker,
another queued task runs during the retry delay, and the original task resumes
with its deferral count intact. It also verifies that unsubmitted rolling work
participates in tail/critical-cost state, transfers one task at a time, and is
fully cancelled on success and exception/destruction paths.
