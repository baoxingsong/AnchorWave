# Parallel anchor organization in `genoAli`

## Scope

`genoAli -t` now limits both anchor-organization workers and sequence-alignment
workers. The first parallel implementation used one task per chromosome. The
second implementation retains that level and additionally schedules independent
novel-anchor gaps within each chromosome. SAM parsing, blacklist
construction/merge, and anchor-file export keep their existing order and
semantics.

The implementation deliberately does not require `-o`, `-f`, or `-b`. A run
with only `-n` stops after writing the anchor file and never enters MAF sequence
alignment.

## Scheduler and task boundaries

`AnchorTaskExecutor` owns one fixed-size worker pool and a shared priority
queue. Both reference self-chaining and query-genome organization submit one
task per usable chromosome. Idle workers always pull the largest queued task;
the initial cost estimate is the square of candidate-anchor count because the
current `longestPath` implementation is quadratic.

Each query chromosome task performs:

1. anchor sorting and the initial longest path;
2. iterative novel-anchor discovery in inter-anchor regions, with independent
   gaps in one iteration submitted to the global executor;
3. final sorting and longest-path selection.

## Second-level gap tasks

All chromosome tasks and gap tasks share the same `AnchorTaskExecutor`; no
secondary thread pool is created. An `AnchorTaskGroup` records the gap tasks
belonging to one chromosome iteration. A chromosome worker waiting for its
group cooperatively executes work from the global priority queue. Consequently,
ten chromosome parents cannot occupy ten workers and deadlock while their gap
children wait in the queue, and the number of operating-system workers never
exceeds `-t`.

Gap cost is estimated as `reference_gap_length * query_gap_length`. Results are
stored in submission-order slots and consumed only after the task-group barrier,
so minimap2 completion order cannot affect anchor order, blacklist updates, or
the next novel-anchor iteration.

Sequence-bearing gap requests can be large. Each chromosome therefore retains
at most `min(4, worker_count)`, with a lower bound of one, unfinished gap
requests. When the window is full, its parent helps the global executor until a
slot is available. In the B73/Mo17 run an unbounded experimental queue used
5,604,880 KiB peak RSS. The selected four-slot window reduced this to
1,557,120 KiB while also giving a shorter wall time.

Workers write only to preallocated task-local result slots. After all tasks
finish, the calling thread merges results in the original `std::map`
chromosome order. Reference-task results are also merged before constructing
the shared gene blacklist. This barrier preserves the serial dependency
between reference-blacklist generation and query-anchor filtering.

FASTA index maps are copied per query chromosome because the legacy
`getSubsequence2` interface uses non-const `std::map::operator[]`. minimap2's
process-global `mm_verbose` is set once before worker execution. Per-task
minimap2 indexes and thread buffers remain private.

## Removal of the macOS progress stream

The vendored `minimap2/ksw2_extz2_sse.c` contained two macOS-only calls inside
the dynamic-programming row loop:

```c
printf("k");
printf("\b");
```

They did not influence alignment and generated hundreds of MiB of terminal
output in minutes. Concurrent calls also serialized workers on the C stdio
lock. The calls were removed. In the full B73/Mo17 anchor-only test, stdout is
now 319 bytes rather than hundreds of MiB.

## Determinism

Anchor files record the complete command in the first line, so files generated
with `-t 1` and `-t 10` necessarily differ in that line. Equivalence is checked
from line 2 onward:

```bash
diff -q \
  <(tail -n +2 t1.bound4.anchors) \
  <(tail -n +2 t10.bound4.anchors)
```

The canonical SHA-256 for the B73/Mo17 test below is:

```text
24fa6b9fac565f5860a0379ea3b51f065c9f0388724b5835fd92dc6ad6866ade
```

It is identical for the new `-t 1` output, new `-t 10` output, and the
pre-parallelization anchor file.

## B73/Mo17 anchor-only benchmark

The benchmark immediately below predates the v1.3.0 correction of the
32-bit `-fa3^2` overflow. It demonstrates parallel equivalence for that older
threshold behavior, but its 109,413-line hash is not the v1.3.0 release golden
result.

Both runs used the complete B73 AGPv4 reference, Mo17 CAU 1.0 assembly,
31,310 reference CDS features and the installation-workflow SAM files. The
commands used the same biological options (`-w 100000 -IV`) and specified only
`-n` as an output; no MAF option was present.

| Implementation | Run | Wall time | User time | Mean CPU | Peak RSS |
|---|---:|---:|---:|---:|---:|
| chromosome tasks only | `-t 1` | 9 min 27.56 s | 543.24 s | 99% | 1,089,504 KiB |
| chromosome tasks only | `-t 10` | 2 min 49.11 s | 620.46 s | 382% | 1,573,328 KiB |
| chromosome + bounded gap tasks | `-t 1` | 9 min 29.86 s | 544.15 s | 99% | 1,216,384 KiB |
| chromosome + bounded gap tasks | `-t 10` | 2 min 29.31 s | 611.94 s | 426% | 1,557,120 KiB |

The second implementation scheduled 37,107 novel-gap tasks and gives a
3.82-fold `-t 1` to `-t 10` speedup. Its ten-thread run is 11.7% faster than the
first implementation and uses 1.0% less peak RSS. All outputs contain 109,413
lines and have the same canonical hash. The remaining gap to tenfold scaling is
mainly serial SAM/GFF work, the dependent longest-path calculation within each
chromosome, and residual task-size imbalance.

First-version artifacts are under
`AnchorWave-wfa2-benchmark-data/installation/anchor-parallel-20260807/` and
second-version artifacts are under
`AnchorWave-wfa2-benchmark-data/installation/anchor-parallel-v2-20260807/`.

## v1.3.0 post-fix release validation

The complete B73/Mo17 inputs were rerun after the `-fa3` correction. Only
anchor output was requested; `-bt 0`, `-M 102`, and otherwise documented
defaults were used. `proali` additionally used `-R 1 -Q 1`.

| command | threads | wall time | peak RSS | lines | canonical body SHA-256 |
|---|---:|---:|---:|---:|---|
| `proali` | 34 | 1:17.31 | 2,247,984 KiB | 66,744 | `48ed1fc723800644fa81b4667b0f686e98a4a7d275f30ee74a835cfe05c86523` |
| `proali` | 1 | 10:03.39 | 527,236 KiB | 66,744 | same |
| `genoAli` | 34 | 1:11.49 | 2,276,688 KiB | 67,842 | `a12dd94fba69babe3d8da539841c221bb0d6a88fddaf9e11f29f21256a102d14` |
| `genoAli` | 1 | 9:35.10 | 527,812 KiB | 67,842 | same |

Within each command, the one-thread and 34-thread files are byte-identical
after omitting the first provenance line. The changed line counts relative to
the historical table are the expected result of restoring the intended
`reference_length * query_length > fa3^2` condition.

## Tests

`anchorwave_anchor_scheduler_unit` checks:

- enforcement of the global worker limit;
- dynamic replenishment of worker slots;
- largest-estimated-cost-first queue order;
- task exception propagation and executor reuse;
- overflow-safe chromosome and gap cost estimation;
- cooperative nested task groups without parent/child deadlock;
- bounded pending work within a nested group.

Real-data validation additionally requires exact canonical anchor equality
between one and multiple threads and between the first and second
implementations. The complete release test suite and the ASan/UBSan suite each
pass 5/5 tests; the scheduler also passes its ThreadSanitizer test.
