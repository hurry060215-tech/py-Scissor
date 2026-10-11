# Full-LUAD Cox recovery after executor-key change

The original tightened Cox execution channel failed on 2026-10-10 with an
executor-key-change transport error. Its lifecycle was not inferred from `ps`:
each command has an isolated process namespace. No terminal exit code was
received, and the old model directory was left unchanged.

At interruption, alpha 0.05 had its complete full path, five folds and refit.
Alpha 0.2 had complete full, fold 1 and fold 2 paths. Its fold 3 was unfinished;
alpha 0.6 had not started. Thus 10 complete nodes could be copied and 11 needed
computation. Every node still uses all applicable samples and all 4,102 cells.

`prepare_isolated_resume.py DESTINATION` snapshots this exact recovery boundary
to a previously nonexistent output root. It checks terminal flags, shapes,
finite values and binary sizes, then makes independent regular-file copies.
Every source/copy hash is checked before and after copying. No hard links or
partial paths are accepted. `isolation-snapshot.json` preserves the old writer's
unknown lifecycle, the source log hash/mtime, all copied hashes, input-manifest
hash and outstanding work. The script is specific to this recorded experiment,
not a general-purpose checkpoint discovery utility.
An optional `--original-session-id` records the local execution channel in the
private snapshot; the public script has no hard-coded channel identifier. The
actual script version used for this experiment is separately preserved with its
hash in the private output's `executed-sources` directory.

`resume_isolated_corrected.R` is a separately named variant of the frozen
ordinary recovery driver. Run `validate` first to independently recompute every
copied coefficient's exact-objective and absolute-KKT certificates against the
same full design, graph, outcomes, folds, alpha/lambda paths and 1e-10 threshold.
Only after that succeeds, run `resume` with the final argument
`confirmed-isolated-new-output` and the exact derived-manifest SHA-256. This
argument does not assert the original writer is dead. Writes remain confined to
the new output root. Existing complete alpha results are preserved bytewise.

Before each missing full/fold/refit node is scheduled, the driver checks the old
source log's hash and mtime. If either advances, it stops new scheduling without
killing the unknown old process. A currently running new node may finish and
save its atomic checkpoint before that check. The unfinished fold is recomputed
from the beginning of the full explicit lambda path; no partial optimizer state
is reused. Results record isolated-resume provenance and reused/new components.

`selftest_isolation_guard.R` extracts the actual scheduling guard from that
driver and verifies unchanged-source acceptance plus refusal of both content
changes and mtime-only changes. The earlier `selftest_resume.R` tests the shared
node reuse/recomputation and numerical equivalence mechanisms.

The isolated alpha 0.2 run subsequently completed every missing node and passed
the unchanged Python/R full-path, five-fold and cold-refit comparison: maximum
coefficient difference 3.0013e-7, no stable sign differences, and identical
lambda.min/lambda.1se. Alpha 0.6 also completed all seven nodes and passed the same comparison, with
maximum coefficient difference 1.4131e-6, no stable sign differences and exact
CV choices. The new recovery driver exited 0 after 5,116.55 seconds (peak RSS
462,114,816 bytes). These are successful new-root results; the original channel
remains of unknown lifecycle.

All new sources are GPL-3.0-or-later; see `../../LICENSE`.
