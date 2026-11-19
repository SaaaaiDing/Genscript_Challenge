# IL4RA: 500 designs on the current node

User requested execution on `p-gc-3039`, using existing allocation `57427235`.
This request supersedes the older campaign's node exclusion. No new allocation is submitted.

The existing RFD4-Proteina generation and screening workflow is retained:
canonical IL4RA ECD, paired generic framework, six regenerated CDRs, two patches,
HCDR3 lengths 12 and 16, native BF16, batch size 1, medium/front_loaded
60 steps / 180 trunk passes, pinned checkpoint revision
`a06e2509c7033b010fc43614e9fd840303588502`.

Twenty serial batches of 25 generate 500 new designs: 125 for each of four
conditions. Seeds are 2026111700 through 2026111719. Eight existing pilot
structures are imported for provenance and do not count toward the 500.
Two existing pilots undergo OpenFold3 operational validation before generation.
Four candidates are refolded after 200 and 400 generated designs; after generation,
the original diversity selection and at-most-80 total refold budget apply.
All original scientific thresholds and stopping rules remain unchanged;
the zero-triage check is evaluated at batch boundaries (first possible stop at 175).
Missing oracle evidence stops the controller; it is never counted as passing.

The prior oracle container startup failed because `cache/of3_jax` did not exist.
That directory and the other declared cache directories were created before launch.
The existing environment passed GPU and checkpoint preflight on this node.

Resources: existing standard/normal allocation, account yjz5549_cr_default,
one A100 40 GB, two CPUs, 48 GiB step memory, 35.25-hour total cap,
12-hour generation cap. The existing CPU calculation remains in its own step.
Raw generation extrapolates to 6.1 hours from the pilot, excluding model loading,
triage, and refolding. The 500-design target may be cut short by original stop conditions.

`status.json` records batch-boundary progress; individual outputs appear within
`runs/batch_XX/` while a batch is running. `controller.log` is the launch log.
`oracles/` retains refolds and confidence files. `reports/all_designs.json`
and `exports/computational_review_queue.json` retain screening evidence.

Create `STOP` in this folder to stop at the next controller boundary.
Do not relaunch over existing state without reconciling outputs and the ledger.
This campaign does not submit candidates to the challenge.
