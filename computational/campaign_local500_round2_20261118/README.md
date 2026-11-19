# Second 500 designs and combined selection

User requested another 500 designs on the current node, followed by unified
screening of the 1000 designs. Run on p-gc-3039 using existing allocation 57427235.
One A100, two CPUs, 48 GiB step memory, 16-hour step limit; generation is capped
at 12 hours. All model work is serial. The earlier CPU calculation is preserved.

Twenty batches of 25 use seeds 2026111900–2026111919, disjoint from round 1.
Each of the original four conditions receives another 125 samples. Model weights,
canonical target, paired framework, all six CDR conditions, BF16, batch size 1,
60 steps / 180 trunk passes and sampling parameters are unchanged. This is full
de novo generation, not the conditional partial-diffusion pilot.

The previous continuation completed its six diagnostic refolds; the positive
control failed calibration, so its 32-variant branch was never started. This new
user-requested generation campaign does not depend on that failed calibration.
No additional OpenFold3 run is scheduled as a mandatory gate. Existing oracle
results remain attached to the first-round candidates as uncalibrated evidence.

`status.json` refreshes during each batch, including completed batch counts and
the number of newly written structure files. A written file is not a validated
candidate until its batch is collected. `runs/batch_XX/execution/process.log`
contains model progress. `campaign.sqlite`, `events.jsonl`, resolved configs,
input hashes and `launch_sources/` retain provenance. Create `STOP` here to stop
before the next batch. Existing state must be reconciled before restarting.

After exactly 500 new designs are accounted for, the controller automatically
runs `scripts/screen_combined1000.py`. It excludes imported pilots and controls,
preserves all failures, deduplicates VH/VL pairs across rounds and applies the
same original geometry/sequence triage to both rounds. Zero CDR glycosylation
motifs and no frozen-reference similarity flag are required for the review pool.
The uncalibrated oracle is not a veto or reward, and untested is not labeled failed.

Within each condition, candidates are prioritized by CDR contact fraction, CDR
contact count, paired-framework RMSD and glycan-site clearance. These are review
heuristics, not a calibrated prediction of affinity. Selection cycles through
the four conditions with six-CDR normalized edit similarity <0.8, yielding up to
20 priority review candidates and 20 reserves, possibly fewer. Duplicate and
excluded rows remain visible with reasons. No official criterion is invented.

Outputs: `../combined1000_20261118/all1000_ranked.csv`, `summary.json`,
`priority_review_20.csv`, `reserves_20.csv`, and complete JSON evidence.
Any selected VH/VL still needs broader novelty, developability, glycan/Fc context,
and full-IgG assembly review before submission. No challenge order is placed.
