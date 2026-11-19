# IL4RA paired-antibody computational preparation

## Current run: second 500, then combined screening of 1000

The user requested another 500 designs and joint screening of both rounds.
`campaign_local500_round2_20261118/` runs 20 batches of 25, with 125 designs
per original condition and new seeds 2026111900–2026111919. Resolved generation
settings differ from round 1 only in seed and output directory. Current node:
p-gc-3039, existing allocation 57427235. See its `status.json` for live counts.

After generation the controller automatically writes a combined 1000-design
report, up to 20 priority review entries and 20 reserves under
`combined1000_20261118/`. Both rounds receive the same original geometry and
sequence screen; duplicates and failures are retained and labeled. Existing
OpenFold3 evidence is retained without acting as a veto because the positive
control did not calibrate. These are provisional review queues, not final IgGs.

## Current continuation: diagnostics and a conditional 32-variant pilot

This diagnostic continuation ended at 06:33 EST on November 18. Six refolds
completed; the positive control failed calibration. The conditional 32-variant
branch did not run. The later user-requested second 500 campaign above is active.

The 500-design run finished at 05:27 EST on 2026-11-18: 500 unique new
VH/VL pairs, 219 new triage passers, and zero complete computational passers
among 24 new designs plus two imported pilots independently refolded.
The original zero-passer stopping rule fired. All 26 target-aligned Fv RMSDs
were independently recalculated and agree with the original implementation.

`continuation_20261118/` now runs an evaluation-only 6WGL positive control
and a second-seed recheck of four de novo candidates on the same node/allocation.
If the control recovers the known pose and interface in both seeds, a bounded
32-variant partial-diffusion pilot follows, starting only from our own de novo
designs. Candidate thresholds are unchanged. See `plan.json` and `status.json`
in that directory. Known IL4RA binders are never used as generation parents.

## Completed run: 500 designs on p-gc-3039

At the user's request, `campaign_local500_20261117/` generated 500 new designs
(125 per condition) on the current node using allocation `57427235`, step
`57427235.1`. See its `README.md`, `status.json`, and `controller.log` for the
execution plan and live state. The original generation settings and screening
criteria are retained. The earlier 1600-design campaign failed during pilot
oracle startup; its outputs remain preserved in `campaign_36h_v2/`.

## Earlier pilot

This workspace contains the English-language RFD4-Proteina design workflow. The initial paired-antibody job has been submitted and started as Slurm job 57445564 on p-gc-3016. See `reports/design_run.md` for the current method and execution record; the earlier feasibility audit is historical.

Four real paired-antibody inputs are validated. The job runs one representative smoke design, then eight exploratory designs if operational validation passes. The existing RFD4 repository and environments remain unchanged. Independent refolding and final submission screening are not yet complete.

## Files

- `reports/design_run.md`: current design method, resource isolation, and automatic pilot workflow.
- `provenance/submission.json`: submitted job ID, resources, and script hash.
- `provenance/paired_input_validation.json`: real input validation and hashes.
- `reports/feasibility_report.md`: findings, limitations, and execution gates.
- `reports/cpu_checks.json`: reproducible CPU-only check results.
- `reports/validation_summary.json`: current results, including the successful configuration recheck.
- `configs/pilot_generation.yaml`: single-sample generation draft with a mandatory input and explicitly supplied output directory.
- `configs/paired_input.template.yaml`: structural-condition template; requires a reviewed target and paired framework.
- `configs/paired_evaluation.yaml`: explicit target-versus-Fv evaluation draft.
- `configs/resource_plan.json`: proposed budgets and scheduler restrictions.
- `scripts/check_feasibility.py`: bounded CPU-only environment, configuration, and metric checks. It never submits jobs, loads model weights, or enables a GPU.

Run the CPU audit from the RFD4 repository root:

```bash
/storage/group/yjz5549/default/RFD4/RFD4-Proteina/.venv/bin/python \
  /storage/group/yjz5549/default/GenScript_challenge/computational/scripts/check_feasibility.py
```

The script limits itself to one permitted CPU, one numerical-library thread, an 8 GiB virtual-address-space ceiling, a 120-second CPU limit, and a 150-second wall limit. It disables GPU visibility, weight downloads, bytecode writes, and online experiment tracking within its own process. These are preparation limits, not generation budgets.

Use `--config-only` to recheck configuration changes without importing the model dependencies again. This writes `reports/config_checks.json` and preserves the initial full audit. Generation and evaluation require `inputs` and either an explicit `out_dir` override or the project-specific `IL4RA_RUN_DIR` environment variable. Unset output paths must not fall back to a shared default directory.

Configuration composition does not establish that an antibody will be generated correctly. Real structure preparation, annotation round-trip validation, a GPU smoke test, and an oracle smoke test remain required. A successful runtime also does not establish experimental binding.
