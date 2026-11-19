# Completed IL4RA paired-antibody campaign

Both 500-design rounds are complete. The combined 1000 designs contain 998
unique VH/VL pairs. Follow-up diagnostics and the final sequence/assembly review
are complete; the final queue contains 20 priorities and 20 reserves.

Start with the [repository overview](../README.md) and the
[final IgG review](igg_review_20261118/README.md).

- [Round 1](campaign_local500_20261117/): 500 generated, 219 original triage passers.
- [Round 2](campaign_local500_round2_20261118/): 500 generated, 186 original triage passers.
- [Diagnostic continuation](continuation_20261118/): six refolds; the known positive
  control failed calibration. The conditional 32-variant generation branch did not run.
- [Combined first-pass screening](combined1000_20261118/): 237 unique review candidates.
- [Additional review and IgG assembly](igg_review_20261118/): 211 candidates after
  conservative additional filters; final 20 priorities and 20 reserves.

Original strict OF3 complete-pass count remains zero. Failed positive-control
calibration means OF3 evidence is retained without a veto or reward. Untested
candidates remain explicitly untested. The final queue is a computational
experimental-testing priority, not proof of binding or full IgG 3D assembly.

Scripts, configurations, process logs, metrics and failure records are preserved.
CIF structures, SQLite ledgers and the large numbering cache are losslessly
stored in ../archives/; restore them with ../repository_tools/restore_artifacts.py.
Credentials, external model weights and environment caches are excluded.

The [previous overview](../docs/history/computational_readme_before_archive.md)
is retained as historical context; its live-status language is superseded here.
