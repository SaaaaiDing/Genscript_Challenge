# Combined review of 1000 IL4RA paired-Fv designs

Planned automatic output after the second 500-design run completes. Until
`summary.json` reports `status: completed`, a combined 1000-design screen has
not finished. Live generation state is in
`../campaign_local500_round2_20261118/status.json`.

The combined pool contains exactly 500 designs from each generation round,
excluding historical pilots and known controls. Unique sequence-pair counts
may be below 1000; cross-round duplicates remain visible and cannot occupy
multiple shortlist slots.

The screen preserves original geometry/sequence triage, requires zero CDR
glycosylation motifs and no frozen-reference similarity flag, and balances
four design conditions with six-CDR diversity. Existing OpenFold3 results are
retained as uncalibrated evidence; neither a missing result nor a failed result
is automatically treated as a biological negative. Original full-gate labels
are preserved separately from the new review-eligibility label.

Artifacts produced automatically:

- `all1000_ranked.csv` and `all_designs.json`: every design, sequence, reason,
  provenance, duplicate label and available oracle evidence.
- `priority_review_20.csv/.json`: at most 20 diverse review priorities.
- `reserves_20.csv/.json`: at most 20 additional diverse reserves.
- `summary.json` and `policy.json`: exact counts and documented decision rules.

The ranking is a geometry/sequence heuristic, not predicted KD or a probability
of experimental binding. Selected pairs still need broader novelty, detailed
developability and full IgG context/assembly review before challenge submission.
