# IL4RA paired-antibody campaign

Submitted Slurm job: **57451616**. This replaces job 57451588, cancelled before GPU oracle startup to correct container GPU mapping. The first attempt is retained in `../campaign_36h/`. This is a computation campaign, not a challenge submission.

## Budget and design rationale

The completed pilot produced one operational smoke sample and eight balanced samples. The eight-sample generation interval was approximately 351 seconds (43.9 seconds per design); the complete Slurm job took 10 minutes 27 seconds. Four of the eight balanced samples passed preliminary geometry and sequence triage. These are not experimentally validated binders.

The production target is at most **1,600 new paired variable-domain designs**: two IL4RA surface patches, two HCDR3 lengths (12 and 16), and 400 samples per condition. All six CDRs are regenerated. Generic paired frameworks are preserved together; the target is canonical human IL4RA residues 26–232. Patch selection is a structural hypothesis, not a verified functional epitope.

Forty serial batches of 40 samples limit recovery cost. Raw generation extrapolates to 19.5 hours, with approximately 21–24 hours allowed for process startup and triage. The 36-hour allocation also includes independent refolding. The controller stops by 35.25 hours and limits the generation stage to the first 24 hours. Two pilot oracle evaluations precede expansion; then four candidates are refolded after every 200 generated samples. Up to 80 total refolds are attempted if time and candidate diversity permit. Actual yield can be lower because of scientific stop conditions, runtime variance, or scheduler preemption.

Resources: account `yjz5549_cr_default`, partition `standard`, QoS `preempt`, one A100 GPU, two CPUs, 48 GiB host memory. Node `p-gc-3039` is excluded. Generation and refolding run serially. No job arrays, automatic requeue, or changes to existing jobs are used. The job can itself be preempted; completed output and provenance remain available.

## Screening protocol v1

1. Validate T/H/L identities, canonical target sequence, all six regenerated CDRs, input-reference integrity, and paired-framework geometry. Check contacts against the intact canonical target aligned to the saved conditioning frame, rather than relying solely on the slightly distorted generated target.
2. Require zero target–antibody heavy-atom pairs below 2 Å, at least two specified hotspots within 5 Å, at least five target contact residues, at least three CDR contact residues, and CDR participation in at least half the contacting antibody residues. Require target-backbone and paired-framework deviations at most 2.5 Å. Require at least 4 Å clearance from annotated glycosylation-site ASN ND2 atoms. This last check is only a glycan-risk proxy; full glycans are not modeled.
3. Exclude additional CDR cysteines. Record CDR N-linked glycosylation motifs; zero motifs are required for the final computational review queue. Deduplicate paired VH/VL sequences and select diverse candidates using all six CDRs, with greedy normalized edit similarity below 0.8 between selected candidates.
4. Compare IMGT CDR sequences with a frozen TheraSAbDab IL4R reference set containing eight numbered therapeutic antibodies. Flag candidates with both six-CDR similarity and HCDR3 similarity at least 0.8. This is a conservative internal similarity rule, not proof of legal novelty or an exhaustive literature/patent search. References are used only for exclusion, not generation.
5. Refold the T/H/L complex with OpenFold3, templating the target only. Retain original and predicted structures and confidence JSON. Recalculate target-aligned paired-Fv RMSD and separate VH/VL confidence. Require paired-Fv RMSD at most 2.5 Å, VH and VL mean CA pLDDT at least 80, and refold geometry passing the same checks. At least one target-contacting chain must have its own target-pair minimum PAE at most 2 Å and pair ipTM at least 60 on the repository's 0–100 scale. Require H–L ipTM at least 60 as a separate pairing check; it cannot substitute for target-interface evidence.
6. Export at most 20 diverse computational review candidates. Before challenge submission, these still require broader novelty review, detailed developability and nonspecificity assessment, glycan and full-IgG context review, and complete HC/LC assembly. No experimental affinity or specificity is inferred from these gates.

Missing or malformed oracle evidence stops expansion rather than being counted as passing. Zero triage passers among 160 new designs, or zero computational passers after at least 24 complete refolds, stops expansion for scientific review. Thresholds are fixed in `criteria_v1.json`; automatic threshold relaxation is prohibited.

## Evidence and operation

- `plan.json`: frozen budget and scheduler settings.
- `submission.json`: Slurm receipt and submission-time script hashes.
- `campaign.sqlite`: input/config/run/design/oracle/metric ledger.
- `events.jsonl`: append-only worker events, ingested idempotently.
- `configs/`: resolved generation configuration for each planned batch.
- `status.json`: live counters, current stage, errors, and stopping reason; created once the job starts.
- `reports/all_designs.json`: individual triage and refold evidence, including failures.
- `exports/computational_review_queue.json`: eligible diverse candidates, possibly empty.
- `runs/` and `oracles/`: command logs, generated structures, refolds, confidence files, and runtime configuration.

Create `computational/STOP` or `campaign_36h_v2/STOP` to stop at the next batch boundary. Existing artifacts must be reconciled before any restart; the controller deliberately refuses silent duplicate initialization or submission.

The campaign uses the repository's working `.venv`, the pinned local RFD4 checkpoint, and existing local OpenFold3 assets. The alternate conda environment was not modified.
