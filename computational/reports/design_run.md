# Initial IL4RA paired-antibody design run

This report supersedes the earlier scheduler hold in `feasibility_report.md`.

## Method

Use the existing RFD4-Proteina runtime and the operational settings from the user's TAS1R2/TAS1R3 campaign: native BF16 inference, batch size 1, no data-loader workers, eager execution, and medium/front-loaded recycling (60 steps, 180 trunk passes). Generation and expensive refolding are separated.

The starting antibody structure is the human germline IGHV3-23/IGKV1-39 Fab, PDB 5I19. Only the paired variable-region framework is retained. All six IMGT CDR intervals are removed and replaced with generated segments. No known IL4RA-binding antibody is used as a starting binder. The source Fab contains a previously used H3; that H3 is explicitly removed along with the other five CDRs, not retained or mutated.

The target is the canonical P24394-1 sequence, residues 26-232, using the complete AlphaFold v6 ECD backbone. Its experimentally observed core aligns to 3BPN at 1.1185 Angstrom C-alpha RMSD. The experimental 3BPN construct has engineered substitutions and incomplete termini, so it is used as a structural cross-check rather than silently treating it as the canonical antigen.

Four input problems combine two exposed surface patches with H3 lengths 12 and 16. Hotspots are native/UniProt residues 107/108/109 and 67/70. The machine-readable condition specification uses the original source numbering; exported contigs renumber residues. Do not interpret exported residue IDs directly as UniProt positions without applying the mapping.

Surface selection uses residue SASA, target pLDDT, distance from known glycosylation sites, and a terminal-position filter. This is an initial accessibility heuristic, not explicit glycan/Fc modeling or an established optimal epitope. Framework/CDR novelty and independent complex refolding remain required.

Sources: [5I19](https://www.rcsb.org/structure/5I19), [human germline framework library](https://pmc.ncbi.nlm.nih.gov/articles/PMC4968113/), [3BPN](https://www.rcsb.org/structure/3BPN), [UniProt P24394](https://www.uniprot.org/uniprotkb/P24394/entry), [AlphaFold model](https://alphafold.ebi.ac.uk/entry/P24394).

## Input validation

Four annotated CIFs passed the actual RFD4 pipeline validator, index-group/molecule consistency checks, residue-gap checks, and six-CDR condition checks after reload. All six CDR placeholders are sequence-unconditioned and internally coordinate-unconditioned. The paired framework uses one shared internal-coordinate group, separate from the target group. Files and hashes are recorded in `provenance/paired_input_validation.json`.

The authoring pipeline reports copied auxiliary annotations on newly inserted sentinel atoms. These are placeholders, not native residue mappings; scientific chain/CDR provenance is retained separately. The model validation path passed. This does not establish exhaustive chemistry correctness or successful antibody generation.

## Execution and resource protection

Slurm job **57445564** was submitted on 2026-11-17 and started on **p-gc-3016** at 17:12:24 EST. It explicitly excludes the protected node **p-gc-3039** and requests one A100, two CPUs, 24 GiB host memory, and one hour. Numerical-library threads are restricted to one. The existing MMGBSA processes and their job step were not changed.

The job uses QoS `preempt`, which has no configured preemption targets. It may itself be preempted. No elevated priority, overlapping GPU sharing, automatic requeue, or job array is used. Although account association output listed only `normal`, the cluster does not enforce QoS associations through `AccountingStorageEnforce`; a non-submitting scheduler test accepted `preempt`, and the actual job was accepted and started with that QoS.

The upstream launcher rendered the batch script through a non-submitting shim. The project-owned script was then checked, given an explicit QoS directive and runtime guard, changed to one numerical-library thread, and submitted through the real Slurm command. This avoids modifying the shared RFD4 repository. The final script hash and job ID are retained in `provenance/submission.json`.

## Automatic sequence of work

1. Run the repository GPU/checkpoint preflight inside the new allocation.
2. Generate one representative H3-16 paired-antibody smoke candidate.
3. Require finite coordinates, correct T/H/L chains, plausible VH/VL lengths, generated antibody positions, and canonical target identity.
4. If those operational checks pass, generate eight exploratory designs: four inputs, two replicates each, serially at batch size 1.

Any failed preflight or smoke check stops the job before expansion. A `computational/STOP` file prevents starting the next major phase; it does not interrupt an already running sample. It is not a live per-sample stop mechanism.

Outputs are under `runs/representative_smoke/` and `runs/initial_designs/`. Logs are under `runs/initial_designs/slurm_run_outputs/`. Raw generated targets and saved conditioned reference coordinates must be distinguished, as established in the source sweetness campaign. Predicted target collapse is a scientific quality failure, not evidence that the reference input was corrupted.

Generated VH/VL pairs are preliminary candidates, not final IgG submissions. Independent paired-chain refolding, native-target geometry/glycan checks, developability screening, CDR novelty assessment, and constant-region assembly remain to be completed before selecting challenge submissions.
