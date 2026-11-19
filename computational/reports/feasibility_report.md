# RFD4-Proteina feasibility and resource-isolation audit

Audit date: 2026-11-17. Scope: computational preparation only.

**Decision:** prepare paired-Fv generation and explicit paired-chain evaluation. Do not launch GPU work yet. The input template still needs a reviewed paired framework and target, and the current account has no confirmed non-preempting submission route.

No model weights were loaded, no candidate sequences were generated, no Slurm jobs were submitted, and no existing processes were stopped or reconfigured. All new project artifacts are under `computational/`. Existing Chinese notes predate the English-artifact requirement and were preserved.

## Environment and assets

Repository: `/storage/group/yjz5549/default/RFD4/RFD4-Proteina`.

Audited commit: `aaccdd2ad0720c180468df5f9cf037a713d403e5`.

| Item | Finding | Consequence |
|---|---|---|
| User-specified conda environment | `/storage/group/yjz5549/default/RFD4/env`, Python 3.12.16; no distribution metadata for RFD4, Torch, AtomWorks, Hydra, Biotite, or ANARCII | It is not currently a verified RFD4 runtime. No packages were installed or changed. |
| Existing repository environment | `.venv`, Python 3.12.10; RFD4 0.0.3, Torch 2.14.1+cu126, AtomWorks 3.0.4, Biotite 1.6.0, Hydra 1.3.8, ANARCII 2.0.8 | This is the candidate runtime used explicitly for CPU audit, not an implicit switch. |
| Released model bundle | Cached under `/storage/group/yjz5549/default/RFD4/huggingface/hub/models--nvidia--RFD4-Proteina/snapshots/a06e2509c7033b010fc43614e9fd840303588502` | Set the cache path explicitly and pin this revision. The default home cache did not contain it. |
| Native checkpoint directories | `data_path/flow_checkpoints` and `data_path/ae_checkpoints` absent | Keep the cached HF route; do not set `hf_ckpt_path=null`. |
| OpenFold3 | Existing CUDA 12 SIF and compressed converted weights; `/usr/bin/apptainer` available | Assets exist, but container execution, GPU fit, and multichain oracle behavior are not yet tested. |
| Local scheduler profile | `scripts/.user_info_local` absent | A project-specific scheduler configuration remains necessary before using the launcher. |

The HF bundle has a 1,742,660,596-byte `model.safetensors` file. OpenFold3 has a 1,369,238,320-byte compressed weight file. Presence and size checks are not integrity checks or successful model-load tests. Reuse these assets without copying or redownloading them.

Repository scripts prefer an active environment if one is supplied. Activating the user-specified conda path before running the wrapper can therefore prevent fallback to the populated `.venv`. Future execution must name and validate the intended runtime explicitly.

Existing untracked files `configs/inference_experiment/pdl1_of3_search_local.yaml` and `smoke_pdl1.yaml` were not modified. The local `.env.local` was read only for selected non-secret runtime paths.

## Resource isolation

At audit time, `p-gc-3039` had two A100 GPUs, both allocated by Slurm. The current desktop allocation `57427235` owned one GPU, 16 CPUs, and 128 GiB. Other work was active on the node, including CPU processes belonging to the same user. A visible GPU utilization of zero is only a snapshot, not an authorization to share it.

The cluster uses `select/cons_tres`, `CR_CORE_MEMORY`, and cgroup task tracking. However, `PreemptType=preempt/qos`, and QoS `normal` lists `preempt` as a preemption target. The queried account association permits `normal` only. QoS `preempt` has no listed preemption targets but is not currently shown as authorized for this account.

Therefore, the previous assumption that a normal-queue submission cannot evict other jobs is false on this cluster. Excluding the current node protects that node, but a normal-QoS job could still preempt a job on another node. Neither lowering Unix priority nor selecting an apparently idle GPU solves scheduler preemption.

Before submission, establish an authorized QoS or scheduler-enforced allocation route that cannot preempt other jobs. Do not silently choose `preempt`, switch accounts, request elevated priority, or rely on a transient idle-node snapshot. If only the current node's tasks must be protected, narrowing the requirement is a user decision; this audit retains the broader non-disruption requirement.

Provisional requests in `configs/resource_plan.json` are 1 A100, 4 CPUs, 32 GiB, and 1 hour for generation; 1 A100, 4 CPUs, 48 GiB, and 1 hour for a separate OpenFold3 pilot. These are request ceilings to test, not measured requirements or promises of fit. Generation and refolding must not overlap. No arrays, broad sweeps, or automatic OOM retries are planned.

The repository launcher hardcodes `OMP_NUM_THREADS=4` and `OPENBLAS_NUM_THREADS=4` into generated jobs. It may also add a Foldseek-dependent aggregation job and spawn refolding daemons when configured. For the pilot, disable automatic aggregation and Foldseek, remove the entire metrics group from generation, and inspect actual generated SBATCH/runtime settings before any submission. Merely exporting one thread in the calling shell would not override the launcher's later assignments.

## Paired-antibody input design

The shipped antibody tutorial and eight campaign inputs are VHH examples. They demonstrate framework preservation, CDR generation, hotspot sampling, and refinement, but are not a validated paired VH/VL recipe.

Use chain roles `T=IL4RA`, `H=VH`, and `L=VL`, with an explicit identity mapping retained through serialization and oracle renaming. Model the ECD plus Fv initially; add full-chain/Fab compatibility checks later. Do not substitute a VHH, an scFv linker, or a generic binder for the required paired chains.

Planned conditions:

1. Preserve the IL4RA target sequence and its reviewed backbone geometry in internal-coordinate group 1.
2. Preserve a documented generic human VH/VL framework and the relative VH/VL framework orientation in shared internal-coordinate group 2.
3. Remove all six source CDR segments from the input contig and replace them with newly generated segments. Do not retain an existing IL4RA binder's CDR sequences or binding pose.
4. Leave antibody-to-target rigid-body placement free. Avoid absolute coordinate conditions that inadvertently lock the input pose.
5. Add reviewed target hotspots only after examining accessibility, glycan masking, and full-domain context.
6. Use explicit antibody fold conditioning as support, not as evidence of success.

`paired_input.template.yaml` deliberately leaves source paths and the contig mandatory and unset. Its selection expressions are a draft, not a materialized or round-trip-validated CIF. No germline scaffold, CDR boundary, loop-length range, or epitope residue has been invented to make the configuration look runnable.

Target preparation still requires a canonical P24394-1 26-232 sequence, structure-to-UniProt mapping, missing-residue and mutation review, and glycan annotation. The challenge's 3BPN reference is a cytokine-receptor complex, not an antibody framework. A full ECD plus ordinary VH/VL domains is roughly 430-460 protein residues before optional glycans; residue count alone does not predict peak GPU memory.

## Evaluation findings and corrections

### Single-chain selection

`src/rfproteina/utils/chain.py:get_chain_of_interest` returns one chain, chosen by the lowest sequence-conditioned fraction and then coordinate-conditioned fraction. It cannot represent the complete paired Fv as a set. Ties may be marked ambiguous.

**Consequence:** do not use automatic COI export, single-chain Foldseek clusters, single-chain apo transforms, or implicit COI masks as the sole paired-antibody analysis.

### False-positive interface confidence

`src/rfproteina/metrics/oracles/oracle_metric.py:extract_confidence_metrics` computes `min_ipae_coi` as the minimum over a selected chain's interfaces with every other chain. A confident H-L interface can dominate poor H-T and L-T interfaces. The built-in antibody pass criterion uses this COI minimum.

The CPU check includes a synthetic counterexample: H-L error 0.5, H-T error 18, and L-T error 20. The default COI minimum reports 0.5. This is a test of score semantics, not biological data.

**Correction:** retain the H-T, L-T, and H-L pairs separately, map them by verified chain identity, and evaluate target-facing contacts. Do not require both antibody chains to contribute equally; do not treat a favorable single scalar as proof of binding. Preserve directional PAE values when available rather than only a symmetric summary.

### Refolding template leakage

The stock PPI oracle templates any protein atoms carrying absolute or internal coordinate conditions. Framework-conditioned VH/VL can therefore enter the template as well as the antigen.

**Correction:** `paired_evaluation.yaml` explicitly templates target chain T only. Framework and CDR regions are not supplied as antibody templates. This is an independent-prediction check, although it is not statistically independent experimental evidence.

### Target alignment

The stock target-aligned binder RMSD uses absolute `C_CRD` masks and a single-chain COI, while the recommended target input uses internal coordinate groups. Depending on the input and annotation handling, this can leave an unsuitable or empty alignment subset.

**Correction:** explicitly align on T backbone atoms and measure H, L, and the combined H/L backbone separately. The draft does not reuse built-in antibody pass thresholds. CDR-specific and paired-orientation checks should be added after numbering and atom correspondence are validated.

### Missing results and diversity

Oracle output presence, chain mapping, and expected metric completeness must be checked explicitly. Missing, NaN, or skipped results are not passes. Failed generations and failed evaluations remain in the attempt denominator.

Assess diversity using all six CDR sequences plus interface/pose families. Common immunoglobulin frameworks make whole-fold-only clustering insufficient. Refinement of independently generated candidates may follow a successful diverse pilot, with parent-child provenance retained.

## Validation and next gates

The CPU audit successfully checked environment metadata, asset presence, imports of Torch/AtomWorks/RFD4 condition handling, and the synthetic paired-confidence counterexample. CUDA remained uninitialized. The full check used approximately 940 MiB peak RSS and took 85.6 seconds with one permitted CPU and single-threaded numerical libraries.

The first configuration check exposed an OmegaConf merge detail: a child `out_dir: ???` did not replace an inherited `out_dir: null`. Both drafts now require `${oc.env:IL4RA_RUN_DIR}` unless an explicit `out_dir` override is supplied. The focused recheck passed for generation and evaluation in 0.53 seconds, with approximately 258 MiB peak RSS. Generation has no metrics/oracles; evaluation contains only the intended oracle and three explicit target-aligned metrics, with no stock pass criteria or custom metrics inherited.

`cpu_checks.json` preserves the original check, including the corrected configuration failure. `config_checks.json` records the successful focused recheck. `validation_summary.json` combines current statuses without rewriting the original history. Composed configurations are saved alongside them. These results do not verify GPU kernels, actual antibody input annotations, or the OpenFold3 container.

Remaining gates, in order:

1. CPU environment and configuration checks are complete; no currently unresolved failures remain in this check scope.
2. Prepare a provenance-reviewed target and generic paired framework; author and round-trip validate the actual annotated input, including all six regenerated CDRs.
3. Establish scheduler-enforced non-preemption and a separate allocation excluding `p-gc-3039`.
4. Run the repository GPU preflight inside that allocation, then generate one representative longest input at batch size 1.
5. Verify T/H/L identities, generated regions, framework pairing, output integrity, and measured GPU/host-memory/time peaks.
6. Run one separate OpenFold3 evaluation and verify explicit target-versus-Fv metrics and sequence/atom mapping.
7. Only then choose a small multi-seed pilot size based on measured cost and proceed to diversified sampling.

No performance estimate, computational pass, or generated model should be described as confirmed binding.
