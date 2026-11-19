# Challenge continuation, 2026-11-18

The original 500-design objective is complete. This continuation follows the
user's instruction to keep working toward novel paired IgG IL4RA binders.

## Evidence from the completed run

| Condition | Generated | Initial triage pass |
| --- | ---: | ---: |
| Patch 1, H3 length 12 | 125 | 69 |
| Patch 1, H3 length 16 | 125 | 73 |
| Patch 2, H3 length 12 | 125 | 42 |
| Patch 2, H3 length 16 | 125 | 35 |
| Total | 500 | 219 |

All 500 VH/VL sequence pairs are unique. The previous controller's 223 initial
passers include four imported pilot passers; 219 is the new-design denominator.
It refolded 24 new designs and two pilots. None met the complete eligibility
criteria. Of the 26, 23 had both VH/VL pLDDT >=80 and 22 had H-L ipTM >=60;
none met target-interface or refold-geometry requirements. Strong pairing alone
does not establish target binding. The remaining unrefolded candidates have not
been demonstrated to fail independent refolding.

`reports/campaign500_audit.json` independently verifies every target-aligned Fv
RMSD using Biotite superposition; all agree within 0.001 Angstrom with the
original implementation. `reports/all500_sequences.csv` retains all sequences
and initial screening labels. Neither file is a final submission shortlist.

## Next computation

1. Refold an isolated positive control constructed from the known 6WGL Fv and
   the canonical 207-residue target used by the design workflow, with two seeds.
2. Recheck four original de novo designs with the smallest pose deviations among
   those maintaining H-L pairing, using a second seed. This tests reproducibility,
   not evidence of measured binding.
3. Only if both positive-control runs recover the interface and pose under the
   original thresholds, generate eight partial-diffusion variants of each of
   those four de novo designs (32 total), using half of the existing integration
   trajectory. This is a bounded hypothesis test of local sequence/backbone
   refinement. The same target, framework conditions, checkpoint, and recycle
   preset are retained; only the initialization/trajectory budget changes.
4. Independently refold up to four eligible unique variants per parent; retain
   all generated failures in the denominator and apply the unchanged complete
   candidate criteria. At most 20 diverse computational review candidates can
   be exported, possibly zero. Review the outcome before any further scaling.

The known antibody is evaluation-only. Its files are confined to `benchmark/`
and diagnostic oracle inputs. The refinement code asserts that every parent
comes from the completed de novo campaign and matches its recorded SHA-256.
The 6WGL receptor has L at canonical position 207, whereas the challenge target
has C. This benchmark deliberately uses the canonical target, not the crystal
mutation. Source details, structural alignment, cropping, removed glycans, and
round-trip validation are retained in `benchmark/provenance.json`.

Source: [RCSB PDB 6WGL](https://www.rcsb.org/structure/6WGL).

## Operation and remaining challenge work

The controller runs on p-gc-3039 under existing allocation 57427235, with one
A100, two CPUs, 48 GiB step memory, and a four-hour bound. It serializes model
generation and refolding. `status.json`, `controller.log`, `oracles/*/summary.json`,
and `reports/refinement_designs.json` record progress and outcomes. Create `STOP`
here to prevent the next refinement batch. Do not restart over existing state.

Any review candidates still need broader CDR novelty checking, developability,
glycan/Fc accessibility and full-IgG context review, full HC/LC assembly, and
official template validation. This workflow does not place an order or submit
challenge candidates. It does not infer experimental affinity from model scores.
