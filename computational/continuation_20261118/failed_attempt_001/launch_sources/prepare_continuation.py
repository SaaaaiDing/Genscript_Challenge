#!/usr/bin/env python3
"""Audit the completed 500 and prepare a strictly evaluation-only positive control."""
import csv
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import rfproteina.conditions
from atomworks.io import parse
from atomworks.io.utils.atom_array_plus import as_atom_array_plus, concatenate_atom_array_plus
from atomworks.io.utils.io_utils import to_cif_file
from atomworks.io.transforms.atom_array import add_id_and_entity_annotations, add_iid_annotations
from atomworks.ml.conditions import C_SEQ
from biotite.structure import superimpose, rmsd
from rfproteina.conditions import C_CRDGRP
from rfproteina.conditions_reporting import run_pipeline_validation, check_index_group_molecule_consistency, check_index_group_res_id_gaps
from rfproteina.utils.parse_presets import RCSB_WITH_METALS, DISTILLATION, ANNOTATED_CIF
from quality_screen import ROOT, sequence

BASE = ROOT / 'continuation_20261118'
OLD = ROOT / 'campaign_local500_20261117'


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def load(path, preset):
    return as_atom_array_plus(parse(str(path), config=preset)['assemblies']['1'][0])


def main():
    records = json.loads((OLD / 'reports/all_designs.json').read_text())
    new = [r for r in records if not r['design_id'].startswith('pilot')]
    folded = [r for r in records if 'oracle' in r]
    assert len(new) == 500
    independent = (json.loads((BASE / 'reports/campaign500_audit.json').read_text())['independent_rmsd_audit']
                   if '--benchmark-only' in sys.argv else [])
    for r in ([] if '--benchmark-only' in sys.argv else folded):
        files = list((OLD / 'oracles').glob(f"set_*/results/evaluation_samples/{r['design_id']}_*/*_refolded.cif"))
        assert len(files) == 1
        a = load(r['path'], ANNOTATED_CIF)
        b = load(files[0], ANNOTATED_CIF)
        keys = lambda aa: list(zip(aa.chain_id, aa.res_id, aa.atom_name))
        lookup = {k: i for i, k in enumerate(keys(b))}
        b = b[[lookup[k] for k in keys(a)]]
        bb = np.isin(a.atom_name, ['N', 'CA', 'C', 'O'])
        target = (a.chain_id == 'T') & bb
        _, transform = superimpose(a.coord[target], b.coord[target])
        aligned = transform.apply(b.coord)
        value = float(rmsd(a.coord[bb & (a.chain_id != 'T')], aligned[bb & (a.chain_id != 'T')]))
        assert abs(value - r['oracle']['target_aligned_fv_rmsd']) < 0.001
        independent.append({'design_id': r['design_id'], 'independent_target_aligned_fv_rmsd': value})
    summary = dict(generated=500, unique_pairs=len({r['sequence_hash'] for r in new}),
                   triage_passers=sum(r['triage_pass'] for r in new),
                   per_condition={t: {'generated': sum(r['task'] == t for r in new),
                                      'triage_passers': sum(r['task'] == t and r['triage_pass'] for r in new)}
                                  for t in sorted({r['task'] for r in new})},
                   refolded_total=len(folded), refolded_new=sum(not r['design_id'].startswith('pilot') for r in folded),
                   complete_passers=sum(r.get('computational_pass', False) for r in new),
                   oracle_gate_pass_counts={k: sum(r['oracle'][k] != 0 for r in folded)
                                            for k in ['target_interface_pass', 'refold_geometry_pass']},
                   vh_vl_confidence_both_at_least_80=sum(min(r['oracle']['vh_plddt'], r['oracle']['vl_plddt']) >= 80 for r in folded),
                   pairing_iptm_at_least_60=sum(r['oracle']['pairing_iptm'] >= 60 for r in folded),
                   independent_rmsd_audit=independent)
    dump(BASE / 'reports/campaign500_audit.json', summary)
    fields = ['design_id', 'task', 'VH', 'VL', 'triage_pass', 'cdr_glycosylation_motifs', 'known_reference_similar', 'path']
    with (BASE / 'reports/all500_sequences.csv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(new)

    # Use crystal antibody coordinates only in the isolated benchmark directory.
    # Substitute the exact canonical AF target used in generation, aligned to the
    # resolved crystal receptor core. No antibody residues are redesigned here.
    raw = BASE / 'benchmark/raw/6WGL.cif'
    preset = RCSB_WITH_METALS.replace(add_missing_atoms=False, fix_arginines=False,
                                    convert_mse_to_met=False, altloc='first', return_atom_array_plus=True)
    crystal = load(raw, preset)
    target = load(ROOT / 'inputs/raw/IL4RA_AF.cif', DISTILLATION)
    target = target[(target.res_id >= 26) & (target.res_id <= 232)].copy()
    native = crystal[(crystal.chain_id == 'C') & (crystal.atom_name == 'CA')]
    moving = np.array([target.coord[(target.res_id == int(n) + 25) & (target.atom_name == 'CA')][0] for n in native.res_id])
    mismatches = [{'uniprot_residue': int(n) + 25, 'crystal': str(name),
                   'canonical': str(target.res_name[(target.res_id == int(n) + 25) & (target.atom_name == 'CA')][0])}
                  for n, name in zip(native.res_id, native.res_name)
                  if target.res_name[(target.res_id == int(n) + 25) & (target.atom_name == 'CA')][0] != name]
    assert mismatches == [{'uniprot_residue': 207, 'crystal': 'LEU', 'canonical': 'CYS'}], mismatches
    fitted, transform = superimpose(native.coord, moving)
    target.coord = transform.apply(target.coord)
    target.chain_id[:] = 'T'
    target.res_id -= 25
    refs = json.loads((ROOT / 'inputs/raw/known_il4ra_antibodies.json').read_text())['rows']
    ref = next(r for r in refs if r['Therapeutic'] == 'Dupilumab')
    parts = [target]
    for source, role, key in [('A', 'H', 'HeavySequence'), ('B', 'L', 'LightSequence')]:
        chain = crystal[(crystal.chain_id == source) & (crystal.res_id <= len(ref[key])) & ~crystal.hetero].copy()
        assert sequence(chain) == ref[key]
        chain.chain_id[:] = role
        parts.append(chain)
    control = concatenate_atom_array_plus(parts)
    control = add_id_and_entity_annotations(control, overwrite=True)
    control = add_iid_annotations(control, overwrite=True)
    C_SEQ.set_annotation_from_ground_truth(control, control.chain_id == 'T')
    groups = np.full(len(control), -1, dtype=int)
    groups[(control.chain_id == 'T') & np.isin(control.atom_name, ['N', 'CA', 'C', 'O'])] = 1
    C_CRDGRP.set_annotation(control, groups)
    checks = (run_pipeline_validation, check_index_group_molecule_consistency, check_index_group_res_id_gaps)
    errors = [x for check in checks for x in check(control)]
    assert not errors, errors
    path = BASE / 'benchmark/inputs/positive_control_6wgl.cif'
    assert not path.exists()
    to_cif_file(control, str(path), save_standard_annotations=True, chain_disambiguation='chain_iid')
    reloaded = load(path, ANNOTATED_CIF)
    errors = [x for check in checks for x in check(reloaded)]
    assert not errors, errors
    for chain in ['T', 'H', 'L']:
        assert sequence(control[control.chain_id == chain]) == sequence(reloaded[reloaded.chain_id == chain])
    assert np.allclose(control.coord, reloaded.coord, atol=0.002)
    dump(BASE / 'benchmark/provenance.json', dict(source_url='https://www.rcsb.org/structure/6WGL',
         source_sha256=hashlib.sha256(raw.read_bytes()).hexdigest(), role='evaluation_only; never a generation input',
         parser='RCSB_WITH_METALS without missing-atom addition, arginine repair or MSE conversion; first altloc',
         construct='Canonical 207-residue IL4RA AF target aligned onto 6WGL; exact reference VH/VL cropped from A/B',
         excluded='Fab constant domains, crystal target chain, waters and glycans',
         crystal_target_differences=mismatches,
         target_core_alignment_rmsd=float(rmsd(native.coord, fitted)),
         control_sha256=hashlib.sha256(path.read_bytes()).hexdigest(), validation='pipeline, molecule, gaps, sequence and coordinate round-trip'))
    best = sorted([r for r in folded if not r['design_id'].startswith('pilot') and r['oracle']['pairing_iptm'] >= 60],
                  key=lambda r: r['oracle']['target_aligned_fv_rmsd'])[:4]
    dump(BASE / 'reports/recheck_candidates.json', best)
    print(json.dumps({'summary': {k:v for k,v in summary.items() if k!='independent_rmsd_audit'},
                      'positive_control': str(path), 'recheck_ids': [r['design_id'] for r in best]}, indent=2))


if __name__ == '__main__':
    main()
