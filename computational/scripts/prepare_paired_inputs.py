#!/usr/bin/env python3
"""Prepare canonical IL4RA plus germline paired Fv with all six IMGT CDRs removed."""
import os
import resource
from pathlib import Path
import json
import hashlib

os.environ['CUDA_VISIBLE_DEVICES'] = ''
os.sched_setaffinity(0, {17} & os.sched_getaffinity(0) or {min(os.sched_getaffinity(0))})
os.nice(10)
resource.setrlimit(resource.RLIMIT_AS, (8 * 1024**3, 8 * 1024**3))

import numpy as np
import torch
from anarcii import Anarcii
from biotite.sequence import ProteinSequence
from biotite.structure import get_residue_starts, sasa, superimpose, rmsd
import rfproteina.conditions
from atomworks.io import parse
from atomworks.io.utils.atom_array_plus import as_atom_array_plus, concatenate_atom_array_plus
from atomworks.io.utils.io_utils import to_cif_file
from atomworks.io.transforms.atom_array import add_id_and_entity_annotations, add_iid_annotations
from atomworks.ml.conditions import C_SEQ
from rfproteina.conditions import C_CRDGRP, C_HOT
from rfproteina.conditions_reporting import run_pipeline_validation, check_index_group_molecule_consistency, check_index_group_res_id_gaps
from rfproteina.inference_utils.condition_spec import resolve_condition_spec_inputs
from rfproteina.utils.parse_presets import DISTILLATION, RCSB_WITH_METALS, ANNOTATED_CIF

ROOT = Path(__file__).resolve().parents[1]
torch.set_num_threads(1)

def write_json(path, data):
    path.write_text(json.dumps(data, indent=2) + '\n')

def load(name, preset):
    config = preset.replace(add_missing_atoms=False, fix_arginines=False,
                            convert_mse_to_met=False, altloc='first', return_atom_array_plus=True)
    parsed = parse(str(ROOT / 'inputs/raw' / name), config=config)
    return as_atom_array_plus(parsed['assemblies']['1'][0])

def sequence(aa):
    indices = get_residue_starts(aa)
    return ''.join(ProteinSequence.convert_letter_3to1(str(x)) for x in aa.res_name[indices])

def main():
    target = load('IL4RA_AF.cif', DISTILLATION)
    target = target[(target.res_id >= 26) & (target.res_id <= 232)].copy()
    target.chain_id[:] = 'T'
    canonical = json.loads((ROOT / 'inputs/raw/P24394.json').read_text())
    expected = canonical['sequence']['value'][25:232]
    assert sequence(target) == expected and len(get_residue_starts(target)) == 207
    assert np.isfinite(target.coord).all()
    (ROOT / 'inputs/IL4RA_ECD.fasta').write_text('>IL4RA_P24394_1_26_232\n' + expected + '\n')

    # Compare the complete canonical prediction with the experimentally observed core.
    crystal = load('3BPN.cif', RCSB_WITH_METALS)
    crystal = crystal[(crystal.chain_id == 'B') & (crystal.res_id >= 5) & (crystal.atom_name == 'CA')]
    native_ids = crystal.res_id + 22
    predicted = np.array([target.coord[(target.res_id == n) & (target.atom_name == 'CA')][0] for n in native_ids])
    fitted, _ = superimpose(crystal.coord, predicted)
    comparison_rmsd = float(rmsd(crystal.coord, fitted))
    assert comparison_rmsd < 3.0, f'Unexpected target core divergence: {comparison_rmsd}'

    fab = load('5I19.cif', RCSB_WITH_METALS)
    numberer = Anarcii(cpu=True, ncpu=1, verbose=False)
    seqs = {str(c): sequence(fab[fab.chain_id == c]) for c in np.unique(fab.chain_id)}
    numberer.number(seqs)
    numbered = numberer.to_scheme('imgt')
    parts, chain_data, tokens = [target], {}, ['T26-232']
    cdr_names = ['HCDR1','HCDR2','HCDR3','LCDR1','LCDR2','LCDR3']
    for source_chain, record in numbered.items():
        assert not record.get('error'), record
        role = 'H' if record['chain_type'] == 'H' else 'L'
        assert role not in chain_data
        chain = fab[fab.chain_id == source_chain].copy()
        starts = get_residue_starts(chain)
        entries = [(number, letter) for number, letter in record['numbering'] if letter != '-']
        offset = int(record['query_start'])
        assert seqs[source_chain][offset:offset+len(entries)] == ''.join(letter for _,letter in entries)
        ids = chain.res_id[starts[offset:offset+len(entries)]]
        chain = chain[np.isin(chain.res_id, ids)].copy()
        chain.chain_id[:] = role
        # Source coordinates are used only for the generic paired framework.
        parts.append(chain)
        regions = []
        for index, ((number, insertion), letter) in enumerate(entries):
            region = next((f'{role}CDR{i}' for i,(lo,hi) in enumerate([(27,38),(56,65),(105,117)],1) if lo <= number <= hi), 'FR')
            regions.append(region)
        chain_tokens, donor_cdrs = [], {}
        pos = 0
        while pos < len(ids):
            end = pos+1
            while end < len(ids) and regions[end] == regions[pos]: end += 1
            region = regions[pos]
            if region == 'FR':
                chain_tokens.append(f'{role}{ids[pos]}-{ids[end-1]}')
            else:
                donor_cdrs[region] = ''.join(letter for _,letter in entries[pos:end])
                length = end-pos
                chain_tokens.append(f'{region}:{length}-{length}')
            pos = end
        assert len(donor_cdrs) == 3
        chain_data[role] = {'source_chain':source_chain,'source_variable_sequence':sequence(chain),
                            'removed_cdrs':donor_cdrs,'tokens':chain_tokens}
        tokens += ['/0'] + chain_tokens

    assert set(chain_data) == {'H','L'}
    combined = concatenate_atom_array_plus(parts)
    combined = add_id_and_entity_annotations(combined, overwrite=True)
    combined = add_iid_annotations(combined, overwrite=True)
    to_cif_file(combined, str(ROOT/'inputs/paired_source.cif'), save_standard_annotations=True, chain_disambiguation='chain_iid')

    # Select exposed patches away from known glycosylation sites and the Fc-proximal terminus.
    glyco = [f['location']['start']['value'] for f in canonical['features']
             if f['type']=='Glycosylation' and 26 <= f['location']['start']['value'] <= 232]
    ca = target[target.atom_name == 'CA']
    glyco_coords = np.array([ca.coord[ca.res_id == n][0] for n in glyco])
    surface = sasa(target, point_number=200)
    candidates = []
    for i,n in enumerate(ca.res_id):
        exposure = float(np.nansum(surface[target.res_id == n]))
        glyco_distance = float(np.linalg.norm(glyco_coords-ca.coord[i],axis=1).min())
        confidence = float(ca.b_factor[i])
        if exposure > 60 and glyco_distance > 12 and confidence >= 85 and 35 <= n <= 215:
            candidates.append((int(n),i,exposure,glyco_distance,confidence))
    assert len(candidates) >= 4, 'Insufficient accessible nonglycan candidate residues.'
    patches = []
    for center in sorted(candidates, key=lambda v:-v[2]):
        if patches and min(np.linalg.norm(ca.coord[center[1]]-p['center']) for p in patches) < 18: continue
        neighbours = [c for c in candidates if np.linalg.norm(ca.coord[c[1]]-ca.coord[center[1]]) < 10]
        if len(neighbours)<2: continue
        neighbours.sort(key=lambda v:-v[2])
        patches.append({'center':ca.coord[center[1]],'hotspots':[x[0] for x in neighbours[:3]],'evidence':neighbours[:3]})
        if len(patches)==2:break
    assert len(patches)==2

    framework = "(sel('H') or sel('L')) and not sel('*/["+','.join(cdr_names)+"]')"
    spec = {}
    for patch_index,patch in enumerate(patches,1):
        for h3_length in (12,16):
            contig = ', '.join(f'HCDR3:{h3_length}-{h3_length}' if t.startswith('HCDR3:') else t for t in tokens)
            name=f'il4ra_p{patch_index}_h3_{h3_length}'
            spec[name]={'input':str(ROOT/'inputs/paired_source.cif'),'input_as_annotated_cif':True,'contig':contig,
                        'conditions':{'C_SEQ':{True:[{'select':"sel('T')"},{'select':framework}]},
                                      'C_CRDGRP':{1:[{'select':"sel('T') and is_protein_backbone"}],
                                                  2:[{'select':framework+' and is_protein_backbone'}]},
                                      'C_CATH':{'2.60.40':[{'select':"sel('H') or sel('L')"}]},
                                      'C_HOT':{True:[{'select':"sel('T/*/["+','.join(map(str,patch['hotspots']))+"]/CA')"}]}}}
    spec_path=ROOT/'configs/paired_designs.json'
    write_json(spec_path,spec)
    outputs=list(resolve_condition_spec_inputs(spec_path,out_dir=ROOT/'inputs/annotated',validate=True))
    validation=[]
    for path in outputs:
        loaded=as_atom_array_plus(parse(str(path),config=ANNOTATED_CIF)['assemblies']['1'][0])
        violations=[x for check in (run_pipeline_validation,check_index_group_molecule_consistency,check_index_group_res_id_gaps) for x in check(loaded)]
        assert not violations, violations
        assert set(loaded.chain_id)=={'T','H','L'}
        cdr=np.isin(loaded.res_name,cdr_names)
        assert len(set(loaded.res_name[cdr]))==6
        assert not C_SEQ.mask(loaded)[cdr].any()
        assert not C_CRDGRP.is_internal(loaded)[cdr].any()
        validation.append({'path':str(path),'sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest(),'checks':'pipeline, molecule, gaps, all six CDR conditions','hotspot_atoms':int(C_HOT.mask(loaded).sum())})
    write_json(ROOT/'provenance/paired_input_validation.json',validation)
    write_json(ROOT/'provenance/design_method.json',{
        'target':'P24394-1 residues 26-232, canonical sequence, AlphaFold v6 backbone',
        'experimental_core_rmsd_angstrom':comparison_rmsd,
        'framework':'5I19 human germline IGHV3-23/IGKV1-39; retain framework only',
        'cdr_numbering':'IMGT; regenerate positions 27-38, 56-65, 105-117 in both chains',
        'framework_details':chain_data,'glycosylation_sites':glyco,
        'patches':[{'hotspots':p['hotspots'],'evidence':p['evidence']} for p in patches],
        'limits':['Glycan-distance filtering is a heuristic, not explicit glycan modeling.','No experimental binding claim.','Framework and CDR novelty require final reference comparisons.']})
    print('PAIRED_INPUTS_READY',len(outputs),flush=True)

if __name__=='__main__':main()
