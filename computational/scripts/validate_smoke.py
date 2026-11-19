#!/usr/bin/env python3
"""Require a real paired-chain smoke output before expanding the generation batch."""
from pathlib import Path
import json
import os
os.environ['CUDA_VISIBLE_DEVICES']=''
import numpy as np
import rfproteina.conditions
from atomworks.io import parse
from atomworks.io.utils.atom_array_plus import as_atom_array_plus
from atomworks.ml.conditions import C_SEQ
from biotite.structure import get_residue_starts
from biotite.sequence import ProteinSequence
from rfproteina.utils.parse_presets import ANNOTATED_CIF
root=Path(__file__).resolve().parents[1]
files=list((root/'runs/representative_smoke/inference_samples').rglob('*.cif'))
assert len(files)==1, f'Expected one smoke structure; found {len(files)}.'
aa=as_atom_array_plus(parse(str(files[0]),config=ANNOTATED_CIF)['assemblies']['1'][0])
assert set(aa.chain_id)=={'T','H','L'}
assert np.isfinite(aa.coord).all()
sequences={}
for chain_id in ('T','H','L'):
    chain=aa[aa.chain_id==chain_id]
    sequence=''.join(ProteinSequence.convert_letter_3to1(str(r)) for r in chain.res_name[get_residue_starts(chain)])
    assert set(sequence)<=set('ACDEFGHIKLMNPQRSTVWY')
    if chain_id in {'H','L'}:
        assert 85 <= len(sequence) <= 150
        assert (~C_SEQ.mask(chain)).any(), 'No generated residues found.'
    sequences[chain_id]=sequence
canonical=json.loads((root/'inputs/raw/P24394.json').read_text())['sequence']['value'][25:232]
assert sequences['T']==canonical, 'Target identity differs from canonical IL4RA.'
(root/'provenance/smoke_validation.json').write_text(json.dumps({'status':'passed','structure':str(files[0]),'sequences':sequences,'scope':'Operational integrity only; binding and developability untested.'},indent=2)+'\n')
print('Paired-chain smoke structure passed operational checks.')
