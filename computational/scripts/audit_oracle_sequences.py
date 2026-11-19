#!/usr/bin/env python3
"""Check exact sequences at the model's CIF decoder boundary, without inference."""
from pathlib import Path
import json
import os
import subprocess
import rfproteina.conditions
from atomworks.io import parse
from rfproteina.utils.parse_presets import STANDARD_ANNOTATIONS_ONLY
from rfproteina.metrics.oracles.af3_oracle import AF3Oracle
from rfproteina.metrics.oracles.oracle_metric import OracleJobInput
from quality_screen import ROOT, sequence

base=ROOT/'continuation_20261118'
records=json.loads((ROOT/'campaign_local500_20261117/reports/all_designs.json').read_text())
paths=[base/'failed_attempt_001/positive_control_6wgl.cif',base/'benchmark/inputs/positive_control_6wgl.cif']
paths += [Path(next(r['path'] for r in records if r['design_id']==did)) for did in ['b00_d013','b07_d014']]
jobs=[]
oracle=AF3Oracle(oracle_name='of3_ppi',template_query="chain_id == 'T' and is_protein",seed=20261118)
for i,path in enumerate(paths):
    aa=parse(str(path),config=STANDARD_ANNOTATIONS_ONLY)['assemblies']['1'][0]
    entries,_=oracle._build_manifest([OracleJobInput(atom_array=aa,example_id=f'sequence_audit_{i}',out_dir=base/'reports')],base/'reports')
    jobs.append({'source':str(path),'cif':entries[0]['cif'],
                 'expected':{str(c):sequence(aa[aa.chain_id==c]) for c in sorted(set(aa.chain_id))}})
input_path=base/'reports/oracle_sequence_audit_inputs.json'
input_path.write_text(json.dumps(jobs))
code='''import json,sys
from alphafold3.common import folding_input
from alphafold3.constants import chemical_components
jobs=json.load(open(sys.argv[1]));out=[]
for job in jobs:
    inp=folding_input.Input.from_mmcif(job['cif'],ccd=chemical_components.Ccd())
    actual={c.id:c.sequence for c in inp.chains}
    out.append({'source':job['source'],'expected':job['expected'],'actual':actual,'exact_match':actual==job['expected']})
with open(sys.argv[2],'w') as f:json.dump(out,f,indent=2)
for item in out:print(item['source'],item['exact_match'])
'''
subprocess.run(['apptainer','exec','--bind','/storage/group/yjz5549/default,/storage/home/sqd5856/Desktop/default',
                '--env','JAX_PLATFORMS=cpu',os.environ['AF3_IMAGE'],'/alphafold3_venv/bin/python','-c',code,
                str(input_path),str(base/'reports/oracle_sequence_boundary_audit.json')],check=True)
