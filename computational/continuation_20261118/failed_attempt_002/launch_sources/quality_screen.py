#!/usr/bin/env python3
"""Paired-antibody triage against the intact canonical target, with explicit provenance."""
import csv
import hashlib
import json
import re
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
from biotite.sequence import ProteinSequence
from biotite.structure.io.pdbx import CIFFile, get_structure

ROOT=Path(__file__).resolve().parents[1]
METHOD=json.loads((ROOT/'provenance/design_method.json').read_text())
SOURCE=get_structure(CIFFile.read(ROOT/'inputs/raw/IL4RA_AF.cif'),model=1,use_author_fields=True)
SOURCE=SOURCE[(SOURCE.res_id>=26)&(SOURCE.res_id<=232)]
CANONICAL=json.loads((ROOT/'inputs/raw/P24394.json').read_text())['sequence']['value'][25:232]
SPECS=json.loads((ROOT/'configs/paired_designs.json').read_text())

def fit(x,y):
    u,_,vt=np.linalg.svd((x-x.mean(0)).T@(y-y.mean(0)))
    if np.linalg.det(u@vt)<0:u[:,-1]*=-1
    rot=u@vt;return rot,y.mean(0)-x.mean(0)@rot

def sequence(a):
    return ''.join(ProteinSequence.convert_letter_3to1(str(r)) for r in a.res_name[a.atom_name=='CA'])

def analyze(path,task):
    path=Path(path);f=CIFFile.read(path);a=get_structure(f,model=1,use_author_fields=True)
    assert set(a.chain_id)=={'T','H','L'} and np.isfinite(a.coord).all()
    assert sequence(a[a.chain_id=='T'])==CANONICAL
    assert np.array_equal(f.block['atom_site']['id'].as_array(int),np.arange(1,len(a)+1))
    conditioned=f.block['atom_site']['condition_sequence_1_residue'].as_array()
    cdr=conditioned=='<M>'
    reference=np.full(a.coord.shape,np.nan)
    cat=f.block['annotation_coords-with-insertion_1_atom']
    reference[cat['idx0'].as_array(int)]=np.column_stack([cat[k].as_array(float) for k in ('val0','val1','val2')])
    lookup={(int(n),str(atom)):i for i,(n,atom) in enumerate(zip(SOURCE.res_id,SOURCE.atom_name))}
    out=np.flatnonzero((a.chain_id=='T')&np.isin(a.atom_name,['N','CA','C','O']))
    src=np.array([lookup[(int(a.res_id[i])+25,str(a.atom_name[i]))] for i in out])
    x=SOURCE.coord[src].astype(float);y=reference[out]
    assert np.isfinite(y).all()
    rot,translation=fit(x,y);native=SOURCE.coord@rot+translation
    integrity=float(np.sqrt(np.mean(np.sum((native[src]-y)**2,axis=1))))
    assert integrity<0.1, 'Saved target reference differs from canonical input.'
    raw_target_rmsd=float(np.sqrt(np.mean(np.sum((a.coord[out]-y)**2,axis=1))))
    binder=np.isin(a.chain_id,['H','L'])&~np.isin(a.element,['H','D'])
    heavy=~np.isin(SOURCE.element,['H','D'])
    distances=cdist(native[heavy],a.coord[binder])
    clash_pairs=int((distances<2.0).sum())
    nearest=distances.min(axis=1)
    native_res=SOURCE.res_id[heavy]
    contacted_target_residues=len(np.unique(native_res[nearest<=5.0]))
    hotspot_query=SPECS[task]['conditions']['C_HOT']['true'][0]['select']
    hotspots=[int(n) for n in re.search(r'\[([0-9,]+)\]',hotspot_query).group(1).split(',')]
    hotspot_distances={str(n):float(nearest[native_res==n].min()) for n in hotspots}
    hotspot_contacts=sum(d<=5 for d in hotspot_distances.values())
    ab_nearest=distances.min(axis=0);binder_indices=np.flatnonzero(binder)
    contacting=binder_indices[ab_nearest<=5]
    contact_res={(str(a.chain_id[i]),int(a.res_id[i])) for i in contacting}
    cdr_contacts={(str(a.chain_id[i]),int(a.res_id[i])) for i in contacting if cdr[i]}
    fraction=len(cdr_contacts)/len(contact_res) if contact_res else 0.0
    glyco_atoms=heavy&np.isin(SOURCE.res_id,METHOD['glycosylation_sites'])&(SOURCE.atom_name=='ND2')
    glyco_min=float(cdist(native[glyco_atoms],a.coord[binder]).min())
    framework=binder&~cdr&np.isin(a.atom_name,['N','CA','C','O'])
    ref=reference[framework];pred=a.coord[framework]
    assert np.isfinite(ref).all()
    r,t=fit(ref,pred)
    framework_rmsd=float(np.sqrt(np.mean(np.sum((ref@r+t-pred)**2,axis=1))))
    seqs={c:sequence(a[a.chain_id==c]) for c in ('H','L')}
    cdrs={}
    for chain in ('H','L'):
        ca=np.flatnonzero((a.chain_id==chain)&(a.atom_name=='CA'))
        chunks=[];chunk=''
        for i in ca:
            if cdr[i]:chunk+=ProteinSequence.convert_letter_3to1(str(a.res_name[i]))
            elif chunk:chunks.append(chunk);chunk=''
        if chunk:chunks.append(chunk)
        assert len(chunks)==3, 'Expected three regenerated CDRs per variable domain.'
        cdrs.update({chain+str(i+1):s for i,s in enumerate(chunks)})
    extra_cys=sum(s.count('C') for s in cdrs.values())
    glyco_motifs=sum(len(re.findall(r'N[^P][ST]',s)) for s in cdrs.values())
    donor_identical=sum(cdrs[c+str(i)]==METHOD['framework_details'][c]['removed_cdrs'][c+'CDR'+str(i)] for c in ('H','L') for i in (1,2,3))
    geometry_pass=(clash_pairs==0 and hotspot_contacts>=2 and contacted_target_residues>=5 and
                   len(cdr_contacts)>=3 and fraction>=0.5 and raw_target_rmsd<=2.5 and framework_rmsd<=2.5 and glyco_min>=4.0)
    sequence_pass=extra_cys==0 and donor_identical<6
    return {'path':str(path),'task':task,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'sequence_hash':hashlib.sha256((seqs['H']+'|'+seqs['L']).encode()).hexdigest(),
            'VH':seqs['H'],'VL':seqs['L'],'cdrs':cdrs,'geometry_pass':bool(geometry_pass),
            'sequence_pass':sequence_pass,'triage_pass':bool(geometry_pass and sequence_pass),
            'target_reference_rmsd':integrity,'raw_target_rmsd':raw_target_rmsd,'paired_framework_rmsd':framework_rmsd,
            'clash_pairs_below_2A':clash_pairs,'target_contact_residues':contacted_target_residues,
            'hotspot_contacts':hotspot_contacts,'hotspot_distances':hotspot_distances,
            'cdr_contact_residues':len(cdr_contacts),'cdr_contact_fraction':fraction,'glycan_site_min_distance':glyco_min,
            'extra_cdr_cysteines':extra_cys,'cdr_glycosylation_motifs':glyco_motifs,'donor_identical_cdr_count':donor_identical,
            'novelty_status':'external_reference_comparison_pending'}

def main():
    files=[ROOT/'runs/initial_designs/inference_metrics/per_example_metrics_job0.csv',ROOT/'runs/representative_smoke/inference_metrics/per_example_metrics.csv']
    results=[]
    for file in files:
        for row in csv.DictReader(file.open()):
            try:results.append(analyze(row['designed_structure_path'],row['task']))
            except Exception as e:results.append({'path':row['designed_structure_path'],'task':row['task'],'error':repr(e)})
    (ROOT/'reports/pilot_quality.json').write_text(json.dumps(results,indent=2)+'\n')
    for r in results:print(json.dumps({k:v for k,v in r.items() if k not in {'path','VH','VL','cdrs','sha256','sequence_hash'}}))

if __name__=='__main__':main()
