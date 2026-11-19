"""Independent Fv checks and a rigid, source-Fab constant-domain placement check.

This is a steric proxy for one observed Fab elbow, not a full IgG prediction.
"""
import sys,json,hashlib
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
from biotite.structure.io.pdbx import CIFFile,get_structure
BASE=Path(__file__).resolve().parents[1];ROOT=BASE.parent
sys.path.insert(0,str(ROOT/'scripts'))
from quality_screen import SOURCE,fit,sequence

def main():
    source=get_structure(CIFFile.read(ROOT/'inputs/raw/5I19.cif'),model=1,use_author_fields=True)
    spans={'H':[(25,33),(50,58),(96,108)],'L':[(26,32),(49,52),(88,97)]}
    varlen={'H':119,'L':107}
    source_fr=[];const=[]
    for ch in ['H','L']:
        ca=source[(source.chain_id==ch)&(source.atom_name=='CA')]
        fr=np.ones(varlen[ch],dtype=bool)
        for lo,hi in spans[ch]:fr[lo:hi]=False
        source_fr.extend(ca.coord[:varlen[ch]][fr])
        const.append(source[(source.chain_id==ch)&(source.res_id>int(ca.res_id[varlen[ch]-1]))&~np.isin(source.element,['H','D'])])
    const=const[0]+const[1];source_fr=np.array(source_fr)
    records=json.loads((ROOT/'combined1000_20261118/all_designs.json').read_text())
    pool=[r for r in records if r['review_eligible'] and not r['duplicate_of']]
    results={}
    for ix,r in enumerate(pool):
        path=Path(r['path']);assert hashlib.sha256(path.read_bytes()).hexdigest()==r['sha256']
        f=CIFFile.read(path);a=get_structure(f,model=1,use_author_fields=True)
        mask=f.block['atom_site']['condition_sequence_1_residue'].as_array()=='<M>'
        pred_fr=np.concatenate([a.coord[(a.chain_id==ch)&(a.atom_name=='CA')&~mask]for ch in ['H','L']])
        assert pred_fr.shape==source_fr.shape
        rot,t=fit(source_fr,pred_fr);proj=const.coord@rot+t
        elbow_fit=float(np.sqrt(np.mean(np.sum((source_fr@rot+t-pred_fr)**2,axis=1))))
        ref=np.full(a.coord.shape,np.nan);cat=f.block['annotation_coords-with-insertion_1_atom']
        ref[cat['idx0'].as_array(int)]=np.column_stack([cat[k].as_array(float)for k in ['val0','val1','val2']])
        ti=np.flatnonzero((a.chain_id=='T')&np.isin(a.atom_name,['N','CA','C','O']))
        lookup={(int(n),str(at)):i for i,(n,at)in enumerate(zip(SOURCE.res_id,SOURCE.atom_name))}
        si=[lookup[(int(a.res_id[i])+25,str(a.atom_name[i]))]for i in ti]
        tr,tt=fit(SOURCE.coord[si],ref[ti]);target=SOURCE.coord@tr+tt
        d=cdist(target[~np.isin(SOURCE.element,['H','D'])],proj)
        heavy={ch:a[(a.chain_id==ch)&~np.isin(a.element,['H','D'])] for ch in ['H','L']}
        hl=cdist(heavy['H'].coord,heavy['L'].coord)
        ss={}
        for ch in ['H','L']:
            assert sequence(a[a.chain_id==ch])==r['V'+ch]
            sg=a.coord[(a.chain_id==ch)&(a.res_name=='CYS')&(a.atom_name=='SG')]
            ss[ch]=float(np.linalg.norm(sg[0]-sg[1])) if len(sg)==2 else None
        results[r['candidate_id']]=dict(source_sha256=r['sha256'],source_fab='5I19',
            source_framework_CA_fit_rmsd_A=elbow_fit,
            source_Fab_constant_target_min_distance_A=float(d.min()),
            source_Fab_constant_target_clashes_below_2A=int((d<2).sum()),
            Fv_HL_clashes_below_2A=int((hl<2).sum()),Fv_HL_atom_contacts_below_5A=int((hl<5).sum()),
            Fv_intrachain_SG_distance_A=ss,
            full_H2L2_3D_assembly_validated=False,
            limitation='Rigid 5I19 Fab-elbow placement only; no hinge/Fc, second arm, relaxation or explicit glycans.')
        if (ix+1)%25==0:print(f'{ix+1}/{len(pool)} geometries reviewed',flush=True)
    (BASE/'reports/geometry.json').write_text(json.dumps(results,indent=2)+'\n')
    print('Completed',len(results),flush=True)
if __name__=='__main__':main()
