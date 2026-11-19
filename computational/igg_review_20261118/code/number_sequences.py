"""Number all design and frozen therapeutic-reference variable sequences on CPU."""
import os
os.environ['CUDA_VISIBLE_DEVICES']=''
os.environ['OMP_NUM_THREADS']='2'
import csv, hashlib, json, time
from pathlib import Path
import torch
from anarcii import Anarcii

BASE=Path(__file__).resolve().parents[1]
ROOT=BASE.parent
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def key(s):return hashlib.sha256(s.encode()).hexdigest()
def main():
    torch.set_num_threads(2)
    records=json.loads((ROOT/'combined1000_20261118/all_designs.json').read_text())
    seqs={key(r[k]):r[k] for r in records for k in ['VH','VL']}
    refs=list(csv.DictReader((ROOT/'inputs/raw/TheraSAbDab.csv').open(encoding='utf-8-sig')))
    pairs=[]; excluded=[]
    aa=set('ACDEFGHIKLMNPQRSTVWY')
    for i,r in enumerate(refs):
        for h,l,arm in [('HeavySequence','LightSequence','primary'),('HeavySequence(ifbispec)','LightSequence(ifbispec)','secondary')]:
            hs,ls=r[h].strip(),r[l].strip()
            if hs and ls and set(hs)<=aa and set(ls)<=aa:
                seqs[key(hs)]=hs;seqs[key(ls)]=ls
                pairs.append(dict(name=r['Therapeutic'],arm=arm,target=r['Target'],VH=hs,VL=ls))
            else:excluded.append(dict(name=r['Therapeutic'],arm=arm,reason='missing_or_nonstandard_sequence'))
    dump(BASE/'reports/reference_pairs.json',pairs)
    dump(BASE/'reports/reference_excluded_inputs.json',excluded)
    dump(BASE/'reports/numbering_input_sequences.json',seqs)
    cachepath=BASE/'reports/numbering.json'
    cache=json.loads(cachepath.read_text()) if cachepath.exists() else {}
    todo=[(k,v)for k,v in seqs.items() if k not in cache]
    print(f'Numbering {len(todo)} new unique chains; {len(pairs)} reference pairs',flush=True)
    model=Anarcii(cpu=True,ncpu=2,verbose=False)
    for i in range(0,len(todo),64):
        model.number(dict(todo[i:i+64]));out=model.to_scheme('imgt')
        cache.update(out);dump(cachepath,cache)
        print(f'{len(cache)}/{len(seqs)} chains numbered',flush=True)
    dump(BASE/'reports/numbering_summary.json',dict(unique_chains=len(seqs),reference_rows=len(refs),reference_pairs=len(pairs),errors=sum(bool(v.get('error'))for v in cache.values())))
if __name__=='__main__':main()
