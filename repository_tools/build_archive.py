"""Archive large raw artifacts without changing original calculation files."""
import csv,gzip,hashlib,io,json,tarfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'archives'
def sha(b):return hashlib.sha256(b).hexdigest()
def main():
    OUT.mkdir(exist_ok=True);files=[]
    for p in sorted((ROOT/'computational').rglob('*')):
        if not p.is_file()or p.is_symlink():continue
        if any(x in p.parts for x in ['cache','tmp','sockets','__pycache__']):continue
        if p.suffix in ['.cif','.sqlite']or p.name=='numbering.json':files.append(p)
    unique={};rows=[]
    for p in files:
        data=p.read_bytes();key=sha(data)
        if key not in unique:unique[key]=p
        rows.append(dict(original_path=p.relative_to(ROOT).as_posix(),sha256=key,bytes=len(data)))
    entries=list(unique.items());where={};parts=[]
    for n,start in enumerate(range(0,len(entries),50),1):
        path=OUT/f'raw_artifacts_{n:03}.tar.gz'
        with path.open('wb')as raw,gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0,compresslevel=6)as gz,tarfile.open(fileobj=gz,mode='w')as tar:
            for key,p in entries[start:start+50]:
                data=p.read_bytes();assert sha(data)==key
                member=f'objects/{key}{p.suffix}';info=tarfile.TarInfo(member);info.size=len(data);info.mode=0o644;info.mtime=0
                tar.addfile(info,io.BytesIO(data));where[key]=(path.name,member)
        assert path.stat().st_size<90*1024**2,'Archive must fit GitHub per-file limit.'
        parts.append(dict(path=path.relative_to(ROOT).as_posix(),bytes=path.stat().st_size,sha256=sha(path.read_bytes())))
        print(f'Wrote {path.name}: {path.stat().st_size/1024**2:.1f} MiB',flush=True)
    for r in rows:r['archive'],r['member']=where[r['sha256']]
    with (OUT/'file_map.csv').open('w',newline='')as f:
        w=csv.DictWriter(f,fieldnames=['original_path','sha256','bytes','archive','member']);w.writeheader();w.writerows(rows)
    summary=dict(original_file_count=len(rows),unique_objects=len(entries),original_bytes=sum(r['bytes']for r in rows),
                 unique_uncompressed_bytes=sum(p.stat().st_size for p in unique.values()),compressed_bytes=sum(p['bytes']for p in parts),parts=parts,
                 file_map_sha256=sha((OUT/'file_map.csv').read_bytes()),
                 exclusions='Credentials, caches, temporary files, model weights and external environments are not published. Non-CIF duplicated content-addressed artifacts are omitted; source/config/report copies remain tracked.')
    (OUT/'manifest.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items()if k!='parts'},indent=2))
if __name__=='__main__':main()
