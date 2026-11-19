"""Verify every archived object or restore original relative artifact paths.

Usage: python3 repository_tools/restore_artifacts.py --verify-only
       python3 repository_tools/restore_artifacts.py --destination /path/to/clone
Optional --prefix computational/campaign_local500_20261117 limits restoration.
Existing files are accepted only when their SHA256 already matches.
"""
import argparse,csv,hashlib,json,tarfile
from collections import defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def digest(data):return hashlib.sha256(data).hexdigest()
def main():
    ap=argparse.ArgumentParser();g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--verify-only',action='store_true');g.add_argument('--destination',type=Path)
    ap.add_argument('--prefix',default='');args=ap.parse_args()
    manifest=json.loads((ROOT/'archives/manifest.json').read_text())
    mapping=ROOT/'archives/file_map.csv';assert digest(mapping.read_bytes())==manifest['file_map_sha256']
    rows=list(csv.DictReader(mapping.open()));byarchive=defaultdict(lambda:defaultdict(list))
    for r in rows:
        if not args.verify_only and not r['original_path'].startswith(args.prefix):continue
        rel=Path(r['original_path']);assert not rel.is_absolute()and '..'not in rel.parts
        byarchive[r['archive']][r['member']].append(r)
    count=0;restored=0
    for part in manifest['parts']:
        p=ROOT/part['path'];assert digest(p.read_bytes())==part['sha256'],str(p)
        wanted=byarchive[p.name]
        with tarfile.open(p,'r:gz')as tar:
            seen=set()
            for member in tar:
                if member.name not in wanted:continue
                assert member.isfile();data=tar.extractfile(member).read();rr=wanted[member.name]
                assert all(digest(data)==r['sha256']and len(data)==int(r['bytes'])for r in rr)
                seen.add(member.name);count+=1
                if args.destination:
                    dest=args.destination.resolve()
                    for r in rr:
                        target=dest/r['original_path'];assert dest in target.resolve().parents
                        if target.exists():assert digest(target.read_bytes())==r['sha256'],f'Refusing to overwrite changed file: {target}'
                        else:target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
                        restored+=1
            assert seen==set(wanted),(p.name,set(wanted)-seen)
    if args.verify_only:assert count==manifest['unique_objects']
    print(json.dumps(dict(status='passed',objects_verified=count,paths_restored_or_confirmed=restored)))
if __name__=='__main__':main()
