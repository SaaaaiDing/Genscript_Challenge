"""Cross-format checks, chain connectivity accounting and frozen-source checks."""
import json,sys,hashlib,zipfile,re
from pathlib import Path
from collections import Counter
from finalize_review import BASE,ROOT,sheet_values,sim,cdrseq,AA
def main():
    c=json.loads((BASE/'reports/constant_regions.json').read_text())
    allrows=[]
    for name in ['priority20','reserve20']:
        rows=json.loads((BASE/'exports'/f'{name}.json').read_text());assert len(rows)==20
        allrows.extend(rows)
        seqs={};key=None
        for line in (BASE/'exports'/f'{name}_full_HC_LC.fasta').read_text().splitlines():
            if line.startswith('>'):key='|'.join(line[1:].split('|')[:2]);assert key not in seqs;seqs[key]=''
            else:seqs[key]+=line.strip()
        assert len(seqs)==40
        with zipfile.ZipFile(BASE/'exports'/f'{name}_official_template.xlsx')as z:
            assert z.testzip()is None
            v=sheet_values(z,'xl/worksheets/sheet1.xml')
            assert b'fullCalcOnLoad="1"'in z.read('xl/workbook.xml')
            assert 'IgG Example'not in ''.join(v.values())
            for i,r in enumerate(rows,2):
                a=r['IgG'];ident=r['candidate_id']
                assert hashlib.sha256(Path(r['path']).read_bytes()).hexdigest()==r['sha256']
                assert seqs[ident+'|HC']==a['HC']==v[f'C{i}']+v[f'E{i}'].rstrip('*')
                assert seqs[ident+'|LC']==a['LC']==v[f'F{i}']+v[f'H{i}'].rstrip('*')
                assert a['HC']==r['VH']+c['CH']and a['LC']==r['VL']+c['CL']
                assert set(a['HC']+a['LC'])<=AA
                assert a['HC_length']in [449,453]and a['LC_length']==214
                assert len(re.findall(r'N[^P][ST]',a['HC']))==1 and not re.search(r'N[^P][ST]',a['LC'])
                assert r['review_v2_eligible']and not r['additional_review_reasons']
                used=[]
                for kind,bonds in [('H',a['HC_intrachain_disulfide_positions']),('L',a['LC_intrachain_disulfide_positions'])]:
                    for copy in ['1','2']:
                        used.extend((kind+copy,p)for bond in bonds for p in bond)
                for bond in a['H2L2_interchain_disulfide_positions']:
                    used.extend((bond[f'chain{n}'],bond[f'position{n}'])for n in [1,2])
                expected={(ch,p)for ch in ['H1','H2','L1','L2']for p,aa in enumerate(a['HC']if ch[0]=='H'else a['LC'],1)if aa=='C'}
                assert len(used)==32 and len(set(used))==32 and set(used)==expected
    assert len({r['sequence_hash']for r in allrows})==40
    maxsim=max(sim(cdrseq(a),cdrseq(b))for i,a in enumerate(allrows)for b in allrows[i+1:])
    assert maxsim<.8
    assert all(r['sequence_liabilities']['chemical_motif_count']==0 for r in allrows[:20])
    result=dict(status='passed',candidate_count=40,workbooks=2,FASTA_records=80,
                full_H2L2_expected_disulfides_per_candidate=16,all_cysteines_accounted_once=True,
                max_selected_six_CDR_similarity=maxsim,source_CIF_hashes_verified=40,
                independent_CSV_JSON_workbook_FASTA_consistency=True,
                actual_full_IgG_3D_validation=False)
    # Also independently parse each CSV and compare the sequence-bearing columns.
    import csv
    for name in ['priority20','reserve20']:
        csvrows=list(csv.DictReader((BASE/'exports'/f'{name}_review.csv').open()))
        js=json.loads((BASE/'exports'/f'{name}.json').read_text());assert len(csvrows)==len(js)
        for a,b in zip(csvrows,js):
            assert a['candidate_id']==b['candidate_id']
            for k in ['HC','LC']:assert a[k]==b['IgG'][k]
            for k in ['VH','VL']:assert a[k]==b[k]
    (BASE/'reports/validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
