"""Auditable sequence screening, exact vendor-constant assembly and workbook export."""
import collections,csv,hashlib,json,re,sys,zipfile
from pathlib import Path
from xml.etree import ElementTree as ET
from xml.sax.saxutils import escape
BASE=Path(__file__).resolve().parents[1];ROOT=BASE.parent
AA=set('ACDEFGHIKLMNPQRSTVWY');ORDER=('H1','H2','H3','L1','L2','L3')
NS={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def seqkey(s):return hashlib.sha256(s.encode()).hexdigest()
def cdrseq(r):return '|'.join(r['cdrs'][k]for k in ORDER)
def distance(a,b):
    # Myers bit-vector Levenshtein; tested against ordinary dynamic programming.
    if not a:return len(b)
    peq={}
    for i,c in enumerate(a):peq[c]=peq.get(c,0)|(1<<i)
    pv=~0;mv=0;score=len(a);last=1<<(len(a)-1)
    for c in b:
        eq=peq.get(c,0);xv=eq|mv;xh=(((eq&pv)+pv)^pv)|eq
        ph=mv|~(xh|pv);mh=pv&xh
        if ph&last:score+=1
        if mh&last:score-=1
        ph=(ph<<1)|1;mh<<=1;pv=mh|~(xv|ph);mv=ph&xv
    return score
def sim(a,b):return 1-distance(a,b)/max(len(a),len(b),1)
def sheet_values(z,name):
    ss=[''.join(n.itertext()) for n in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si',NS)]
    out={}
    for c in ET.fromstring(z.read(name)).findall('.//m:sheetData/m:row/m:c',NS):
        v=c.find('m:v',NS)
        if c.attrib.get('t')=='inlineStr':s=''.join(c.find('m:is',NS).itertext())
        elif v is None or v.text is None:s=''
        else:s=ss[int(v.text)] if c.attrib.get('t')=='s' else v.text
        out[c.attrib['r']]=s
    return out
def constants():
    p=BASE/'raw/official_template.xlsx'
    with zipfile.ZipFile(p)as z:v=sheet_values(z,'xl/worksheets/sheet2.xml')
    assert v['A3']=='Human IgG1' and v['D3']=='Human Ig Kappa'
    ch=v['B3'].rstrip('*');cl=v['E3'].rstrip('*')
    assert v['B3']==ch+'**' and v['E3']==cl+'**'
    assert ch==''.join((BASE/'raw/P01857-1.fasta').read_text().splitlines()[1:])
    assert cl==json.loads((BASE/'raw/P01834.json').read_text())['sequence']['value']
    from biotite.sequence import NucleotideSequence
    assert str(NucleotideSequence(v['C3']).translate(complete=True))==v['B3']
    assert str(NucleotideSequence(v['F3']).translate(complete=True))==v['E3']
    assert len(ch)==330 and len(cl)==107 and set(ch+cl)<=AA
    x=dict(CH_Type=v['A3'],CL_Type=v['D3'],CH=ch,CL=cl,CH_template=v['B3'],CL_template=v['E3'],
           template_url='https://www.genscript.com/gsfiles/techfiles/TurboCHO_IgG_AA_Seperate_Region.xlsx',
           template_sha256=sha(p),CH_source_cell='Constant Region Seq!B3',CL_source_cell='Constant Region Seq!E3',
           independent_identity_checks=['UniProt P01857-1 secreted isoform','UniProt P01834','Translation of vendor constant DNA'],
           stop_handling='Keep vendor ** in workbook constant formula/cache; remove terminal ** from mature-protein FASTA.',
           signal_peptide='None; mature amino acid sequences, vendor supplies expression signal.',
           stoichiometry='H2L2: two identical heavy chains and two identical kappa light chains; no linker or tag')
    dump(BASE/'reports/constant_regions.json',x);return x
def numbered(seq,cache):
    n=cache[seqkey(seq)]
    if n.get('error') or not n.get('numbering'):return dict(error=n.get('error')or'missing_numbering')
    entries=[(p,a)for p,a in n['numbering']if a!='-']
    start=int(n['query_start']);letters=''.join(a for _,a in entries)
    assert seq[start:start+len(letters)]==letters
    cdrs={str(i):''.join(a for (p,_),a in entries if lo<=p<=hi)for i,(lo,hi)in enumerate([(27,38),(56,65),(105,117)],1)}
    framework=''.join(a for (p,_),a in entries if not any(lo<=p<=hi for lo,hi in [(27,38),(56,65),(105,117)]))
    cys=[i+start+1 for i,((p,_),a)in enumerate(entries)if p in [23,104]and a=='C']
    return dict(chain_type=n['chain_type'],start=start,end=start+len(entries),cdrs=cdrs,framework=framework,conserved_C_positions=cys)
def liabilities(r):
    hits=[];runs=[];charge=0
    for k,s in r['cdrs'].items():
        charge+=sum(s.count(c)for c in 'KR')-sum(s.count(c)for c in 'DE')
        for pattern,name in [(r'N[GS]','deamidation_motif'),(r'D[GS]','isomerization_motif'),(r'DP','cleavage_motif'),(r'[MW]','oxidation_susceptible_residue')]:
            for m in re.finditer('(?=('+pattern+'))',s):hits.append(dict(cdr=k,position=m.start()+1,motif=m.group(1),kind=name))
        runs.extend([dict(cdr=k,position=m.start()+1,sequence=m.group())for m in re.finditer(r'[AVILMFWY]{5,}',s)])
    chem=sum(h['kind']in ['deamidation_motif','isomerization_motif']for h in hits)
    return dict(motif_hits=hits,chemical_motif_count=chem,hydrophobic_runs_ge5=runs,CDR_KR_minus_DE=charge,
                risk_tier=0 if chem==0 else 1,
                interpretation='Sequence flags are prioritization heuristics, not measured degradation, aggregation or expression.')
def broad_novelty(r,refs):
    query=cdrseq(r);scores=[];flags=[]
    for ref in refs:
        six=sim(query,ref['six']);h3=sim(r['cdrs']['H3'],ref['cdrs']['H3'])
        scores.append((six,h3,ref['name'],ref['arm']))
        if six>=.8 and h3>=.8:flags.append(dict(name=ref['name'],arm=ref['arm'],six_cdr_similarity=six,H3_similarity=h3))
    best=max(scores)
    return dict(reference_pairs=len(refs),flagged_matches=flags,nearest_name=best[2],nearest_arm=best[3],
                max_six_CDR_similarity=best[0],H3_similarity_at_nearest=best[1],
                scope='Frozen local TheraSAbDab CSV (clinical-status column Feb 2025), paired sequences with successful IMGT numbering; not exhaustive novelty or patent clearance.')
def select(pool,limit):
    selected=[]
    # Exhaust motif-free candidates before admitting motif-bearing candidates.
    for tier in [0,1]:
        groups={t:sorted([r for r in pool if r['task']==t and r['sequence_liabilities']['risk_tier']==tier],key=lambda r:(
            -r['cdr_contact_fraction'],-r['cdr_contact_residues'],r['sequence_liabilities']['chemical_motif_count'],
            r['paired_framework_rmsd'],-r['glycan_site_min_distance'],r['candidate_id']))for t in sorted({r['task']for r in pool})}
        while len(selected)<limit:
            progress=False
            for t,rows in groups.items():
                while rows:
                    r=rows.pop(0)
                    if all(sim(cdrseq(r),cdrseq(s))<.8 for s in selected):
                        selected.append(r);progress=True;break
                if len(selected)>=limit:break
            if not progress:break
        if len(selected)>=limit:break
    return selected
def workbook(rows,dest,c):
    with zipfile.ZipFile(BASE/'raw/official_template.xlsx')as zin,zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED)as zout:
        for info in zin.infolist():
            data=zin.read(info.filename)
            if info.filename=='xl/worksheets/sheet1.xml':
                s=data.decode()
                for i,r in enumerate(rows,2):
                    cells=[]
                    for col,style,value in [('B',9,r['candidate_id']),('C',9,r['VH']),('D',10,c['CH_Type']),('F',9,r['VL']),('G',16,c['CL_Type'])]:
                        cells.append((col,f'<c r="{col}{i}" s="{style}" t="inlineStr"><is><t>{escape(value)}</t></is></c>'))
                    cells.append(('A',f'<c r="A{i}" s="15"><f>IF(B{i}="","",COUNTA($B$2:B{i}))</f><v>{i-1}</v></c>'))
                    for col,typ,ran,style,value in [('E','D','A$3:$C$6',1,c['CH_template']),('H','G','D$3:$F$6',11,c['CL_template'])]:
                        formula=f'IF({typ}{i}="","",VLOOKUP({typ}{i},\'Constant Region Seq\'!${ran},2,0))'
                        cells.append((col,f'<c r="{col}{i}" s="{style}" t="str"><f>{escape(formula)}</f><v>{value}</v></c>'))
                    replacement=''.join(v for _,v in sorted(cells))
                    pattern=rf'(<row\b[^>]*\br="{i}"[^>]*>).*?(</row>)'
                    s,count=re.subn(pattern,lambda m:m.group(1)+replacement+m.group(2),s,flags=re.S)
                    assert count==1,(i,count)
                data=s.encode()
            elif info.filename=='xl/workbook.xml':
                data=re.sub(rb'<calcPr\b[^>]*/>',b'<calcPr calcId="191029" fullCalcOnLoad="1" forceFullCalc="1"/>',data)
            zout.writestr(info,data)
    with zipfile.ZipFile(dest)as z:
        v=sheet_values(z,'xl/worksheets/sheet1.xml')
        assert z.read('xl/worksheets/sheet2.xml')==zipfile.ZipFile(BASE/'raw/official_template.xlsx').read('xl/worksheets/sheet2.xml')
        for i,r in enumerate(rows,2):
            assert [v[f'{k}{i}']for k in 'BCDFG']==[r['candidate_id'],r['VH'],c['CH_Type'],r['VL'],c['CL_Type']]
            assert v[f'E{i}']==c['CH_template']and v[f'H{i}']==c['CL_template']
        assert sum(bool(value)for k,value in v.items()if re.fullmatch('B[0-9]+',k)and k!='B1')==len(rows)
def assembly(r,c):
    h,l=r['VH']+c['CH'],r['VL']+c['CL'];nh,nl=len(r['VH']),len(r['VL'])
    assert set(h+l)<=AA and h[:nh]==r['VH']and l[:nl]==r['VL']
    assert h[nh:]==c['CH']and l[nl:]==c['CL']
    hdis=[r['numbering']['H']['conserved_C_positions']]+[[nh+a,nh+b]for a,b in [(27,83),(144,204),(250,308)]]
    ldis=[r['numbering']['L']['conserved_C_positions'],[nl+27,nl+87]]
    inter=[dict(chain1='H1',position1=nh+103,chain2='L1',position2=nl+107),
           dict(chain1='H2',position1=nh+103,chain2='L2',position2=nl+107),
           dict(chain1='H1',position1=nh+109,chain2='H2',position2=nh+109),
           dict(chain1='H1',position1=nh+112,chain2='H2',position2=nh+112)]
    for seq,bonds in [(h,hdis),(l,ldis)]:
        for a,b in bonds:assert seq[a-1]==seq[b-1]=='C'
    for bond in inter:
        for n in [1,2]:assert (h if bond[f'chain{n}'].startswith('H')else l)[bond[f'position{n}']-1]=='C'
    assert h.count('C')==11 and l.count('C')==5
    assert h[nh+179:nh+182]=='NST'
    return dict(HC=h,LC=l,HC_length=len(h),LC_length=len(l),VH_CH_junction=r['VH'][-10:]+'|'+c['CH'][:10],
                VL_CL_junction=r['VL'][-10:]+'|'+c['CL'][:10],HC_intrachain_disulfide_positions=hdis,
                LC_intrachain_disulfide_positions=ldis,H2L2_interchain_disulfide_positions=inter,
                Fc_N_glycosylation_N_position=nh+180,sequence_assembly_pass=True,
                disulfide_interpretation='Expected connectivity from conserved sequence; not proof of actual oxidation or assembly.',
                full_H2L2_3D_assembly_validated=False)
def main():
    c=constants();cache=json.loads((BASE/'reports/numbering.json').read_text())
    geometry=json.loads((BASE/'reports/geometry.json').read_text())
    records=json.loads((ROOT/'combined1000_20261118/all_designs.json').read_text())
    method=json.loads((ROOT/'provenance/design_method.json').read_text())
    # Derive the exact expected framework from the original retained fragments.
    expected={}
    for ch in ['H','L']:
        s=method['framework_details'][ch]['source_variable_sequence'];spans={'H':[(25,33),(50,58),(96,108)],'L':[(26,32),(49,52),(88,97)]}[ch]
        expected[ch]=''.join(a for i,a in enumerate(s)if not any(lo<=i<hi for lo,hi in spans))
    refs=[];excluded=[]
    for r in json.loads((BASE/'reports/reference_pairs.json').read_text()):
        h,l=numbered(r['VH'],cache),numbered(r['VL'],cache)
        if h.get('error')or l.get('error')or h.get('chain_type')!='H'or l.get('chain_type')not in ['K','L']:
            excluded.append(dict(name=r['name'],arm=r['arm'],H=h,L=l));continue
        r['cdrs']={ch+str(i):n['cdrs'][str(i)]for ch,n in [('H',h),('L',l)]for i in [1,2,3]}
        r['six']=cdrseq(r);refs.append(r)
    dump(BASE/'reports/reference_numbering_exclusions.json',excluded)
    dump(BASE/'reports/numbered_reference_cdrs.json',refs)
    for i,r in enumerate(records):
        reasons=[];r['numbering']={}
        for ch in ['H','L']:
            s=r['V'+ch];n=numbered(s,cache);r['numbering'][ch]=n
            if n.get('error'):reasons.append(ch+'_numbering_error');continue
            checks={ch+'_wrong_chain_type':n['chain_type']==('H'if ch=='H'else'K'),
                    ch+'_incomplete_variable_domain':n['start']==0 and n['end']==len(s),
                    ch+'_framework_changed':n['framework']==expected[ch],
                    ch+'_conserved_cysteine':len(n['conserved_C_positions'])==2 and s.count('C')==2,
                    ch+'_CDR_boundary_mismatch':all(n['cdrs'][str(k)]==r['cdrs'][ch+str(k)]for k in [1,2,3]),
                    ch+'_variable_glycosylation':not re.search(r'N[^P][ST]',s)}
            reasons.extend(k for k,v in checks.items()if not v)
        r['sequence_liabilities']=liabilities(r)
        if r['sequence_liabilities']['hydrophobic_runs_ge5']:reasons.append('CDR_hydrophobic_run_ge5')
        if abs(r['sequence_liabilities']['CDR_KR_minus_DE'])>6:reasons.append('CDR_absolute_KR_minus_DE_gt6')
        if r['candidate_id']in geometry:
            g=geometry[r['candidate_id']];r['assembly_geometry']=g
            if g['Fv_HL_clashes_below_2A']:reasons.append('Fv_HL_clash')
            if g['source_Fab_constant_target_clashes_below_2A']:reasons.append('source_Fab_constant_target_clash')
            if any(x is None or not 1.5<=x<=3 for x in g['Fv_intrachain_SG_distance_A'].values()):reasons.append('Fv_disulfide_geometry')
            r['broader_novelty']=broad_novelty(r,refs)
            if r['broader_novelty']['flagged_matches']:reasons.append('broader_reference_similarity')
        r['additional_review_reasons']=reasons
        r['review_v2_eligible']=r['review_eligible']and not r['duplicate_of']and not reasons
        if (i+1)%100==0:print('Reviewed',i+1,flush=True)
    pool=[r for r in records if r['review_v2_eligible']]
    selected=select(pool,40);assert len(selected)>=20,'Not enough candidates: inspect review before further relaxation.'
    for i,r in enumerate(selected):r['review_v2_rank']=i+1;r['review_v2_queue']='priority'if i<20 else'reserve';r['IgG']=assembly(r,c)
    for name,rows in [('priority20',selected[:20]),('reserve20',selected[20:])]:
        dump(BASE/'exports'/f'{name}.json',rows)
        workbook(rows,BASE/'exports'/f'{name}_official_template.xlsx',c)
        (BASE/'exports'/f'{name}_full_HC_LC.fasta').write_text(''.join(f'>{r["candidate_id"]}|{ch}|mature|Human_IgG1_kappa\n{r["IgG"][ch]}\n'for r in rows for ch in ['HC','LC']))
        with (BASE/'exports'/f'{name}_review.csv').open('w')as f:
            fields=['candidate_id','task','risk_tier','chemical_motif_count','cdr_contact_fraction','cdr_contact_residues','paired_framework_rmsd','glycan_site_min_distance','oracle_status','original_full_gate_pass','HC_length','LC_length','VH','VL','HC','LC']
            w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader()
            for r in rows:w.writerow(dict(r,**r['sequence_liabilities'],**r['IgG']))
    dump(BASE/'reports/all1000_reviewed.json',records)
    old={r['candidate_id']for r in json.loads((ROOT/'combined1000_20261118/priority_review_20.json').read_text())}
    summary=dict(status='completed',total_generated=len(records),unique_pairs=len({r['sequence_hash']for r in records}),
                 original_unique_review_pool=len(geometry),review_v2_pool=len(pool),priority_count=len(selected[:20]),reserve_count=len(selected[20:]),
                 reference_source_rows=1133,reference_pairs_screened=len(refs),reference_pairs_numbering_excluded=len(excluded),
                 numbering_errors_in_designs=sum(any(n.get('error')for n in r['numbering'].values())for r in records),
                 additional_reasons_in_original_pool=dict(collections.Counter(x for r in records if r['candidate_id']in geometry for x in r['additional_review_reasons'])),
                 pool_risk_tiers=dict(collections.Counter(r['sequence_liabilities']['risk_tier']for r in pool)),
                 priority_risk_tiers=dict(collections.Counter(r['sequence_liabilities']['risk_tier']for r in selected[:20])),
                 reserve_risk_tiers=dict(collections.Counter(r['sequence_liabilities']['risk_tier']for r in selected[20:])),
                 priority_per_condition=dict(collections.Counter(r['task']for r in selected[:20])),
                 changed_priority_candidates=len({r['candidate_id']for r in selected[:20]}-old),
                 sequence_assembly_pass_count=len(selected),full_H2L2_3D_assembly_validated=False,
                 geometry_thresholds_relaxed=False,chemical_motif_preference_relaxed_in_priority=any(r['sequence_liabilities']['risk_tier']for r in selected[:20]),
                 experimental_binding_validated=False,original_strict_complete_pass_count=sum(r.get('original_full_gate_pass')is True for r in records))
    dump(BASE/'reports/summary.json',summary);print(json.dumps(summary,indent=2),flush=True)
    dump(BASE/'reports/policy.json',dict(
        inherited='Keep combined1000 v1 geometry, glycan clearance, sequence uniqueness and IL4RA novelty criteria.',
        additional_hard_checks=['complete IMGT VH and kappa VL','unchanged generic framework','CDRs match generated annotations','two conserved cysteines per variable chain','no full-variable N-X-S/T sequon','no Fv H/L heavy-atom clash below 2 A','intrachain Fv SG distances 1.5-3 A','no target clash below 2 A with rigid 5I19 Fab constants','no six-CDR>=0.8 AND H3>=0.8 similarity to any numbered therapeutic reference pair'],
        additional_conservative_heuristics=['no CDR hydrophobic run >=5 from AVILMFWY','absolute CDR K+R-D-E <=6'],
        soft_tiers='Tier 0: no CDR N[GS]/D[GS]. Tier 1: motif-bearing candidates admitted only after diverse tier-0 pool is exhausted; flags retained, no sequence edits.',
        ranking='Within tier: condition round-robin; CDR contact fraction descending, CDR contact count descending, chemical motif count ascending, framework RMSD ascending, glycan clearance descending, ID.',
        diversity='six-CDR normalized Levenshtein similarity <0.8 across all selected 40',
        oracle='No veto/reward due failed positive-control calibration; retain original pass labels and all evidence.',
        assembly='Exact official Human IgG1 + Human Ig Kappa constants; H2L2 connectivity sequence audit. No claim of full IgG 3D prediction.',
        novelty='Only local frozen therapeutic reference snapshot, not exhaustive novelty or patent clearance.'))
if __name__=='__main__':main()
