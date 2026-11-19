#!/usr/bin/env python3
"""Uniform, transparent review ranking of two 500-design rounds."""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

from campaign_local500 import ROOT, dump, novelty, similarity, cdr_sequence

ROUND1 = ROOT/'campaign_local500_20261117'
ROUND2 = ROOT/'campaign_local500_round2_20261118'
BASE = ROOT/'combined1000_20261118'
POLICY = {
    'id': 'combined1000_geometry_sequence_review_v1',
    'eligibility': 'Original geometry/sequence triage AND zero CDR glycosylation motifs AND no frozen-reference similarity flag',
    'oracle_policy': 'Evidence only: positive-control calibration failed; no oracle veto, reward, or assumption that untested means failed',
    'rank_within_condition': ['cdr_contact_fraction descending', 'cdr_contact_residues descending',
                              'paired_framework_rmsd ascending', 'glycan_site_min_distance descending', 'candidate_id ascending'],
    'diversity': 'Greedy six-CDR normalized edit similarity <0.8 between all selected entries',
    'balance': 'Round-robin over four conditions; use another condition when one has no eligible diverse entries',
    'queue_limit': 20,
    'reserve_limit': 20,
    'interpretation': 'Provisional computational review priority; not binding probability, experimental affinity, or challenge-ready IgG',
}


def ranking_key(r):
    return (-r['cdr_contact_fraction'], -r['cdr_contact_residues'],
            r['paired_framework_rmsd'], -r['glycan_site_min_distance'], r['candidate_id'])


def triage_reasons(r):
    if 'error' in r:
        return ['analysis_error']
    checks = [
        ('steric_clash', r['clash_pairs_below_2A']==0),
        ('hotspot_coverage', r['hotspot_contacts']>=2),
        ('target_contacts', r['target_contact_residues']>=5),
        ('cdr_contacts', r['cdr_contact_residues']>=3),
        ('cdr_contact_fraction', r['cdr_contact_fraction']>=0.5),
        ('target_deformation', r['raw_target_rmsd']<=2.5),
        ('framework_deformation', r['paired_framework_rmsd']<=2.5),
        ('glycan_site_clearance', r['glycan_site_min_distance']>=4),
        ('extra_cdr_cysteine', r['extra_cdr_cysteines']==0),
        ('all_donor_cdrs_retained', r['donor_identical_cdr_count']<6),
    ]
    reasons=[name for name,passed in checks if not passed]
    assert bool(r['triage_pass']) == (not reasons), 'Cached triage disagrees with common screening policy.'
    if r['cdr_glycosylation_motifs']:
        reasons.append('cdr_glycosylation_motif')
    if r['known_reference_similar']:
        reasons.append('frozen_reference_similarity')
    return reasons


def diverse_balanced(pool, limit=40):
    groups={task: sorted([r for r in pool if r['task']==task],key=ranking_key)
            for task in sorted({r['task'] for r in pool})}
    selected=[]
    while len(selected)<limit:
        progressed=False
        for task in groups:
            while groups[task]:
                r=groups[task].pop(0)
                if all(similarity(cdr_sequence(r),cdr_sequence(other))<0.8 for other in selected):
                    selected.append(r)
                    progressed=True
                    break
            if len(selected)>=limit:
                break
        if not progressed:
            break
    return selected


def build_records(rounds, diagnostic_results=()):
    records=[]
    for label,source_rows in rounds:
        for original in source_rows:
            if original['design_id'].startswith('pilot'):
                continue
            r=dict(original)
            r.update(round=label,candidate_id=label+'_'+original['design_id'],source_design_id=original['design_id'])
            # Reapply the same reference screen to both rounds, preserving the
            # old complete-pass label instead of silently relabeling it.
            r['original_full_gate_pass']=original.get('computational_pass')
            r['oracle_evidence']=[]
            if 'oracle' in original:
                r['oracle_evidence'].append({'seed':20261118,'metrics':original['oracle'],'source':'round1_original'})
            if label=='R1':
                for result in diagnostic_results:
                    if result['design_id']==original['design_id']:
                        r['oracle_evidence'].append({'seed':result['seed'],'metrics':result.get('oracle',result['metrics']),
                                                     'source':'continuation_recheck'})
            r['oracle_status']='evaluated_uncalibrated' if r['oracle_evidence'] else 'not_evaluated'
            if 'error' not in r:
                novelty(r)
            r['screen_reasons']=triage_reasons(r)
            r['review_eligible']=not r['screen_reasons']
            r['duplicate_of']=None
            records.append(r)
    assert len({r['candidate_id'] for r in records})==len(records)
    # Keep the best eligible geometry for any identical VH/VL pair; both rows
    # remain in the 1000-design report and the duplicate is explicitly linked.
    seen={}
    ordered=sorted(records,key=lambda r:(not r['review_eligible'],ranking_key(r) if 'error' not in r else (0,0,0,0,r['candidate_id'])))
    for r in ordered:
        key=r.get('sequence_hash')
        if key is None:
            continue
        if key in seen:
            r['duplicate_of']=seen[key]
        else:
            seen[key]=r['candidate_id']
    return records


def main():
    a=json.loads((ROUND1/'reports/all_designs.json').read_text())
    b=json.loads((ROUND2/'reports/all_designs.json').read_text())
    assert len([r for r in a if not r['design_id'].startswith('pilot')])==500
    assert len(b)==500, 'Unified selection runs only after all 500 new designs are accounted for.'
    diagnostic=json.loads((ROOT/'continuation_20261118/reports/diagnostic_refolds.json').read_text())
    records=build_records([('R1',a),('R2',b)],diagnostic)
    assert len(records)==1000
    for r in records:
        if 'sha256' in r:
            assert hashlib.sha256(Path(r['path']).read_bytes()).hexdigest()==r['sha256'], r['candidate_id']
    pool=[r for r in records if r['review_eligible'] and not r['duplicate_of']]
    selected=diverse_balanced(pool)
    queues={r['candidate_id']:('priority_review' if i<20 else 'reserve',i+1) for i,r in enumerate(selected)}
    for r in records:
        r['queue'],r['selection_order']=queues.get(r['candidate_id'],('unselected',None))
    dump(BASE/'policy.json',POLICY)
    dump(BASE/'all_designs.json',records)
    dump(BASE/'priority_review_20.json',selected[:20])
    dump(BASE/'reserves_20.json',selected[20:])
    fields=['candidate_id','round','source_design_id','task','VH','VL','triage_pass','review_eligible',
            'screen_reasons','duplicate_of','queue','selection_order','original_full_gate_pass','oracle_status',
            'clash_pairs_below_2A','hotspot_contacts','target_contact_residues','cdr_contact_residues',
            'cdr_contact_fraction','raw_target_rmsd','paired_framework_rmsd','glycan_site_min_distance',
            'extra_cdr_cysteines','cdr_glycosylation_motifs','known_reference_similar','path','error']
    for name,rows in [('all1000_ranked.csv',sorted(records,key=lambda r:(not r['review_eligible'],r['duplicate_of'] is not None,
                         ranking_key(r) if 'error' not in r else (0,0,0,0,r['candidate_id'])))),
                      ('priority_review_20.csv',selected[:20]),('reserves_20.csv',selected[20:])]:
        with (BASE/name).open('w') as handle:
            writer=csv.DictWriter(handle,fieldnames=fields,extrasaction='ignore');writer.writeheader()
            for r in rows:
                writer.writerow(dict(r,screen_reasons=';'.join(r['screen_reasons'])))
    summary={'status':'completed','total_generated':len(records),
             'unique_valid_sequence_pairs':len({r['sequence_hash'] for r in records if r.get('sequence_hash')}),
             'analysis_errors':sum('error' in r for r in records),
             'triage_passers':sum(bool(r.get('triage_pass')) for r in records),
             'unique_review_eligible':len(pool),'priority_review_count':len(selected[:20]),'reserve_count':len(selected[20:]),
             'per_round':{label:dict(generated=sum(r['round']==label for r in records),
                                    triage_passers=sum(r['round']==label and bool(r.get('triage_pass')) for r in records))
                          for label in ['R1','R2']},
             'selected_per_condition':dict(Counter(r['task'] for r in selected[:20])),
             'policy':POLICY,'interpretation':'Review shortlist only; broader novelty/developability and full IgG assembly still required.'}
    dump(BASE/'summary.json',summary)
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    main()
