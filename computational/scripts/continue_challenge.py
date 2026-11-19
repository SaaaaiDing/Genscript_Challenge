#!/usr/bin/env python3
"""Bounded diagnostic continuation after the 500-design scientific stop."""
import hashlib
import csv
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time
import traceback
from pathlib import Path

import numpy as np
from biotite.structure.io.pdbx import CIFFile, get_structure
from campaign_local500 import REPO, ROOT, MODEL, dump, stop_process, run_process, refold_metrics, novelty
from quality_screen import analyze, fit

BASE = ROOT / 'continuation_20261118'
STATE = {'status': 'running', 'stage': 'diagnostics', 'generated': 0, 'refolded': 0,
         'job_id': os.environ.get('SLURM_JOB_ID'), 'node': os.uname().nodename,
         'started_at': time.strftime('%Y-%m-%dT%H:%M:%S%z')}


def save(**updates):
    STATE.update(updates, updated_at=time.strftime('%Y-%m-%dT%H:%M:%S%z'))
    dump(BASE / 'status.json', STATE)


def basic_metrics(source, predicted, confidence):
    a = get_structure(CIFFile.read(source), model=1, use_author_fields=True)
    b = get_structure(CIFFile.read(predicted), model=1, use_author_fields=True, extra_fields=['b_factor'])
    keys = lambda x: list(zip(x.chain_id, x.res_id, x.atom_name))
    mapping = {k:i for i,k in enumerate(keys(b))}
    b = b[[mapping[k] for k in keys(a)]]
    mask = (a.chain_id == 'T') & np.isin(a.atom_name, ['N','CA','C','O'])
    rot, translation = fit(b.coord[mask], a.coord[mask])
    aligned = b.coord @ rot + translation
    antibody = (a.chain_id != 'T') & np.isin(a.atom_name, ['N','CA','C','O'])
    metrics = {'target_aligned_fv_rmsd': float(np.sqrt(np.mean(np.sum((aligned[antibody] - a.coord[antibody])**2,axis=1))))}
    for chain in ['H','L']:
        metrics[chain + '_plddt'] = float(b.b_factor[(b.chain_id == chain)&(b.atom_name == 'CA')].mean())
    for key,value in json.loads(Path(confidence).read_text()).items():
        match = re.fullmatch(r'(min_ipae|iptm)_([THL])(?:_\d+)?_([THL])(?:_\d+)?',key)
        if match:
            metric,c,d = match.groups()
            metrics[f'{metric}_{"".join(sorted([c,d]))}'] = float(value)
    assert all(np.isfinite(v) for v in metrics.values())
    return metrics


def evaluate(chosen, tag, seed):
    folder = BASE / 'oracles' / tag
    inputs = folder / 'inputs'
    inputs.mkdir(parents=True, exist_ok=False)
    for r in chosen:
        shutil.copy2(r['path'], inputs / (r['design_id'] + '.cif'))
    env = os.environ.copy()
    socket_dir = Path('/tmp') / ('il4ra_continue_' + os.environ['SLURM_JOB_ID'])
    socket_dir.mkdir(exist_ok=True)
    sock = socket_dir / 'o.sock'
    assert not sock.exists(), 'Unexpected live or stale continuation socket; inspect before reuse.'
    env.update(TMPDIR=str(socket_dir), AF3_LAUNCHER_SOCKET_TEMPLATE=str(sock),
               APPTAINER_BIND_PATHS='/storage/group/yjz5549/default,/storage/home/sqd5856/Desktop/default')
    env['PYTHONPATH'] = str(ROOT/'scripts') + os.pathsep + env.get('PYTHONPATH', '')
    env['IL4RA_RAW_ORACLE_ARCHIVE'] = str(folder/'raw_outputs')
    args = [sys.executable, str(REPO / 'src/rfproteina/af3/af3_daemon.py'), 'serve',
            '--socket-path',str(sock), '--gpu-device',os.environ['CUDA_VISIBLE_DEVICES'], '--lazy-warmup']
    save(stage=tag)
    with (folder/'daemon.log').open('w') as log:
        daemon = subprocess.Popen(args, stdout=log, stderr=subprocess.STDOUT, env=env, start_new_session=True)
        try:
            for _ in range(600):
                if sock.exists():
                    break
                assert daemon.poll() is None
                time.sleep(.1)
            assert sock.exists()
            run_process(['rfd4-evaluate','--config-path',str(ROOT/'configs'),'--config-name','paired_evaluation',
                         f'inputs={inputs}',f'out_dir={folder}/results','compress_outputs=false',
                         'metrics.refolding_oracles.oracle_ppi._target_=il4ra_oracle.PairedAF3Oracle',
                         f'+metrics.refolding_oracles.oracle_ppi.seed={seed}'], folder/'execution',7200,env)
        finally:
            stop_process(daemon)
            if sock.exists():
                sock.unlink()
    results=[]
    for r in chosen:
        predicted=list((folder/'results').rglob(r['design_id']+'*_of3_ppi_refolded.cif'))
        conf=list((folder/'results').rglob(r['design_id']+'*_of3_ppi_confidences.json'))
        assert len(predicted)==len(conf)==1, f'Missing oracle evidence: {r["design_id"]}'
        result={'design_id':r['design_id'],'seed':seed,'predicted':str(predicted[0]),'confidence':str(conf[0]),
                'role':r.get('role','de_novo_candidate_recheck'),
                'metrics':basic_metrics(r['path'],predicted[0],conf[0])}
        if r.get('task'):
            metrics,geometry=refold_metrics(r,predicted[0],conf[0])
            result.update(oracle=metrics,refold_geometry=geometry)
        results.append(result)
    STATE['refolded']+=len(results)
    dump(folder/'summary.json',results)
    save()
    return results


def meets_original_gates(record, result):
    m=result['oracle']
    return bool(record['triage_pass'] and not record['known_reference_similar'] and
                record['cdr_glycosylation_motifs']==0 and m['refold_geometry_pass'] and
                m['target_aligned_fv_rmsd']<=2.5 and m['vh_plddt']>=80 and
                m['vl_plddt']>=80 and m['target_interface_pass'] and m['pairing_iptm']>=60)


def refinement_pilot(parents):
    """Test 32 partial-diffusion variants of our own four de novo designs."""
    from campaign_local500 import similarity, cdr_sequence
    originals=json.loads((ROOT/'campaign_local500_20261117/reports/all_designs.json').read_text())
    original_hashes={r['sequence_hash'] for r in originals if 'sequence_hash' in r}
    records=[]
    seen=set(original_hashes)
    for i,parent in enumerate(parents):
        if (BASE/'STOP').exists() or (ROOT/'STOP').exists():
            save(status='stopped',stop_reason='STOP marker')
            return
        source=Path(parent['path'])
        assert source.is_relative_to(ROOT/'campaign_local500_20261117/runs')
        assert hashlib.sha256(source.read_bytes()).hexdigest()==parent['sha256']
        folder=BASE/'refinement'/f'parent_{i:02}'
        args=['rfd4-generate','--config-path',str(ROOT/'configs'),'--config-name','pilot_generation',
              f'inputs={source}',f'out_dir={folder}',f'hf_ckpt_path={MODEL}',
              'num_replicates=8','batch_size=1','num_workers=0',f'seed={2026111840+i}',
              'compress_outputs=false','ncpus_=2','generation.args.partial_trajectory_fraction_steps=0.5']
        save(stage=f'partial_diffusion_parent_{i:02}')
        run_process(args,folder/'execution',3600)
        files=list((folder/'inference_metrics').glob('per_example_metrics*.csv'))
        assert len(files)==1
        rows=list(csv.DictReader(files[0].open()))
        assert len(rows)==8
        for j,row in enumerate(rows):
            r=analyze(row['designed_structure_path'],parent['task'])
            novelty(r)
            r.update(design_id=f'refine_p{i:02}_d{j:02}',parent_design_id=parent['design_id'],
                     parent_path=str(source),seed=2026111840+i,example_id=row['example_id'],
                     duplicate_pair=r['sequence_hash'] in seen)
            seen.add(r['sequence_hash'])
            records.append(r)
        STATE['generated']+=len(rows)
        dump(BASE/'reports/refinement_designs.json',records)
        save(triage_passers=sum(r['triage_pass'] for r in records))
    # Screen up to four distinct variants per parent, retaining the original
    # novelty and sequence gates and recording every generated failure.
    chosen=[]
    for parent in parents:
        pool=[r for r in records if r['parent_design_id']==parent['design_id'] and r['triage_pass'] and
              not r['known_reference_similar'] and not r['cdr_glycosylation_motifs'] and not r['duplicate_pair']]
        pool.sort(key=lambda r:(-r['cdr_contact_fraction'],r['paired_framework_rmsd']))
        chosen.extend(pool[:4])
    if chosen:
        results=evaluate(chosen,'refinement_seed20261120',20261120)
        by_id={r['design_id']:r for r in records}
        for result in results:
            record=by_id[result['design_id']]
            record['oracle']=result['oracle']
            record['refold_geometry']=result['refold_geometry']
            record['computational_pass']=meets_original_gates(record,result)
        dump(BASE/'reports/refinement_designs.json',records)
    eligible=sorted([r for r in records if r.get('computational_pass')],key=lambda r:r['oracle']['target_aligned_fv_rmsd'])
    diverse=[]
    for r in eligible:
        if all(similarity(cdr_sequence(r),cdr_sequence(other))<.8 for other in diverse):
            diverse.append(r)
        if len(diverse)==20:
            break
    dump(BASE/'exports/computational_review_queue.json',diverse)
    save(status='completed',stage='refinement_pilot_complete',computational_passers=len(eligible),
         diverse_review_candidates=len(diverse),
         note='32-variant hypothesis test complete. Review results before further scaling; challenge submission is not performed.')


def main():
    assert os.environ.get('SLURM_JOB_ID')=='57427235'
    assert os.uname().nodename.split('.')[0]=='p-gc-3039'
    assert not (BASE/'status.json').exists(), 'Reconcile an existing continuation before restart.'
    save()
    control={'design_id':'positive_control_6wgl','path':str(BASE/'benchmark/inputs/positive_control_6wgl.cif'),
             'role':'evaluation_only_positive_control_not_a_submission_candidate'}
    candidates=json.loads((BASE/'reports/recheck_candidates.json').read_text())
    results=evaluate([control],'control_seed20261118',20261118)
    results+=evaluate([control]+candidates,'recheck_seed20261119',20261119)
    dump(BASE/'reports/diagnostic_refolds.json',results)
    controls=[r['metrics'] for r in results if r['design_id']==control['design_id']]
    # Establish operational sensitivity on the known interaction before testing
    # a new generation hypothesis. These are the original interface thresholds.
    calibrated=all(m['target_aligned_fv_rmsd']<=2.5 and m['iptm_HL']>=60 and
                   min(m['H_plddt'],m['L_plddt'])>=80 and
                   any(m[f'iptm_{c}T']>=60 and m[f'min_ipae_{c}T']<=2 for c in ['H','L'])
                   for m in controls)
    save(positive_control_recovered_both_seeds=calibrated)
    if not calibrated:
        save(status='diagnostics_completed_needs_review',stage='oracle_calibration_review',
             note='The known positive control did not satisfy the original interface/pose thresholds in both seeds. Reassess evaluator before more generation.')
        return
    refinement_pilot(candidates)


if __name__=='__main__':
    def interrupted(signum, frame):
        raise KeyboardInterrupt(f'Signal {signum}')
    signal.signal(signal.SIGTERM,interrupted)
    try:
        main()
    except BaseException as exc:
        save(status='failed',error=repr(exc),traceback=traceback.format_exc())
        raise
