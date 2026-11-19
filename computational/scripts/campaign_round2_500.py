#!/usr/bin/env python3
"""Generate 500 additional designs, then automatically screen the combined 1000."""
import csv
import hashlib
import json
import os
import shlex
import signal
import subprocess
import sys
import time
import traceback
import uuid
from pathlib import Path

from campaign_local500 import ROOT, REPO, MODEL, dump, ledger, stop_process, analyze, novelty

BASE=ROOT/'campaign_local500_round2_20261118'
DB=BASE/'campaign.sqlite'
SEED_BASE=2026111900
STATE={}


def event(kind,**fields):
    stream=BASE/'events.jsonl'
    with stream.open('a') as handle:
        handle.write(json.dumps(dict(event_id=str(uuid.uuid4()),type=kind,**fields),allow_nan=False)+'\n')
        handle.flush();os.fsync(handle.fileno())
    with ledger.connect(DB) as conn:
        ledger.ingest_events(conn,DB,stream)


def command(task,i):
    return ['rfd4-generate','--config-path',str(BASE/'configs'),'--config-name','pilot_generation',
            f'inputs={ROOT}/inputs/annotated/{task}.cif',f'out_dir={BASE}/runs/batch_{i:02}',
            f'hf_ckpt_path={MODEL}','num_replicates=25','batch_size=1','num_workers=0',
            f'seed={SEED_BASE+i}','compress_outputs=false','ncpus_=2']


def initialize():
    import rfproteina
    from hydra import compose,initialize_config_dir
    from omegaconf import OmegaConf
    assert not DB.exists()
    (BASE/'configs/resolved').mkdir(parents=True,exist_ok=True)
    source=ROOT/'configs/pilot_generation.yaml'
    (BASE/'configs/pilot_generation.yaml').write_bytes(source.read_bytes())
    tasks=sorted(p.stem for p in (ROOT/'inputs/annotated').glob('*.cif'))
    assert len(tasks)==4
    hashes={str(ROOT/f'inputs/annotated/{t}.cif'):hashlib.sha256((ROOT/f'inputs/annotated/{t}.cif').read_bytes()).hexdigest() for t in tasks}
    plan=dict(target_new_designs=500,combined_target=1000,batches=20,designs_per_batch=25,tasks=tasks,
              seeds=list(range(SEED_BASE,SEED_BASE+20)),node='p-gc-3039',allocation_id='57427235',
              account='yjz5549_cr_default',partition='standard',qos='normal',cpus=2,gpu='a100:1',memory='48G',walltime='16:00:00',
              input_sha256=hashes,generation_hours_cap=12,
              method='Same four inputs and generation configuration as round 1; new random seeds.',
              evaluation_policy='Generate the requested additional 500, then uniformly screen all 1000. Uncalibrated OpenFold3 is retained as evidence, not an expansion gate or eligibility veto.',
              stopping_policy='Stop on infrastructure failure, incomplete batch, changed inputs, STOP marker or time cap. No automatic threshold changes.',
              combined_output=str(ROOT/'combined1000_20261118'))
    dump(BASE/'plan.json',plan)
    ledger.initialize(BASE,'IL4RA second 500 on current node',objective='Generate 500 additional paired Fvs and uniformly screen both rounds of 1000',
                      rfproteina_revision='aaccdd2ad0720c180468df5f9cf037a713d403e5',atomworks_revision='3.0.4')
    with ledger.connect(DB) as conn:
        for task in tasks:
            ledger.register_input(conn,DB,task,ROOT/f'inputs/annotated/{task}.cif',input_id=task,
                                  intent='Same validated canonical IL4RA and paired framework; six regenerated CDRs')
        conn.commit()
        for i in range(20):
            args=command(tasks[i%4],i)
            with initialize_config_dir(config_dir=str(BASE/'configs'),version_base=None):
                cfg=compose(config_name='pilot_generation',overrides=args[5:])
            path=BASE/f'configs/resolved/batch_{i:02}.yaml';path.write_text(OmegaConf.to_yaml(cfg,resolve=True))
            cid=ledger.register_config(conn,DB,f'batch_{i:02}',path,config_id=f'config_{i:02}',source_config=BASE/'configs/pilot_generation.yaml',overrides=args[5:])
            ledger.plan_run(conn,input_id=tasks[i%4],config_id=cid,phase='baseline',run_id=f'batch_{i:02}',command=shlex.join(args),
                            out_dir=str(BASE/f'runs/batch_{i:02}'),cluster='local',partition='standard',account='yjz5549_cr_default',qos='normal',
                            gres='gpu:a100:1',runtime='16:00:00',mem='48G',cpus_per_task='2',slurm_extra_options='Existing allocation 57427235 on p-gc-3039')
        conn.commit()
    print('Initialized 20 x 25 additional designs with non-overlapping seeds.',flush=True)


def save(**updates):
    STATE.update(updates,updated_at=time.strftime('%Y-%m-%dT%H:%M:%S%z'))
    dump(BASE/'status.json',STATE)


def run_generation(args,folder,timeout):
    execution=folder/'execution';execution.mkdir(parents=True,exist_ok=True)
    dump(execution/'command.json',dict(argv=args,timeout_seconds=timeout))
    with (execution/'process.log').open('w') as handle:
        proc=subprocess.Popen(args,stdout=handle,stderr=subprocess.STDOUT,start_new_session=True)
        end=time.monotonic()+timeout
        try:
            while proc.poll() is None:
                remaining=end-time.monotonic()
                if remaining<=0:
                    raise TimeoutError('Generation batch time budget exhausted.')
                try:
                    proc.wait(timeout=min(30,remaining))
                except subprocess.TimeoutExpired:
                    written=len(list((folder/'inference_samples').glob('*/*_generated.cif')))
                    save(current_batch_written=written,new_structures_written=STATE['generated']+written)
            if proc.returncode:
                raise RuntimeError(f'Generation exited {proc.returncode}; see {execution}/process.log')
        finally:
            stop_process(proc)


def main():
    assert os.environ.get('SLURM_JOB_ID')=='57427235'
    assert os.uname().nodename.split('.')[0]=='p-gc-3039'
    assert not (BASE/'status.json').exists(),'Reconcile existing outputs before restarting.'
    assert os.environ.get('CUDA_VISIBLE_DEVICES','').isdigit()
    plan=json.loads((BASE/'plan.json').read_text())
    dump(BASE/'allocation.json',{k:os.environ.get(k) for k in ['SLURM_JOB_ID','SLURM_STEP_ID','SLURM_JOB_NODELIST','CUDA_VISIBLE_DEVICES','SLURM_CPUS_PER_TASK']})
    STATE.update(status='running',job_id=os.environ['SLURM_JOB_ID'],step_id=os.environ.get('SLURM_STEP_ID'),node=os.uname().nodename,
                 started_at=time.strftime('%Y-%m-%dT%H:%M:%S%z'),generated=0,target_new_designs=500,combined_target=1000)
    save()
    deadline=time.monotonic()+plan['generation_hours_cap']*3600
    records=[]
    for i in range(20):
        if (ROOT/'STOP').exists() or (BASE/'STOP').exists():
            save(status='stopped',stop_reason='STOP marker');return
        if deadline-time.monotonic()<2400:
            save(status='stopped',stop_reason='Generation time budget');return
        task=plan['tasks'][i%4];run_id=f'batch_{i:02}';folder=BASE/'runs'/run_id
        source=ROOT/f'inputs/annotated/{task}.cif'
        assert hashlib.sha256(source.read_bytes()).hexdigest()==plan['input_sha256'][str(source)]
        save(stage=run_id,current_batch_written=0,new_structures_written=STATE['generated'])
        event('run_status',run_id=run_id,status='running')
        failure=None
        try:
            run_generation(command(task,i),folder,min(5400,deadline-time.monotonic()))
        except Exception as exc:
            failure=exc
        files=list((folder/'inference_metrics').glob('per_example_metrics*.csv'))
        rows=list(csv.DictReader(files[0].open())) if len(files)==1 else []
        for j,row in enumerate(rows):
            assert row['task']==task
            path=Path(row['designed_structure_path']);did=f'b{i:02}_d{j:03}'
            event('design',design_id=did,run_id=run_id,cif_path=str(path),chain_of_interest='H,L',seed=SEED_BASE+i,example_id=row['example_id'])
            try:
                r=analyze(path,task);novelty(r)
            except Exception as exc:
                r=dict(path=str(path),task=task,triage_pass=False,error=repr(exc))
            r.update(design_id=did,seed=SEED_BASE+i,example_id=row['example_id'])
            records.append(r)
            values={k:int(v) if isinstance(v,bool) else v for k,v in r.items() if isinstance(v,(int,float,bool))}
            event('metrics',design_id=did,evaluator='geometry_sequence_v1',values=values)
        STATE['generated']+=len(rows)
        event('run_status',run_id=run_id,status='failed' if failure or len(rows)!=25 else 'completed')
        dump(BASE/'reports/all_designs.json',records)
        save(triage_passers=sum(r['triage_pass'] for r in records),analysis_errors=sum('error' in r for r in records),new_structures_written=STATE['generated'])
        if failure:
            raise failure
        assert len(rows)==25,'Incomplete batch; preserve output and stop.'
    assert STATE['generated']==500
    save(stage='combined1000_screening')
    with (BASE/'combined_screening.log').open('w') as handle:
        subprocess.run([sys.executable,'-B',str(ROOT/'scripts/screen_combined1000.py')],stdout=handle,stderr=subprocess.STDOUT,check=True,timeout=7200)
    summary=json.loads((ROOT/'combined1000_20261118/summary.json').read_text())
    save(status='completed',stage='combined1000_screening_complete',combined_summary=summary)


if __name__=='__main__':
    if '--initialize' in sys.argv:
        initialize()
    else:
        def interrupted(signum,frame):
            raise KeyboardInterrupt(f'Signal {signum}')
        signal.signal(signal.SIGTERM,interrupted)
        try:
            main()
        except BaseException as exc:
            if STATE:
                save(status='failed',error=repr(exc),traceback=traceback.format_exc())
            raise
