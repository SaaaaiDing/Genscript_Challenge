#!/usr/bin/env python3
"""Bounded, serial paired-Fv campaign with immutable events and conservative gates."""
import csv
import hashlib
import json
import os
import re
import shlex
import shutil
import signal
import subprocess
import sys
import time
import traceback
import uuid
from pathlib import Path

import numpy as np
from biotite.structure.io.pdbx import CIFFile, get_structure
from quality_screen import ROOT, analyze, fit

REPO = Path('/storage/group/yjz5549/default/RFD4/RFD4-Proteina')
BASE = ROOT / 'campaign_local500_20261117'
BATCHES = 20
DESIGNS_PER_BATCH = 25
SEED_BASE = 2026111700
DB = BASE / 'campaign.sqlite'
MODEL = '/storage/group/yjz5549/default/RFD4/huggingface/hub/models--nvidia--RFD4-Proteina/snapshots/a06e2509c7033b010fc43614e9fd840303588502'
sys.path.insert(0, str(REPO / 'skills/rfd4-proteina-campaign/scripts'))
import campaign_ledger as ledger


def dump(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)


def event(kind, **fields):
    stream = BASE / 'events.jsonl'
    payload = dict(event_id=str(uuid.uuid4()), type=kind, **fields)
    with stream.open('a') as f:
        f.write(json.dumps(payload, allow_nan=False) + '\n')
        f.flush()
        os.fsync(f.fileno())
    with ledger.connect(DB) as conn:
        ledger.ingest_events(conn, DB, stream)


def command(task, index):
    return ['rfd4-generate', '--config-path', str(ROOT / 'configs'), '--config-name', 'pilot_generation',
            f'inputs={ROOT}/inputs/annotated/{task}.cif', f'out_dir={BASE}/runs/batch_{index:02}',
            f'hf_ckpt_path={MODEL}', f'num_replicates={DESIGNS_PER_BATCH}', 'batch_size=1', 'num_workers=0',
            f'seed={SEED_BASE + index}', 'compress_outputs=false', 'ncpus_=2']


def initialize():
    import rfproteina
    from hydra import compose, initialize_config_dir
    from omegaconf import OmegaConf
    assert not DB.exists(), 'Campaign already exists; do not initialize twice.'
    BASE.mkdir(parents=True, exist_ok=True)
    tasks = sorted(p.stem for p in (ROOT / 'inputs/annotated').glob('*.cif'))
    assert len(tasks) == 4
    plan = dict(target_new_designs=500, batches=BATCHES, designs_per_batch=DESIGNS_PER_BATCH, tasks=tasks,
                generation_hours_cap=12, total_hours_cap=35.25, oracle_count_cap=80,
                zero_passer_stop_after_oracles=24, zero_triage_stop_after_designs=160,
                pilot_seconds_per_design=43.875, expected_generation_hours=7.5,
                account='yjz5549_cr_default', partition='standard', qos='normal',
                node='p-gc-3039', allocation_id='57427235', gpu='a100:1', cpus=2, memory='48G', walltime='35:15:00',
                execution='Existing allocation on the current node, explicitly requested by user',
                novelty_scope='Frozen TheraSAbDab IL4R reference set; not exhaustive patent clearance',
                stopping_policy='Stop on missing oracle evidence or infrastructure failure; never loosen gates automatically',
                selection='Balanced four conditions; sequence deduplication and greedy six-CDR diversity; at most 20 review candidates')
    dump(BASE / 'plan.json', plan)
    criteria = {'id': 'paired_fv_v1', 'name': 'Conservative computational review eligibility', 'criteria': [
        {'metric': key, 'operator': op, 'threshold': val, 'missing': 'fail'}
        for key, op, val in [('triage_pass', '==', 1), ('refold_geometry_pass', '==', 1),
                            ('target_aligned_fv_rmsd', '<=', 2.5), ('vh_plddt', '>=', 80),
                            ('vl_plddt', '>=', 80), ('target_interface_pass', '==', 1),
                            ('pairing_iptm', '>=', 60), ('cdr_glycosylation_motifs', '==', 0),
                            ('known_reference_similar', '==', 0)]]}
    dump(BASE / 'criteria_v1.json', criteria)
    ledger.initialize(BASE, 'IL4RA paired-Fv local 500 campaign',
                      objective='Up to 500 de novo six-CDR paired Fv designs and a diverse computational review queue',
                      rfproteina_revision='aaccdd2ad0720c180468df5f9cf037a713d403e5', atomworks_revision='3.0.4')
    with ledger.connect(DB) as conn:
        for task in tasks:
            ledger.register_input(conn, DB, task, ROOT / f'inputs/annotated/{task}.cif', input_id=task,
                                  intent='Canonical IL4RA ECD; two-chain generic framework; all six CDRs regenerated')
        conn.commit()
        pilot_config = ledger.register_config(conn, DB, 'Completed pilot configuration',
            ROOT / 'runs/initial_designs/hydra/.hydra/config.yaml', config_id='pilot_config',
            source_config=ROOT / 'configs/pilot_generation.yaml')
        conn.commit()
        ledger.register_criteria(conn, criteria)
        for i in range(BATCHES):
            task = tasks[i % 4]
            args = command(task, i)
            with initialize_config_dir(config_dir=str(ROOT / 'configs'), version_base=None):
                resolved = compose(config_name='pilot_generation', overrides=args[5:])
            config = BASE / f'configs/batch_{i:02}.yaml'
            config.parent.mkdir(exist_ok=True)
            config.write_text(OmegaConf.to_yaml(resolved, resolve=True))
            config_id = ledger.register_config(conn, DB, f'batch_{i:02}', config, config_id=f'config_{i:02}',
                                               source_config=ROOT / 'configs/pilot_generation.yaml', overrides=args[5:])
            ledger.plan_run(conn, input_id=task, config_id=config_id, phase='baseline', run_id=f'batch_{i:02}',
                            command=shlex.join(args), out_dir=str(BASE / f'runs/batch_{i:02}'),
                            cluster='local', partition='standard', account='yjz5549_cr_default', qos='normal',
                            gres='gpu:a100:1', runtime='35:15:00', mem='48G', cpus_per_task='2',
                            slurm_extra_options='Existing allocation 57427235 on p-gc-3039; no sbatch submission')
        conn.commit()
    print('Campaign initialized: 20 batches of 25, 125 designs per condition.', flush=True)


def similarity(a, b):
    """One minus normalized Levenshtein distance; no framework contribution."""
    prev = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        cur = [i]
        for j, y in enumerate(b, 1):
            cur.append(min(cur[-1] + 1, prev[j] + 1, prev[j-1] + (x != y)))
        prev = cur
    return 1 - prev[-1] / max(len(a), len(b), 1)


def cdr_sequence(record):
    return '|'.join(record['cdrs'][k] for k in ('H1', 'H2', 'H3', 'L1', 'L2', 'L3'))


def novelty(record):
    refs = json.loads((ROOT / 'inputs/raw/known_il4ra_cdrs.json').read_text())
    assert refs, 'Empty novelty reference set.'
    scores = [(similarity(cdr_sequence(record), cdr_sequence(ref)),
               similarity(record['cdrs']['H3'], ref['cdrs']['H3']), ref['name']) for ref in refs]
    record['known_reference_similar'] = int(any(a >= .8 and b >= .8 for a, b, _ in scores))
    record['closest_known_reference'] = max(scores)[2]
    record['max_known_six_cdr_similarity'] = max(s[0] for s in scores)
    record['novelty_status'] = 'frozen_reference_screen_only; broader_novelty_review_required'


def stop_process(proc):
    if proc.poll() is None:
        os.killpg(proc.pid, signal.SIGTERM)
        try:
            proc.wait(timeout=20)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait()


def run_process(args, folder, timeout, env=None):
    folder.mkdir(parents=True, exist_ok=True)
    dump(folder / 'command.json', dict(argv=args, timeout_seconds=timeout))
    with (folder / 'process.log').open('w') as log:
        proc = subprocess.Popen(args, stdout=log, stderr=subprocess.STDOUT, env=env, start_new_session=True)
        try:
            code = proc.wait(timeout=max(1, timeout))
            if code:
                raise RuntimeError(f'Command failed with exit {code}: {folder}/process.log')
        finally:
            stop_process(proc)


def refold_metrics(design, predicted, confidence):
    original = get_structure(CIFFile.read(design['path']), model=1, use_author_fields=True)
    folded = get_structure(CIFFile.read(predicted), model=1, use_author_fields=True, extra_fields=['b_factor'])
    assert set(folded.chain_id) == {'T', 'H', 'L'}
    keys = lambda a: list(zip(a.chain_id, a.res_id, a.atom_name))
    mapping = {k: i for i, k in enumerate(keys(folded))}
    folded = folded[[mapping[k] for k in keys(original)]]
    bb = np.isin(original.atom_name, ['N', 'CA', 'C', 'O'])
    target = (original.chain_id == 'T') & bb
    r, t = fit(folded.coord[target], original.coord[target])
    coord = folded.coord @ r + t
    m = {}
    for label, chains in [('fv', ['H', 'L']), ('vh', ['H']), ('vl', ['L'])]:
        mask = np.isin(original.chain_id, chains) & bb
        m[f'target_aligned_{label}_rmsd'] = float(np.sqrt(np.mean(np.sum((coord[mask] - original.coord[mask])**2, axis=1))))
    for label, chain in [('vh', 'H'), ('vl', 'L')]:
        m[label + '_plddt'] = float(folded.b_factor[(folded.chain_id == chain) & (folded.atom_name == 'CA')].mean())
    conf = json.loads(Path(confidence).read_text())
    pairs = {}
    for key, value in conf.items():
        match = re.fullmatch(r'(min_ipae|iptm)_([THL])(?:_\d+)?_([THL])(?:_\d+)?', key)
        if match:
            metric, a, b = match.groups()
            pairs[(metric, ''.join(sorted([a, b])))] = float(value)
    for pair in ['HT', 'LT', 'HL']:
        for metric in ['min_ipae', 'iptm']:
            value = pairs[(metric, pair)]
            assert np.isfinite(value), f'Missing confidence for {metric} {pair}'
            m[f'{metric}_{pair}'] = value
    # A target-contacting chain must satisfy BOTH interface confidence thresholds.
    from scipy.spatial.distance import cdist
    target_heavy = (folded.chain_id == 'T') & ~np.isin(folded.element, ['H', 'D'])
    contact = {c: bool((cdist(folded.coord[target_heavy], folded.coord[(folded.chain_id == c) & ~np.isin(folded.element, ['H', 'D'])]) <= 5).any()) for c in ['H', 'L']}
    m['target_interface_pass'] = int(any(contact[c] and m[f'min_ipae_{c}T'] <= 2 and m[f'iptm_{c}T'] >= 60 for c in ['H', 'L']))
    m['pairing_iptm'] = m['iptm_HL']
    geometry = analyze(predicted, design['task'])
    m['refold_geometry_pass'] = int(geometry['geometry_pass'])
    assert all(np.isfinite(v) for v in m.values())
    return m, geometry


def main():
    assert os.environ.get('SLURM_JOB_ID'), 'GPU work requires a scheduler allocation.'
    assert os.uname().nodename.split('.')[0] == 'p-gc-3039', 'User requested the current node.'
    assert os.environ['SLURM_JOB_ID'] == '57427235', 'Use the existing allocation only.'
    assert not (BASE / 'status.json').exists(), 'Existing run requires explicit reconciliation, not blind restart.'
    visible = os.environ.get('CUDA_VISIBLE_DEVICES', '').split(',')
    assert len(visible) == 1 and visible[0].isdigit(), 'Require exactly one scheduler-masked numeric GPU; do not guess a device.'
    gpu_device = visible[0]
    dump(BASE / 'allocation.json', {key: os.environ.get(key) for key in
         ['SLURM_JOB_ID', 'SLURM_JOB_NODELIST', 'SLURM_JOB_GPUS', 'CUDA_VISIBLE_DEVICES', 'SLURM_CPUS_PER_TASK']})
    plan = json.loads((BASE / 'plan.json').read_text())
    start = time.monotonic()
    deadline = start + plan['total_hours_cap'] * 3600
    generation_deadline = start + plan['generation_hours_cap'] * 3600
    records, tested, selected_history = [], set(), []
    state = {'status': 'running', 'job_id': os.environ['SLURM_JOB_ID'], 'node': os.uname().nodename,
             'target_new_designs': 500, 'stage': 'pilot_oracle_validation',
             'started_at': time.strftime('%Y-%m-%dT%H:%M:%S%z'),
             'generated': 0, 'refolded': 0, 'computational_passers': 0}
    dump(BASE / 'status.json', state)
    def available():
        return deadline - time.monotonic()
    def save():
        state.update(triage_passers=sum(r.get('triage_pass', False) for r in records),
                     computational_passers=sum(r.get('computational_pass', False) for r in records),
                     elapsed_hours=(time.monotonic() - start) / 3600,
                     updated_at=time.strftime('%Y-%m-%dT%H:%M:%S%z'))
        dump(BASE / 'status.json', state)
        dump(BASE / 'reports/all_designs.json', records)
        eligible = sorted([r for r in records if r.get('computational_pass')], key=lambda r: r['oracle']['target_aligned_fv_rmsd'])
        diverse = []
        for r in eligible:
            if all(similarity(cdr_sequence(r), cdr_sequence(x)) < .8 for x in diverse):
                diverse.append(r)
            if len(diverse) == 20:
                break
        dump(BASE / 'exports/computational_review_queue.json', diverse)
    def add(path, task, run_id, design_id, seed=None, example_id=None):
        event('design', design_id=design_id, run_id=run_id, cif_path=str(path), chain_of_interest='H,L', seed=seed, example_id=example_id)
        try:
            record = analyze(path, task)
            novelty(record)
        except Exception as exc:
            record = dict(path=str(path), task=task, triage_pass=False, error=repr(exc))
        record['design_id'] = design_id
        records.append(record)
        values = {k: int(v) if isinstance(v, bool) else v for k, v in record.items() if isinstance(v, (int, float, bool))}
        event('metrics', design_id=design_id, evaluator='geometry_sequence_v1', values=values)
    def evaluate(chosen):
        if not chosen:
            return
        idx = state['refolded']
        folder = BASE / f'oracles/set_{idx:03}'
        inputs = folder / 'inputs'
        inputs.mkdir(parents=True)
        for r in chosen:
            shutil.copy2(r['path'], inputs / (r['design_id'] + '.cif'))
        env = os.environ.copy()
        socket_dir = Path('/tmp') / ('il4ra500_' + os.environ['SLURM_JOB_ID'])
        socket_dir.mkdir(exist_ok=True)
        sock = socket_dir / 'o.sock'
        if sock.exists():
            sock.unlink()
        env['TMPDIR'] = str(socket_dir)
        env['AF3_LAUNCHER_SOCKET_TEMPLATE'] = str(sock)
        env['APPTAINER_BIND_PATHS'] = '/storage/group/yjz5549/default,/storage/home/sqd5856/Desktop/default'
        daemon_args = [sys.executable, str(REPO / 'src/rfproteina/af3/af3_daemon.py'), 'serve', '--socket-path', str(sock), '--gpu-device', gpu_device, '--lazy-warmup']
        with (folder / 'daemon.log').open('w') as log:
            daemon = subprocess.Popen(daemon_args, stdout=log, stderr=subprocess.STDOUT, env=env, start_new_session=True)
            try:
                for _ in range(600):
                    if sock.exists():
                        break
                    assert daemon.poll() is None, 'Oracle daemon failed to start.'
                    time.sleep(.1)
                assert sock.exists(), 'Oracle daemon socket timeout.'
                args = ['rfd4-evaluate', '--config-path', str(ROOT / 'configs'), '--config-name', 'paired_evaluation',
                        f'inputs={inputs}', f'out_dir={folder}/results', 'compress_outputs=false', '+metrics.refolding_oracles.oracle_ppi.seed=20261118']
                run_process(args, folder / 'execution', min(available(), 7200), env)
            except BaseException as exc:
                for r in chosen:
                    event('oracle_evaluation', evaluation_id=r['design_id'] + '_of3_failed', design_id=r['design_id'],
                          backend='of3', role='on-target', status='failed', error=repr(exc))
                raise
            finally:
                stop_process(daemon)
        for r in chosen:
            did = r['design_id']
            predicted = list((folder / 'results').rglob(did + '*_of3_ppi_refolded.cif'))
            confidence = list((folder / 'results').rglob(did + '*_of3_ppi_confidences.json'))
            assert len(predicted) == len(confidence) == 1, f'Missing or ambiguous oracle artifacts for {did}'
            eid = did + '_of3'
            event('oracle_evaluation', evaluation_id=eid, design_id=did, backend='of3', role='on-target', status='completed')
            event('oracle_prediction', prediction_id=eid + '_prediction', evaluation_id=eid, cif_path=str(predicted[0]), confidence_path=str(confidence[0]), seed=20261118)
            metrics, geometry = refold_metrics(r, predicted[0], confidence[0])
            r['oracle'] = metrics
            r['refold_geometry'] = geometry
            r['computational_pass'] = bool(r['triage_pass'] and not r['known_reference_similar'] and r['cdr_glycosylation_motifs'] == 0 and
                metrics['refold_geometry_pass'] and metrics['target_aligned_fv_rmsd'] <= 2.5 and metrics['vh_plddt'] >= 80 and
                metrics['vl_plddt'] >= 80 and metrics['target_interface_pass'] and metrics['pairing_iptm'] >= 60)
            event('metrics', design_id=did, prediction_id=eid + '_prediction', evaluator='paired_of3_v1', values=metrics)
            tested.add(r['sequence_hash'])
            selected_history.append(r)
            state['refolded'] += 1
        save()
    def choose(n):
        pool = [r for r in records if r.get('triage_pass') and r.get('sequence_hash') not in tested and not r.get('known_reference_similar')]
        pool.sort(key=lambda r: (sum(x['task'] == r['task'] for x in selected_history), r['cdr_glycosylation_motifs'], -r['cdr_contact_fraction'], r['paired_framework_rmsd']))
        chosen = []
        for r in pool:
            if all(r['sequence_hash'] != x['sequence_hash'] and similarity(cdr_sequence(r), cdr_sequence(x)) < .8 for x in chosen + selected_history):
                chosen.append(r)
            if len(chosen) >= n:
                break
        return chosen
    try:
        # Register completed pilot structures with explicit source metadata before metric validation.
        with ledger.connect(DB) as conn:
            for task in plan['tasks']:
                ledger.plan_run(conn, input_id=task, config_id='pilot_config', phase='pilot_import', run_id='pilot_' + task,
                                command='Import existing pilot with original CSV and Hydra provenance retained under runs/initial_designs', out_dir=str(ROOT / 'runs/initial_designs'))
        for i, row in enumerate(csv.DictReader((ROOT / 'runs/initial_designs/inference_metrics/per_example_metrics_job0.csv').open())):
            add(row['designed_structure_path'], row['task'], 'pilot_' + row['task'], f'pilot_{i:02}', example_id=row['example_id'])
        save()
        initial = []
        for patch_tasks in [plan['tasks'][:2], plan['tasks'][2:]]:
            pool = [r for r in records if r.get('triage_pass') and r['task'] in patch_tasks and not r.get('known_reference_similar')]
            assert pool, 'Each target patch requires a pilot for oracle validation.'
            initial.append(min(pool, key=lambda r: (r['cdr_glycosylation_motifs'], -r['cdr_contact_fraction'])))
        evaluate(initial)
        assert state['refolded'] == 2, 'Two complete paired-oracle pilots required before scaling.'
        for i in range(BATCHES):
            if (ROOT / 'STOP').exists() or (BASE / 'STOP').exists() or time.monotonic() + 2400 > generation_deadline or available() < 3600:
                state['stop_reason'] = 'Generation budget or STOP marker reached.'
                break
            task = plan['tasks'][i % 4]
            run_id = f'batch_{i:02}'
            out = BASE / 'runs' / run_id
            state['stage'] = run_id
            save()
            event('run_status', run_id=run_id, status='running')
            failure = None
            try:
                run_process(command(task, i), out / 'execution', min(available(), 5400))
            except Exception as exc:
                failure = exc
            files = list((out / 'inference_metrics').glob('per_example_metrics*.csv'))
            rows = list(csv.DictReader(files[0].open())) if len(files) == 1 else []
            for j, row in enumerate(rows):
                assert row['task'] == task
                add(row['designed_structure_path'], task, run_id, f'b{i:02}_d{j:03}', seed=SEED_BASE + i, example_id=row['example_id'])
            state['generated'] += len(rows)
            event('run_status', run_id=run_id, status='failed' if failure or len(rows) != DESIGNS_PER_BATCH else 'completed')
            save()
            if failure:
                raise failure
            assert len(rows) == DESIGNS_PER_BATCH, 'Incomplete generation batch; inspect retained artifacts.'
            if state['generated'] >= 160 and not any(r.get('triage_pass') for r in records if not r['design_id'].startswith('pilot')):
                state['stop_reason'] = 'No triage passers among 160 new designs.'
                break
            if state['generated'] % 200 == 0:
                evaluate(choose(min(4, 80 - state['refolded'])))
                if state['refolded'] >= 24 and state['computational_passers'] == 0:
                    state['stop_reason'] = 'Zero computational passers after at least 24 complete refolds; redesign required.'
                    break
        if not state.get('stop_reason'):
            while state['refolded'] < 80 and available() > 2400 and not (ROOT / 'STOP').exists() and not (BASE / 'STOP').exists():
                chosen = choose(min(4, 80 - state['refolded']))
                if not chosen:
                    break
                evaluate(chosen)
                if state['refolded'] >= 24 and state['computational_passers'] == 0:
                    state['stop_reason'] = 'Zero computational passers after at least 24 complete refolds; redesign required.'
                    break
        state['status'] = 'completed_with_stop_condition' if state.get('stop_reason') else 'completed'
        save()
    except BaseException as exc:
        state.update(status='failed', error=repr(exc), traceback=traceback.format_exc())
        save()
        raise


if __name__ == '__main__':
    if '--initialize' in sys.argv:
        initialize()
    else:
        main()
