#!/usr/bin/env python3
"""Run bounded CPU-only preparation checks; never sample or submit a job."""

import argparse
import importlib
import importlib.metadata
import json
import os
from pathlib import Path
import resource
import signal
import subprocess
import sys
import time

WORKSPACE = Path(__file__).resolve().parents[1]
REPO = Path('/storage/group/yjz5549/default/RFD4/RFD4-Proteina')
REQUESTED_ENV = Path('/storage/group/yjz5549/default/RFD4/env')
HF_CACHE = Path('/storage/group/yjz5549/default/RFD4/huggingface')
REVISION = 'a06e2509c7033b010fc43614e9fd840303588502'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config-only', action='store_true', help='Recheck configuration changes without importing model dependencies.')
    args = parser.parse_args()
    started = time.monotonic()
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    os.nice(10)
    resource.setrlimit(resource.RLIMIT_AS, (8 * 1024**3, 8 * 1024**3))
    resource.setrlimit(resource.RLIMIT_CPU, (120, 120))
    signal.alarm(150)
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        os.environ[key] = '1'
    os.environ.update({
        'CUDA_VISIBLE_DEVICES': '', 'NVIDIA_VISIBLE_DEVICES': 'void',
        'HF_HUB_OFFLINE': '1', 'TRANSFORMERS_OFFLINE': '1', 'WANDB_MODE': 'disabled',
        'PYTHONDONTWRITEBYTECODE': '1', 'PROJECT_ROOT': str(REPO),
        'DATA_PATH': str(REPO / 'data_path'), 'HF_HOME': str(HF_CACHE),
        'MPLCONFIGDIR': str(WORKSPACE / 'cache' / 'matplotlib'),
    })
    sys.dont_write_bytecode = True
    os.chdir(REPO)
    sys.path.insert(0, str(REPO / 'src'))
    report = {'status': 'running', 'gpu_enabled': False, 'jobs_submitted': 0, 'checks': {}}
    checks = report['checks']

    def check(name, function):
        try:
            checks[name] = {'status': 'passed', 'result': function()}
        except Exception as exc:
            checks[name] = {'status': 'failed', 'error': f'{type(exc).__name__}: {exc}'}

    def environment():
        names = ('rfd4-proteina', 'torch', 'atomworks', 'biotite', 'hydra-core', 'anarcii')
        code = ('import importlib.metadata as m,json,sys; '
                'names={"rfd4-proteina","torch","atomworks","biotite","hydra-core","anarcii"}; '
                'print(json.dumps({"prefix":sys.prefix,"packages":'
                '{d.metadata["Name"]:d.version for d in m.distributions() '
                'if d.metadata.get("Name","").lower() in names}}))')
        requested = subprocess.run([str(REQUESTED_ENV / 'bin/python'), '-I', '-c', code],
                                   capture_output=True, text=True, check=True, timeout=20)
        return {'active_prefix': sys.prefix, 'active_executable': sys.executable,
                'active_versions': {name: importlib.metadata.version(name) for name in names},
                'requested_environment': json.loads(requested.stdout)}

    def assets():
        root = HF_CACHE / 'hub/models--nvidia--RFD4-Proteina/snapshots' / REVISION
        paths = [root / 'config.json', root / 'model.safetensors',
                 REPO / 'data_path/openfold3/runtime/of3_af3_refolding_env-cuda12.sif',
                 REPO / 'data_path/openfold3/weights/of3_ported_weights.bin.zst']
        records = [{'path': str(p), 'exists': p.is_file(),
                    'bytes': p.stat().st_size if p.is_file() else None} for p in paths]
        if not all(p['exists'] for p in records):
            raise RuntimeError(json.dumps(records))
        return {'files': records, 'verification': 'File presence and size only; weights were not loaded.'}

    def compose_configs():
        from hydra import compose, initialize_config_dir
        from omegaconf import OmegaConf
        summaries = {}
        for name in ('pilot_generation', 'paired_evaluation'):
            with initialize_config_dir(config_dir=str(WORKSPACE / 'configs'), version_base=None):
                cfg = compose(config_name=name)
            assert OmegaConf.is_missing(cfg, 'inputs'), 'Input must remain mandatory.'
            assert OmegaConf.to_container(cfg, resolve=False)['out_dir'] == '${oc.env:IL4RA_RUN_DIR}', 'Output must require an explicit run directory.'
            summaries[name] = {'mandatory_input': 'inputs', 'required_output_environment': 'IL4RA_RUN_DIR',
                               'metrics_tags': list(cfg.metrics_tags)}
            if name == 'pilot_generation':
                assert 'metrics' not in cfg or cfg.metrics is None
                assert cfg.batch_size == cfg.num_replicates == 1
                assert cfg.num_workers == 0 and cfg.compile_lengths is None
                assert cfg.search.checkpoints is None
                assert sum(n * depth for n, depth in cfg.generation.n_recycle_schedule) == 180
                summaries[name]['trunk_passes'] = 180
            else:
                assert not cfg.metrics.get('pass_criteria')
                assert not cfg.metrics.get('inline_metrics')
                assert not cfg.metrics.get('custom_metrics')
                assert set(cfg.metrics.refolding_oracles) == {'oracle_ppi'}
                assert len(cfg.metrics.refolding_metrics) == 3
                summaries[name]['oracle'] = cfg.metrics.refolding_oracles.oracle_ppi.oracle_name
                summaries[name]['template_query'] = cfg.metrics.refolding_oracles.oracle_ppi.template_query
            (WORKSPACE / 'reports' / (name + '.composed.yaml')).write_text(OmegaConf.to_yaml(cfg))
        return summaries

    def imports():
        modules = ('torch', 'atomworks', 'rfproteina.conditions', 'rfproteina.inference_utils.condition_spec')
        loaded = {name: getattr(importlib.import_module(name), '__file__', None) for name in modules}
        import torch
        assert not torch.cuda.is_initialized()
        return {'modules': loaded, 'cuda_initialized': torch.cuda.is_initialized()}

    def confidence_counterexample():
        from rfproteina.metrics.oracles.oracle_metric import extract_confidence_metrics
        # Synthetic confidence values only: a well-paired Fv with poor target interfaces.
        summary = {'chain_pair_pae_min': [[0, 0.5, 18], [0.5, 0, 20], [18, 20, 0]]}
        result = extract_confidence_metrics(summary, 'H', ['H', 'L', 'T'])
        assert result['min_ipae_coi'] == 0.5
        assert result['min_ipae_H_T'] == 18 and result['min_ipae_L_T'] == 20
        return {'synthetic_fixture': True, 'metrics': result,
                'conclusion': 'The COI minimum can select VH-VL and is unsafe as a target-binding criterion.'}

    if not args.config_only:
        check('environment_metadata', environment)
        check('cached_assets', assets)
    check('hydra_composition', compose_configs)
    if not args.config_only:
        check('cpu_imports', imports)
        check('paired_confidence_counterexample', confidence_counterexample)
    report.update({
        'status': 'passed' if all(c['status'] == 'passed' for c in checks.values()) else 'incomplete',
        'elapsed_seconds': round(time.monotonic() - started, 2),
        'peak_rss_mib': round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 2),
        'cpu_affinity': sorted(os.sched_getaffinity(0)),
        'scope': 'CPU preparation only; not GPU preflight, input round-trip, inference, or oracle validation.',
    })
    output = WORKSPACE / 'reports' / ('config_checks.json' if args.config_only else 'cpu_checks.json')
    output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    return 0 if report['status'] == 'passed' else 1


if __name__ == '__main__':
    (WORKSPACE / 'reports').mkdir(parents=True, exist_ok=True)
    raise SystemExit(main())
