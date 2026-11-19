"""Project-local AF3 adapter for crystal inputs with incomplete sidechains."""
import os
import shutil
from pathlib import Path
import numpy as np
from rfproteina.metrics.oracles.af3_oracle import AF3Oracle


class PairedAF3Oracle(AF3Oracle):
    @staticmethod
    def _normalize_af3_output(input_array, af3_array):
        # The upstream sidechain-trimming pass compares residue ordinal i in
        # both arrays. AF3 may emit chains in a different order from its input.
        # Align chain order before trimming; subsequent upstream atom mapping
        # still verifies residue identities and exact atom sets.
        chains, starts = np.unique(input_array.chain_id, return_index=True)
        order = chains[np.argsort(starts)]
        if set(order) != set(af3_array.chain_id):
            raise ValueError('AF3 output chain set differs from input.')
        indices = np.concatenate([np.flatnonzero(af3_array.chain_id == c) for c in order])
        reordered = af3_array[indices]
        return AF3Oracle._normalize_af3_output(input_array, reordered)

    def _parse_af3_predicted_structure(self, af3_output_dir, job_name):
        archive = os.environ.get('IL4RA_RAW_ORACLE_ARCHIVE')
        if archive:
            folder = Path(archive) / job_name
            folder.mkdir(parents=True, exist_ok=True)
            for suffix in ['_model.cif', '_summary_confidences.json']:
                source = Path(af3_output_dir) / (job_name + suffix)
                if source.exists():
                    shutil.copy2(source, folder / source.name)
        return super()._parse_af3_predicted_structure(af3_output_dir, job_name)
