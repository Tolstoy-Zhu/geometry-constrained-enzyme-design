"""Verify the five pinned public AF2-PTM parameter files and record local paths."""
import argparse
import hashlib
import json
from pathlib import Path

EXPECTED = {
    'params_model_1_ptm.npz': '5e564f79af5bcd54ccef6e2a6bb0ff01015d01650ebc41d4575e35f0de9ecc84',
    'params_model_2_ptm.npz': '23645d9a82c4af2ed54cd48a7b3c1c2575dc6aa9fe931adb4d7203ca5f0dc398',
    'params_model_3_ptm.npz': 'b53a724294cdd46a45d530b306bced8e429ffee74a4f9a2a897551b30966c415',
    'params_model_4_ptm.npz': '07d0e285ad1406ce3d91a2f6633a8ce9ee99155a53d24ee7783f8fd1e89f5d3a',
    'params_model_5_ptm.npz': '2f8d0bd193e7d6e837958ba92acdc3f02b23b61bf841fa043fc076ee0e1996a5',
}


def verify(data):
    data = Path(data).resolve()
    for name, expected in EXPECTED.items():
        digest = hashlib.sha256()
        with (data / 'params' / name).open('rb') as handle:
            for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b''):
                digest.update(chunk)
        if digest.hexdigest() != expected:
            raise ValueError('AlphaFold parameter checksum mismatch: ' + name)
    receipt = data / 'parameter-provenance.json'
    if receipt.exists():
        existing = json.loads(receipt.read_text())
        if existing.get('parameter_files') != EXPECTED:
            raise ValueError('Existing parameter provenance differs; retain it and review the change')
    else:
        receipt.write_text(json.dumps({
            'source_url': 'https://storage.googleapis.com/alphafold/alphafold_params_2021-07-14.tar',
            'download_method': 'ColabFold official downloader or supplied parameter files',
            'archive_hash_verified': False,
            'parameter_directory': str(data / 'params'),
            'parameter_files': EXPECTED,
            'member_hashes_verified': True,
            'license': 'CC-BY-4.0 (AlphaFold parameter release)',
        }, indent=2), encoding='utf-8')
    print('Five AF2-PTM parameter checksums verified.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', type=Path, required=True)
    verify(parser.parse_args().data)
