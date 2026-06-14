#!/usr/bin/env bash
set -e
python run_pipeline.py --config config.yaml --skip-foldseek --input-tsv data/example_foldseek_result.tsv
