#!/bin/bash
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
exec caffeinate -i ~/qdb-envs/qenv/bin/python -W ignore p5_kick.py ~/Downloads/umich-db-quantumdb-wt/formulation 4 4 >> log_p5_kick.txt 2>&1
