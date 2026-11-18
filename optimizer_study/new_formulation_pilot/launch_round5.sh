#!/bin/bash
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 CONFIRM=1
exec caffeinate -i ~/qdb-envs/qenv/bin/python -W ignore restarts.py ~/Downloads/umich-db-quantumdb-wt/formulation 8 >> log_round5.txt 2>&1
