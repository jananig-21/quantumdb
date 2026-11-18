#!/bin/bash
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
until grep -q SPIQ_ALL_DONE log_spiq.txt 2>/dev/null; do sleep 30; done
exec caffeinate -i ~/qdb-envs/qenv/bin/python -W ignore run_new.py 2 >> log_new.txt 2>&1
