#!/bin/bash
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
caffeinate -i ~/qdb-envs/senv/bin/python -W ignore spiq_new.py qaoa >> log_spiq.txt 2>&1
caffeinate -i ~/qdb-envs/senv/bin/python -W ignore spiq_new.py q2o >> log_spiq.txt 2>&1
echo SPIQ_ALL_DONE >> log_spiq.txt
