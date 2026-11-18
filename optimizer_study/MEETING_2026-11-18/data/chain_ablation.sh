#!/bin/bash
cd "$(dirname "$0")/../scripts"
until grep -q "ALL DONE" ../data/log_ablation.txt; do sleep 60; done
OMP_NUM_THREADS=1 caffeinate -i ~/qdb-envs/qenv/bin/python -W ignore ablation_random.py 3 > ../data/log_ablation_p4.txt 2>&1
