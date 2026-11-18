#!/bin/bash
# Held-out test: seeds not used in rounds 1-2 (trials 101+), every method at the same 3,000-evaluation budget.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 TRIAL_OFFSET=100
PY=~/qdb-envs/qenv/bin/python
W=~/Downloads/umich-db-quantumdb-wt/formulation
caffeinate -i $PY -W ignore pilot.py $W P2,P4 NFT,COBYLA 20 3000 8 > log_nokick.txt 2>&1
KICKS=0.5,1.0 caffeinate -i $PY -W ignore pilot.py $W P2,P4 NFT,COBYLA 20 3000 8 > log_kicks.txt 2>&1
echo done > round3.done
