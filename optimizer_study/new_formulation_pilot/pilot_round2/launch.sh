#!/bin/bash
# Round 2: kick + optimizer, more trials. Trials 1-5 reuse round 1 seeds where the label matches.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
PY=~/qdb-envs/qenv/bin/python
W=~/Downloads/umich-db-quantumdb-wt/formulation
caffeinate -i $PY -W ignore pilot.py $W P2,P4 NFT,COBYLA 20 1000 8 > log_nokick.txt 2>&1
KICKS=0.2,0.35,0.5,0.75,1.0 caffeinate -i $PY -W ignore pilot.py $W P2 NFT,COBYLA 20 1000 8 > log_p2.txt 2>&1
KICKS=0.2,0.35,0.5,0.75,1.0 caffeinate -i $PY -W ignore pilot.py $W P2 NFT 20 3000 8 > log_p2_b3000.txt 2>&1
KICKS=0.2,0.5,1.0 caffeinate -i $PY -W ignore pilot.py $W P4 NFT,COBYLA 10 1000 8 > log_p4.txt 2>&1
KICKS=0.2,0.5,1.0 caffeinate -i $PY -W ignore pilot.py $W P4 NFT 10 3000 8 > log_p4_b3000.txt 2>&1
echo done > round2.done
