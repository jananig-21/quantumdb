#!/bin/bash
# QAOA with the NFT optimizer from an existing single-start SPIQ point (new formulation) + postprocess.
# Usage: run_nft.sh <input_idx> <week> <spiq_week>   (logs to logs/nft_<week>_input<idx>.log)
set -e
IDX="$1"; WEEK="$2"; SPIQ_WEEK="$3"
REPO=/home/repro/sigmod-repro
JSON=$(ls "${REPO}/spiq_init_outputs/single_${SPIQ_WEEK}"/spiq_initial_point_input"${IDX}"_reps2_*_qaoa.json)
echo "=== start $(date): ${JSON}"
cd "${REPO}/base"
python3 IBMQExperiments.py --trial 1 --reps 2 --optimizer 3 --iterations 1000 --input_idx "${IDX}" \
  --week "${WEEK}" --spiq_json "${JSON}" --spiq_starts 1
PCIRC="${REPO}/base/${WEEK}/ExperimentalAnalysis/IBMQ/QPUPerformance/Results/CPU_Data/iterations_1000/reps_2/NFT/input${IDX}/trial1/pcirc.qpy"
[ -f "${PCIRC}" ] || { echo "=== FAILED: ${PCIRC} missing, QAOA did not run from SPIQ $(date)"; exit 1; }
echo "=== qaoa done $(date)"
cd "${REPO}"
python3 base/postprocess_results.py --trial 1 --reps 2 --optimizer NFT --input_idx "${IDX}" --iterations 1000 --week "${WEEK}"
cd "${REPO}/base"
mv "ExperimentalAnalysis/IBMQ/QPUPerformance/Results/CPU_Data/1000_Iterations/${IDX}_predicates-newQUBO/trial1" "${WEEK}/response_pickle_trial1_qaoa"
echo "=== all done $(date)"
