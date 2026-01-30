#!/bin/bash
set -euo pipefail

EXPECTED_TOTAL=144
START_TIME=$(date +%s)

echo "Starting perturbations pipeline @ $(date)"

mkdir -p logs/perturbations_logs
mkdir -p figures

module load miniconda
conda activate spatial_analysis_env

jid_make=$(qsub make_perturbations.qsub | awk '{print $3}' | cut -d. -f1)
echo "Submitted make_perturbations as job $jid_make"

qsub -P <project_id> -hold_jid "$jid_make" -sync y -b y -o /dev/null -e /dev/null -N wait_make true >/dev/null  # Hard wait until job completes

echo "make_perturbations completed."

# Validate index files: 
echo " Validating perturbation index files:"

IDX_8="perturbations_8gb_indices.txt"
IDX_16="perturbations_16gb_indices.txt"
IDX_32="perturbations_32gb_indices.txt"

N_8GB=0
N_16GB=0
N_32GB=0

[[ -f "$IDX_8"  ]] && N_8GB=$(wc -l < "$IDX_8")
[[ -f "$IDX_16" ]] && N_16GB=$(wc -l < "$IDX_16")
[[ -f "$IDX_32" ]] && N_32GB=$(wc -l < "$IDX_32")

TOTAL=$((N_8GB + N_16GB + N_32GB))

echo "8GB tasks : $N_8GB"
echo "16GB tasks: $N_16GB"
echo "32GB tasks: $N_32GB"
echo "TOTAL     : $TOTAL"

if [[ "$TOTAL" -ne "$EXPECTED_TOTAL" ]]; then
    echo "ERROR: Expected $EXPECTED_TOTAL perturbations, found $TOTAL"
    echo "Aborting pipeline."
    exit 1
fi

echo "Submitting perturbation array jobs:"

HOLD_JIDS=()

if [[ "$N_8GB" -gt 0 ]]; then
    jid8=$(qsub -t 1-"$N_8GB" perturbations_8gb.qsub | awk '{print $3}' | cut -d. -f1)
    echo "8GB perturbations submitted as $jid8"
    HOLD_JIDS+=("$jid8")
fi

if [[ "$N_16GB" -gt 0 ]]; then
    jid16=$(qsub -t 1-"$N_16GB" perturbations_16gb.qsub | awk '{print $3}' | cut -d. -f1)
    echo "16GB perturbations submitted as $jid16"
    HOLD_JIDS+=("$jid16")
fi

if [[ "$N_32GB" -gt 0 ]]; then
    jid32=$(qsub -t 1-"$N_32GB" perturbations_32gb.qsub | awk '{print $3}' | cut -d. -f1)
    echo "32GB perturbations submitted as $jid32"
    HOLD_JIDS+=("$jid32")
fi

HOLD_JIDS_CSV=$(IFS=, ; echo "${HOLD_JIDS[*]}")

echo "Submitting merge_parquets.qsub:"

jid_merge=$(qsub -hold_jid "$HOLD_JIDS_CSV" merge_parquets.qsub | awk '{print $3}' | cut -d. -f1)
echo "merge_parquets submitted as $jid_merge"

echo "Submitting stability clustering and plotting:"

jid_jaccard=$(qsub -hold_jid "$jid_merge" neighbor_jaccard_stability.qsub | awk '{print $3}' | cut -d. -f1)
echo "neighbor_jaccard_stability submitted as $jid_jaccard"

jid_plot=$(qsub -hold_jid "$jid_jaccard" plot_stability.qsub | awk '{print $3}' | cut -d. -f1)
echo "plot_stability submitted as $jid_plot"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

REPORT="logs/pipeline_runtime_report.txt"
{
    echo "Pipeline submission report @ $(date)"
    echo "----------------------------------"
    echo "8GB perturbations : $N_8GB"
    echo "16GB perturbations: $N_16GB"
    echo "32GB perturbations: $N_32GB"
    echo "TOTAL             : $TOTAL"
    echo ""
    echo "make_perturbations: $jid_make"
    echo "merge_parquets    : $jid_merge"
    echo "neighbor_jaccard_stability     : $jid_jaccard"
    echo "plot_stability    : $jid_plot"
    echo ""
    echo "Submission time   : ${ELAPSED}s"
} > "$REPORT"

cat "$REPORT"

echo "=========================================="
echo "Pipeline successfully submitted 🎉"
echo "=========================================="
