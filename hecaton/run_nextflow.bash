#!/usr/bin/env bash

WORK_DIR="/scratch1/jboatw2/propinquum_hecaton"
CONFIG_FILE="${WORK_DIR}/config/nextflow_pbspro.config"

# If still can't find tillers_v100
# /home/jboatw2/github/nextflow/launch.sh

/home/jboatw2/github/nextflow/launch.sh run -w functional_workdir_no_align -c ${CONFIG_FILE} -resume hecaton.no_align.nf \
    --genome_file ${WORK_DIR}/data/reference/Sbicolor_454_v3.0.1.fa \
    --bwa_bams "${WORK_DIR}/alignments/*.bam" \
    --model_file ${WORK_DIR}/models/random_forest_model_concat_A_thaliana_ColxCvi_O_sativa_Suijing18_coverage_10x_insertions_balanced_subsample.pkl \
    --cutoff 0.7 \
    --manta_config ${WORK_DIR}/config/configManta_weight_1.py.ini \
    --extra_filtering true \
    --output_dir ${WORK_DIR}/results
