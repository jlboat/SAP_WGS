#!/bin/bash
#PBS -N bwa_<NUMBER>
#PBS -l select=1:ncpus=8:mem=7gb,walltime=48:00:00
#PBS -o bwa_<NUMBER>.out
#PBS -e bwa_<NUMBER>.err

cd $PBS_O_WORKDIR

module add anaconda3/5.1.0-gcc/8.3.1

BASE_DIR="/scratch1/jboatw2/SAP"
OUTPUT_DIR="${BASE_DIR}/results/alignments"
NUM_CPUS=8
GENOME="${BASE_DIR}/data/references/Sbicolor_454_v3.0.1.fa"

DESIGN_FILE="${BASE_DIR}/doc/All_SampleSheet.fixed.csv"
DESIGN=$(cat ${DESIGN_FILE} | head -n <NUMBER> | tail -n 1)

IFS=',' read -ra ARRAY <<< "${DESIGN}"

# RG_SAMPLE_CODE=${ARRAY[0]}
CUSO=${ARRAY[1]}
SAMPLE=$( echo "${ARRAY[4]}" | cut -f 1-2 -d 'R' | sed 's/.$//' )

R1="${BASE_DIR}/results/fastp/${SAMPLE}.R1.trimmed.paired.fq.gz"
R2="${BASE_DIR}/results/fastp/${SAMPLE}.R2.trimmed.paired.fq.gz"

if [ ! -e ${OUTPUT_DIR} ]
then
    mkdir -p ${OUTPUT_DIR}
fi

# load samtools from aligners env
source activate aligners
singularity run -B /scratch1,/zfs ~/singularity_containers/bwa.simg mem -t ${NUM_CPUS} -R "@RG\\tID:${SAMPLE}\\tSM:${CUSO}\\tPL:Illumina" ${GENOME} ${R1} ${R2} | samtools view -t ${NUM_CPUS} -b - > ${OUTPUT_DIR}/${SAMPLE}.bam
