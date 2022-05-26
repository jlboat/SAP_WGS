#!/usr/bin/env bash
#PBS -N vg_align
#PBS -o vg_align.out
#PBS -e vg_align.err
#PBS -l select=1:ncpus=32:mem=64gb,walltime=140:00:00
#PBS -q tillers

cd /scratch1/jboatw2/merge_gvcfs/recal_gvcfs_all/filter/all

vg map --base-name wg --threads 32 --min-mem 16 -f propinquum369-2_1.fastq.gz -f propinquum369-2_2.fastq.gz > propinquum.gam
