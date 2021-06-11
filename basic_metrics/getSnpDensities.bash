#!/usr/bin/env bash

# 01 02 03 06 07 08 09 10
# 
for i in 01 02 03 04 05 06 07 08 09 10
do
    python ~/scripts/snpDensityPlotFromVcf.py --vcf SAP.imputed.all.Chr${i}.vcf.gz --output Chr${i} > chrom${i}_density.txt
done 
