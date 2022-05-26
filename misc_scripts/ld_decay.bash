#!/usr/bin/env bash

singularity run -B /scratch1 ~/singularity_containers/popLdDecay.sif PopLDdecay \
    -InVCF SAP.imputed.MAF_filtered.vcf.gz -MaxDist 600 -MAF 0.1 -OutStat SAP.stat.gz

for i in 01 02 03 04 05 06 07 08 09 10
do
    singularity run -B /scratch1 ~/singularity_containers/popLdDecay.sif PopLDdecay \
        -InVCF SAP.imputed.MAF.Chr${i}.vcf.gz -MaxDist 600 -MAF 0.1 -OutStat SAP.chr${i}.stat.gz
    singularity run -B /scratch1 ~/singularity_containers/popLdDecay.sif /opt/PopLDdecay/bin/Plot_OnePop.pl -inFile SAP.chr${i}.stat.gz -output chr${i}.ld_decay
done
