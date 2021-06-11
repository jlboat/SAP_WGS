#!/usr/bin/env bash


VCF="SAP.imputed.final.MAF_filtered.vcf.gz"
PLINK_OUTPUT="final_data" # The output name for Plink data
#PHENOTYPES="mean_ama_BLUP_tan.txt" # tab-separate file with header = IID\tFID\tPhenotype1\tPhenotype2\t...PhenotypeN
PHENOTYPES="Tannins.tsv"
MISSING=1
HWE=0
MAF=0.05

TMPDIR=.

# BUILD BASIC GEMMA FILES
~/Plink/plink \
  --allow-no-sex \
  --make-bed \
  --out ${PLINK_OUTPUT} \
  --pheno ${PHENOTYPES} \
  --vcf ${VCF}

# BUILD PED/MAP FOR PRUNING
# ~/Plink/plink \
#     --out ${PLINK_OUTPUT} \
#     --recode \
#     --vcf ${VCF}

# PRUNE DATA
# ~/Plink/plink \
#     --file ${PLINK_OUTPUT} \
#     --indep 50 5 2

# gk=1 (centered relatedness matrix) gk=2 (standardized relatedness matrix)
/zfs/tillers/panicle/software/gemma-0.98.3 -bfile ${PLINK_OUTPUT} -gk 2 \
   -o ${PLINK_OUTPUT} -miss ${MISSING} -hwe ${HWE} -maf ${MAF}

# Run gemma eigenvalue decomposition
/zfs/tillers/panicle/software/gemma-0.98.3 -bfile ${PLINK_OUTPUT} \
    -k ./output/${PLINK_OUTPUT}.sXX.txt \
    -eigen -o ${PLINK_OUTPUT} \
    -miss ${MISSING} -hwe ${HWE} -maf ${MAF}

# FIX HERE!
#awk -F ',' '{print $1,($2+$3)/2}' all_Tannin.csv > mean_tan_13_14.txt
#sed 's/PI/pi/g' mean_ama_BLUP_tan.txt > mid
#mv mid mean_ama_BLUP_tan.txt
#python fixFamMultivariate.py mean_ama_BLUP_tan.txt plink.ama_tan_BLUP.fam > mid
#mv mid plink.ama_tan_BLUP.fam

# cd output
# cut -f 2418-2420 plink.pericarp_testa_cov.eigenU.txt > top_three_eigenvectors.txt
# tail -n 3 plink.pericarp_testa_cov.eigenD.txt > top_three_eigenvalues.txt
