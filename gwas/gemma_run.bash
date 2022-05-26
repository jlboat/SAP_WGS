#!/usr/bin/env bash


# EIGENVALUES="top_three_eigenvalues.txt"
# EIGENVECTORS="top_three_eigenvectors.txt"
PLINK_OUTPUT="pruneddata" # The output name for Plink data
MISSING=1
HWE=0
MAF=0.05
PHENOTYPE=3

# k (specifies relatedness matrix)
# lmm=1 (Wald test) lmm=2 (likelihood ratio test) lmm=3 (score test) lmm=4 (all three tests)
# /panicle/software/gemma-0.98.1 
#for i in {26..33}
#do
    #PHENOTYPE=${i}
    PHENOTYPE=11
    ~/gemma-0.98.3 -bfile ${PLINK_OUTPUT} \
        -lmm 1 -n ${PHENOTYPE} -o "${PLINK_OUTPUT}_indelLengths_${PHENOTYPE}" \
        -k ./output/${PLINK_OUTPUT}.sXX.txt \
        -maf ${MAF} -miss ${MISSING} -hwe ${HWE} #-c maturity.tsv
#done

#   -d ./output/plink.pericarp_testa_cov.eigenD.txt -u ./output/plink.pericarp_testa_cov.eigenU.txt \
