#!/usr/bin/env bash

VCF=SAP.imputed.all.vcf.gz
WINDOW=100000

vcftools --gzvcf ${VCF} --out SAP.tajd --TajimaD ${WINDOW}
vcftools --gzvcf ${VCF} --out SAP.nt_diversity.window --window-pi ${WINDOW} # --site-pi
