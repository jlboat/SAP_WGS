setwd('/zfs/tillers/panicle/HDP_data')

source('scripts/LC_functions.R')

library(data.table)
library(tidyverse)
library(progress)
library(R.utils)


vcf <- fread('VCFs/recalibrated/SAP.imputed.MAF_filtered.vcf.gz')
colnames(vcf)[1] <- 'CHROM'
vcf <- data.frame(vcf)
vcf$CHROM <- gsub('Chr0','Chr',vcf$CHROM)

#count_df <- c()
#freq_df <- c()
He <- c()

## compute the genotype count, allele frequencies and expected heterozygosity for each marker
for (i in 1:nrow(vcf)){
        #count_df <- rbind(count_df,count.genotypes(vcf[i,10:ncol(vcf)]))
        #freq_df <- rbind(freq_df, allele.freq(count_df))
        He <- rbind(He, expected.het(vcf[i,10:ncol(vcf)]))
print(i)
}

colnames(He) <- c('p','n','Hexp')
He <- cbind(vcf[,1:2],He)

rm(vcf)

save.image('Results_selection/Expected_Het_WholeSAP.RData')
