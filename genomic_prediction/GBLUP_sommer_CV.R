library(reshape2)
library(tidyverse)
library(sommer)
library(BGLR)
library(BMTME)

### read in phenotypic data and center data by means within a rep of a year
pheno <- read.csv('/zfs/tillers/panicle/HDP_data/ALL_SAP_BLUPS.CuSo.csv', header=T)
colnames(pheno)[1] <- 'CUSo_Geno'
pheno <- pheno[,-2]

## GRM
G <- read.table('/zfs/tillers/panicle/HDP_data/VCFs/recalibrated/SAP_GBS_G-matrix.txt',header=T)
rownames(G) <- colnames(G)

geno_id <- rownames(G)
## subset phenotypic data to include only the lines with genomic data (for GBLUP)
pheno <- pheno[pheno$CUSo_Geno %in% geno_id,]
pheno$CUSo_Geno <- as.factor(as.character(pheno$CUSo_Geno)) ##to change values in levels(pheno$CUSo_Geno)

K <- G[sort(geno_id),sort(geno_id)]
K <- as.matrix(K)
pheno <- pheno[order(pheno$CUSo_Geno),]

traits <- colnames(pheno)[-c(1,4,6,10,14,16:21,25)]
Pred_value <- c()

for (i in 123:222){
	
	## make cross validation folds
        CVlist = CV.KFold(data.frame(pheno[,1]), K = 10, set_seed = i)
        folds = CVlist$CrossValidation_list
	accuracy <- c()
	
	for (j in 1:length(traits)){
		trait <- traits[j]
		df <- pheno[,c(colnames(pheno)[1],trait)]
		colnames(df) <- c('CUSo_Geno','trait')
	        ## make cross validation folds
	        CVlist = CV.KFold(data.frame(pheno[,1]), K = 10, set_seed = i)
	        folds = CVlist$CrossValidation_list
	
	        for(k in 1:10){
	          	
			pred <- c()
			ability <- c()
			test_id = df$CUSo_Geno[folds[[k]]]
		
		        #Make training (TRN) and testing (TST) dfs
		  	df$y <- df$trait
		  	df$y[which(df$CUSo_Geno %in% test_id)] <- NA #Mask yields for validation set
		        
			## gblup model
		  	ans <- mmer(y~1,
			random= ~ vs(CUSo_Geno,Gu=K), #+ vs(Year:CUSo_Geno,Gu=EA) + vs(Year:Rep:BlockA),
			#random= ~ vs(idd,Gu=Kd) + vs(Year:CUSo_Geno,Gu=EA) + vs(Year:Rep:BlockA),
			rcov=~units,
			data=df)
		
			gblup <- predict(ans, classify="CUSo_Geno")$pvals   
			ability <- cor(gblup$predicted.value[which(gblup$CUSo_Geno %in% test_id)],df$trait[which(df$CUSo_Geno %in% test_id)], use='complete')
			pred <- cbind(trait,'GBS_CV10',ability)
			
			accuracy  <- rbind(accuracy,pred)
	    		}
		}
	    Pred_value <- rbind(Pred_value,accuracy)
	}

Pred_value <- data.frame(Pred_value)
colnames(Pred_value) <- c('Trait','Method','Predictive_ability')
Pred_value$Predictive_ability <- as.numeric(Pred_value$Predictive_ability)

write.csv(Pred_value, file='GBLUP_CV10_GBS_prediction.csv', quote=F,row.names=F)
