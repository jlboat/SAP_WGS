
library(data.table)
library(tidyverse)
library(lokern)

setwd('/zfs/tillers/panicle/HDP_data')

## read in snp matrix 0, 1, 2 format and calculate Fst by allele frequencies

sap <- fread('/zfs/tillers/panicle/ssapkot/HDP/data/Genomic/SAP_MAF_0.05_matrix.txt')

sap <- as.matrix(sap[,-1])
sap <- t(sap)
ind <- read.table('/scratch1/ssapkot/Filtered_phased/SAP_WGS_filtered.012.indv', header=F)
colnames(sap) <- ind$V1

## for racial differentiation
#pop <- fread('/zfs/tillers/panicle/ssapkot/HDP/data/Clustering_races_SAP_CSmanuscript.csv')

## for converted vs bred
conv <- read.table('selection/Fst_TajD_old/Fst_types/converted_cuso.txt', header=F)
bred <- read.table('selection/Fst_TajD_old/Fst_types/historic_selected_cuso.txt', header=F)

## make frequency matrix
M=matrix(NA,nrow(sap),2)
M[,1] <- apply(sap[,which(colnames(sap) %in% conv$V1)],1,function(x) sum(x)/(length(x)*2))
M[,2] <- apply(sap[,which(colnames(sap) %in% bred$V1)],1,function(x) sum(x)/(length(x)*2))
colnames(M) <- c('CONV','BRED')

#for (i in 1:6){
#	df <- pop %>% filter(`K-Cluster`==i)
#	M[,i] <- apply(sap[,which(colnames(sap) %in% df$CUso)],1,function(x) sum(x)/(length(x)*2))
#}
#colnames(M) <- c('Durra','Kafir','Caudatum','Guinea','Milo','Mixed')


# average allele frequency across populations
meansB=rowMeans(M)
alleleVar=meansB*(1-meansB) # p * q variance
meanDevB=M-meansB
FST=(meanDevB^2)/alleleVar

print(FST[1:5,])
snp <- read.table('SAP_0.5_MAF_phased_SNPs_positions.txt',header=F)

smoothFst <- snp
for (i in 1:2){
	smooth <- lokerns(FST[,i], n.out=542) #10000 snp window
	smoothFst <- cbind(smoothFst,smooth$est)
}
#colnames(smoothFst) <- c('CHROM','POS','DURRA','KAFIR','CAUDATUM','GUINEA','MILO','MIXED')
colnames(smoothFst) <- c('CHROM','POS','CONV','BRED')

print(smoothFst[1:5,])
df_fst <- gather(smoothFst,"POP","Smoothed_Fst", -c(CHROM,POS))
threshold <- df_fst %>% group_by(POP) %>% summarize(sig=(mean(Smoothed_Fst) + (2*sd(Smoothed_Fst))))

## copy significantly selected regions

#group <- unique(df_fst$POP)
#chrom <- unique(df_fst$CHROM)
#df_sel <- c()
#
#for (i in group){
#	data <- df_fst %>% filter(POP==i)
#	sig <- as.numeric(threshold[threshold$POP==i,'sig'])	
#	data <- data %>% filter(Smoothed_Fst > sig)
#	
#	for (j in chrom){
#		data2 <- data[data$CHROM==j,]
#		if (nrow(data2) > 0){
#		file <- paste0('selection/Fst_byChrom/Fst_Sig_',i,'_',j,'.png')
#		#png(file,units='in',width=10,height=3,res=300)
#		gg <- ggplot(data2, aes(POS,Smoothed_Fst)) + geom_line() + theme_classic()
#		ggsave(file=file)
#		#dev.off()
#		print(head(data2))
#		}
#	}
#}


save(df_fst, M, pop, snp, file='selection/Fst_ConvVsBred_smoothed_10000snp_SAP_5percent.RData')

## all fst by population averages
png('selection/Fst_ConvBred_Smoothed_10000snp.png', units='in',width=10,height=12,res=300)
gg <- ggplot(df_fst,aes(POS,Smoothed_Fst))  + geom_line(color='mediumblue', size=1) + facet_grid(vars(POP),vars(CHROM),scales="free_y",switch='x')
gg <- gg + geom_hline(data=filter(threshold, POP==POP),aes(yintercept=sig), linetype='dashed', col = 'red')
gg <- gg + theme_minimal() #+ scale_color_manual(values = c('black','navyblue','forestgreen'))
gg <- gg + theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.title.y=element_text(size=15),
                axis.text.y=element_text(size=15),
                strip.text.x=element_text(size=17),
                strip.text.y=element_text(size=17),
                panel.background=element_blank(),
                panel.spacing.x=unit(0,"line"))

#gg <- gg + geom_vline(data=filter(gene_df,CHROM==CHROM), aes(xintercept=BIN_START), colour="gray75", linetype="dashed")
#gg <- gg + geom_text_repel(data=filter(gene_df,STAT=='BetweenRaces'),aes(label=Genes),color='navyblue')
#gg <- gg + geom_text_repel(data=filter(all_df,CHROM=='Chr01'),aes(label=letters[1:3]),color='black',font="bold",size=20)
#gg <- gg + ylab("Tajima's D                  Mean Fst                      Mean Fst")
gg
dev.off()



