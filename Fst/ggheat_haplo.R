
library(tidyverse)


## read haploblock matrix
hap <- read.table('selection/filtered_blocks.csv', sep=',',header=T)
hap <- hap[,-1]
colnames(hap)[1] <- 'Block'

ann <- read.csv('../ssapkot/HDP/data/Clustering_races_SAP_CSmanuscript.csv', header=T)
ann <- data.frame(as.factor(ann[,2]))
colnames(ann) <- 'POP'
rownames(ann) <- colnames(hap)

block <- rownames(hap)
block <- data.frame(substr(block, 1,4)) ## cut only chrx part
rownames(block) <- rownames(hap)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",  "#CC79A7", "darkred", "#F0E442","indianred4", "#0072B2","darkslategray" ,"#D55E00",  "indianred4")

chr <- paste0('Chr',c(1:5,8,10))

for (i in chr) {
	hap_n <- hap[which(grepl(i,rownames(hap))),]
	file <- paste0('selection/heatmap_blocks_Subpop_',i,'.png')	
	gg <- pheatmap(hap_n, show_rownames=FALSE,show_colnames=FALSE, annotation_col=ann)
	ggsave(gg, file=file, units="in",width=10,height=10,dpi=300)
}

## subset haploblocks by selected region coordinates
sel <- read.csv('selection/Fst_significant_pop_1Mb.csv', header=T)
coord <- read.csv('selection/filtered_coords.csv')
coord <- coord[,-1]
colnames(coord) <- c('Locus','BIN_START','BIN_END','MID','CHR')
coord$Index <- rownames(coord)

## example: lets pull kafir chr 5 region
sel2 <- sel %>% filter(POP=='kafir')
sel2
kafir5 <- coord %>% filter(CHR==5&BIN_START>39100000&BIN_END<44700000)

hap_k5 <- hap[kafir5$Index,]

file <- paste0('selection/heatmap_blocks_Subpop_Chr2_Guinea.png')
gg <- pheatmap(hap_G2, show_rownames=FALSE,show_colnames=FALSE, annotation_col=ann)
        ggsave(gg, file=file, units="in",width=8,height=3,dpi=300)


