library(ggplot2)
library(reshape2)
library(tidyverse)

df <- read.csv("results_bincount.txt")
df$Bin.Size = as.factor(df$Bin.Size)
df <- df[,1:3]
df.m = melt(df)
cbPalette <- c("gray0","gray75")
cbPalette <- c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02")
#cbPalette <- c( "#E69F00", "#56B4E9", "#009E73",  "#CC79A7", "darkred", "#F0E442","indianred4", "#999999", "#0072B2","darkslategray" ,"#D55E00",  "indianred4")
print(str(df))
#df.m$variable <- gsub('Zero.')
df.m$value <- as.numeric(df.m$value/1000)
df.m$variable <- gsub('Zero.Counts','GBS bin with zero snps', df.m$variable)
df.m$variable <- gsub('Bin.Count','Total bin count', df.m$variable)

gg <- ggplot(df.m, aes(x=Bin.Size, y=value, fill=variable)) + 
    geom_bar(position = position_dodge(), 
             stat = 'identity', 
             width = 0.8)
gg <- gg + theme_classic()
gg <- gg + theme(axis.title.x=element_text(size=15),
                axis.text.x=element_text(size=12),
                axis.ticks.x=element_blank(),
                axis.title.y=element_text(size=15),
                axis.text.y=element_text(size=12),
                legend.text=element_text(size=15),
                panel.background=element_blank(),
                legend.title=element_blank(),
                legend.position = c(0.6,0.8),
                panel.spacing.x=unit(0,"line"))
gg <- gg + xlab('Bin size') + ylab('Count (thousand)') + scale_fill_manual(values=cbPalette)
    
ggsave(gg, file='Results_figures/Bincount_Barplot.png', units='in', width=4,height=3,dpi=300)

### violin + boxplot for GBLUP results
cbPalette <- c("gray0","gray75", "#56B4E9",  "indianred4", "#E69F00","#0072B2","#F0E442" )
gg <- df2 %>% ggplot(aes(Trait,Predictive_ability, fill=Method)) + geom_violin(width=0.9) + geom_boxplot(width=0.1,color='grey10',alpha=0.2, position=position_dodge(width=0.9))
gg <- gg <- gg + theme_minimal()
gg <- gg + theme(axis.title.x=element_blank(),
                axis.text.x=element_text(size=12),
                axis.title.y=element_text(size=12),
                axis.text.y=element_text(size=12),
                legend.text=element_text(size=12),
                panel.background=element_blank(),
                legend.title=element_blank(),
                panel.spacing.x=unit(0,"line"))
gg <- gg + ylab('Predictive ability (r)')
gg <- gg + scale_fill_manual(values=cbPalette)

ggsave(gg, file='Results_figures/GBLUP_CV10.png', units='in', width=7,height=2,dpi=300)

library(RColorBrewer)
cbPalette <- brewer.pal(6,'Dark2')

gg <- ggplot(pca,aes(PC1,PC2, color=factor(Subpop,levels=unique(Subpop)))) + geom_point(size=3)
gg <- gg + theme_void()
gg <- gg + theme(axis.title.x=element_text(size=20),
                axis.text.x=element_text(size=17),
                axis.ticks.x=element_blank(),
                axis.title.y=element_text(size=20, angle=90),
                axis.text.y=element_text(size=17),
                #legend.text=element_text(size=20),
                legend.position='none',
		panel.background=element_blank(),
                legend.title=element_blank(),
		panel.spacing.x=unit(0,"line"),
		plot.margin = margin(10, 10, 10, 10))

gg <- gg + xlab('PC1 9.36%') + ylab('PC2  7.86%') + guides(color = guide_legend(override.aes = list(size = 5)))
gg <- gg +scale_color_manual(values=cbPalette)
ggsave(gg, file='Results_figures/PC1_PC2_plot.png', units='in', width=4.5,height=4,dpi=300)









