library(CMplot)

args <- commandArgs(trailingOnly=TRUE)

if (is.na(args[1])){
    print("ERROR: No file included.")
    print("Rscript manhattan.R gemma.assoc.txt")
    quit("no")
}

run_name <- gsub(".assoc.txt","", args[1])
df <- read.table(args[1], header=TRUE)
df <- df[,c("rs","chr","ps","p_wald")]
colnames(df) <- c("rs","chr","ps",run_name)

CMplot(df,type="p",plot.type="m",LOG10=TRUE,threshold=0.05/nrow(df),file="jpg",memo="",dpi=600,file.output=TRUE,verbose=TRUE,width=14,height=6,chr.labels.angle=45)

CMplot(df,plot.type="q",box=FALSE,file="jpg",memo="",dpi=600,conf.int=TRUE,conf.int.col=NULL,threshold.col="red",threshold.lty=2,file.output=TRUE,verbose=TRUE,width=5,height=5)

#print(nrow(df))
print(0.05/nrow(df))
