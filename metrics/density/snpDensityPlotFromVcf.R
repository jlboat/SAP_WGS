library(CMplot)

args <- commandArgs(trailingOnly=TRUE)

if (is.na(args[1])){
    print("ERROR: No file included.")
    print("Rscript snpDensityPlot.R vcf_output.txt")
    quit("no")
}

df <- read.table(args[1], header=F)
colnames(df) <- c("SNP", "chr", "pos", "Chr01")

CMplot(df, type="p", plot.type="d", bin.size=500000, 
       chr.den.col=c("darkgreen","yellow","red"),
       file="jpg", memo="", dpi=300, file.output=TRUE, verbose=TRUE,
       width=9, height=6)

