library(RColorBrewer)

args <- commandArgs(trailingOnly=TRUE)

if (is.na(args[1]) | is.na(args[2])){
    print("ERROR: No file included.")
    print("Rscript sort_admixture.R admixture_log.out line_names.txt")
    quit("no")
}

# READ DATA
FILENAME <- args[1]
LINE_NAMES_FILE <- args[2]
CUTOFF <- 0.5000001 # Recommend 0.6 or simple majority

df <- read.table(FILENAME, header=F)
samples <- read.table(LINE_NAMES_FILE, header=F)
rownames(df) <- samples$V1

sort_admixture <- function(df, col){
    return(df[order(df[col], decreasing=T),])
}

dimensions <- c(dim(df)[1])
columns <- colnames(df)
sorted_df <- data.frame()

group_numbers <- c()
for (column in 1:length(columns)){
    df <- sort_admixture(df, column)
    if (column != length(columns)){
        cut_df <- df[df[column] < CUTOFF,]
        dimensions <- c(dimensions, dim(cut_df)[1])
        slice_df <- df[1:(dimensions[column] - dimensions[column + 1]),]
        group_numbers <- c(group_numbers, rep(column, nrow(slice_df)))
        sorted_df <- rbind(sorted_df, slice_df) 
        df <- cut_df
    } else {
        group_numbers <- c(group_numbers, rep(column, dim(samples)[1]-length(group_numbers)))
    }
}
sorted_df <- rbind(sorted_df, df)
# PLOT RESULTS
pdf("admixture.pdf")
barplot(t(as.matrix(sorted_df)), 
        col=brewer.pal(n=length(columns), name="Set1"),
        # col=rainbow(length(columns)), 
        xlab="Individual #", 
        ylab="Ancestry", 
        border=NA,
        space = 0
)
dev.off()
sorted_df$group <- group_numbers
write.csv(sorted_df[,c(1, ncol(sorted_df))], "admix.plot_data.csv")
