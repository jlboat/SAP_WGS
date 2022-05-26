import pandas as pd

tajd = pd.read_csv("/scratch1/jboatw2/SAP/results/gwas/QUALITY/recal/metrics/SAP.tajd.Tajima.D", sep="\t")
tajd["end"] = tajd.BIN_START + 99999
tajd.columns = ["#chrom","start","N_SNPS","value","end"]
tajd = tajd[["#chrom","start","end","value"]]
tajd.to_csv("sorghum_tajd.all.txt", sep="\t",index=False)
