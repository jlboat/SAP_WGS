import pandas as pd

df = pd.read_csv("sorghum_pos_sel.all.txt", sep="\t")
# sig_hits = df[df.value >= 3.215e-10]
sig_hits = df[df.value >= 2]
chroms = set(df["#chrom"])
new_df = pd.DataFrame() 
for chrom in chroms:
    chr_temp = sig_hits[sig_hits["#chrom"] == chrom]
    sorted_hits = chr_temp.sort_values("start")
    sorted_hits["group"] = (sorted_hits["start"]>sorted_hits["end"].shift()).cumsum()
    result = sorted_hits.groupby("group").agg({"#chrom": "first","start":"min", "end": "max"})
    new_df = new_df.append(result)
result = new_df.sort_values(["#chrom","start"])
result["markerName"] = result["#chrom"] + "_" + result["start"].astype(str)
result["position"] = result["start"]
result = result[["markerName","#chrom","position","start","end"]]
result.to_csv("pos_selection.range.csv", index=False)
