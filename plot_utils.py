import re, gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def get_gene_position(gene_list, gff3_path):
        gene_pos = []
        
        if gene_list:
            gene_set = set([g.strip() for g in gene_list.split(",") if g.strip()])

            with gzip.open(gff3_path, "rt") as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    if not any(f"gene_name={g};" in line for g in gene_set):
                        continue

                    fields = line.strip().split("\t")
                    gene_name = re.search(r"gene_name=([^;]+)", fields[8]).group(1)
                    gene_pos.append({
                        "#chrom": fields[0], 
                        "pos": (int(fields[3]) + int(fields[4])) // 2, 
                        "symbol": gene_name
                    })
        
        return pd.DataFrame(gene_pos)


def plot_solution(segment_list, interval_list, chr_chrom, gene, cen_start, cen_end, outdir, sample):
    chrom = chr_chrom[3:] if chr_chrom.startswith("chr") else chr_chrom
    if not gene.empty:
        gene = gene[(gene["#chrom"] == chr_chrom) | (gene["#chrom"] == chrom)]
    
    for rank, (segment, interval) in enumerate(zip(segment_list, interval_list)):
        max_cn = max(interval["cn"].max(), segment["cn"].max(), 2)
        chrom_end = interval["end"].iloc[-1]
        
        fig, ax = plt.subplots(figsize=((chrom_end / 2e7) + 0.6, (max_cn + 1) / 2 + 0.6))
        
        # Background lines (lightgrey)
        for y in range(0, max_cn + 2):
            ax.axhline(y, color="#f0f0f0", linestyle="-", linewidth=0.4, zorder=0)
        
        if cen_start and cen_end:
            ax.axvspan(cen_start, cen_end, color="#f0f0f0", zorder=0)

        # Gene lines and labels (lightblue)
        for _, g in gene.iterrows():
            ax.axvline(g["pos"], color="#c5e4f3", linestyle="-", linewidth=4, zorder=0)
            ax.text(g["pos"], 1.02, g["symbol"], transform=ax.get_xaxis_transform(), fontsize=10, style="italic", ha="center", va="bottom")

        # Segment rectangles (grey)
        for _, row in segment.iterrows():
            ax.add_patch(plt.Rectangle((row["start"], row["cn"] - 0.4), row["end"] - row["start"], 0.8, facecolor="#d3d3d3", edgecolor="grey", linewidth=1))

        # Interval step lines (red)
        interval = pd.concat([interval, interval.tail(1)], ignore_index=True)
        interval.loc[interval.index[-1], "start"] = interval.loc[interval.index[-2], "end"]
        ax.step(interval["start"], interval["cn"], where="post", color="red", linewidth=1.5)

        # Axes
        xticks = np.arange(0, chrom_end + 1, 10_000_000)
        ax.xaxis.set_major_locator(plt.FixedLocator(xticks))
        xlabels = [f"{int(x*1e-6)}" if x % 20_000_000 == 0 else "" for x in xticks]
        ax.set_xticklabels(xlabels)
        ax.set_xlim(0, chrom_end) 
        ax.set_xlabel(f"Chromosome {chrom} (Mb)")

        ax.set_ylim(-0.5, max_cn + 1.5)
        ax.set_yticks(range(0, max_cn + 2))
        ax.set_ylabel("Copy Number")

        ax.set_title(f"{sample}", fontsize=14, fontweight="bold", y=1.1)
        ax.tick_params(labelsize=10)
        fig.tight_layout()

        # Save
        filename = f"{outdir}/{sample}.{chr_chrom}.solution{rank + 1}.pdf"
        plt.savefig(filename)
        plt.close()
