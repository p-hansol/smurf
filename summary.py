import argparse, os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import MaxNLocator
import sv_utils


def main():
    parser = argparse.ArgumentParser(description="Summarize segment sizes and copy numbers.")

    parser.add_argument("-i", "--input", required=True, 
                        help="Path to decision.tsv file.")
    parser.add_argument("-c", "--chrA", required=True, 
                        help="Chromosome of interest (e.g., chr2).")
    parser.add_argument("-o", "--outdir", default=os.getcwd(), 
                        help="Output directory (default: current directory).")
    
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    
    # -------------------------------
    # Load segment data
    # -------------------------------
    decision = pd.read_csv(args.input, sep="\t", dtype=str)
    solutions = []
    for _, row in decision.iterrows():
        chrom = row["#chrom"]
        chr_chrom = sv_utils.add_chr_prefix(chrom)
        solution_num = f"solution{row['solution']}"
        solution_path = args.input.replace("decision.tsv", f"{chr_chrom}.{solution_num}.tsv")
        solution = pd.read_csv(solution_path, sep="\t", dtype={"#chrom": str, "start": int, "end": int, "copy_number": int})
        solutions.append(solution)
    segment = pd.concat(solutions, ignore_index=True)

    # -------------------------------
    # Save summary table
    # -------------------------------
    segment["#chrom"] = segment["#chrom"].apply(sv_utils.add_chr_prefix)
    chr_of_interest = sv_utils.add_chr_prefix(args.chrA)

    bin_edges = [0,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000]
    bin_labels = ["1-10","10-100","100-1K","1K-10K","10K-100K","100K-1M","1M-10M","10M-100M",">100M"]

    chrom_levels = ["Chromosome of interest","Other chromosomes"]
    copy_number_levels = ["0","1","2","3+"]
    row_index = pd.MultiIndex.from_product([chrom_levels, copy_number_levels], names=["#chrom","copy_number"])

    segment["#chrom"] = np.where(segment["#chrom"]==chr_of_interest, chrom_levels[0], chrom_levels[1])
    segment["size"] = segment["end"] - segment["start"]
    segment["copy_number"] = segment["copy_number"].map(lambda x: "0" if x==0 else "1" if x==1 else "2" if x==2 else "3+")
    segment["bin"] = pd.cut(segment["size"], bins=bin_edges, labels=bin_labels, right=True)

    summary = segment.groupby(["#chrom","copy_number","bin"]).size().unstack(fill_value=0)
    summary = summary.reindex(index=row_index, columns=bin_labels, fill_value=0)
    summary_path = os.path.join(args.outdir, os.path.basename(args.input).replace("decision.tsv","summary.tsv"))
    summary.to_csv(summary_path, sep="\t", index=False, header=False)

    # -------------------------------
    # Plot barplot
    # -------------------------------
    n_rows, n_cols = summary.shape
    flat = summary.values.flatten()
    x = np.arange(len(flat))
    row_ids = np.repeat(np.arange(n_rows), n_cols)

    fig, ax = plt.subplots(figsize=(14,4))
    colors = plt.cm.tab20c(np.linspace(0,1,n_rows))
    
    for r in range(n_rows):
        ax.bar(x[row_ids==r], flat[row_ids==r], color=colors[r])

    ax.set_xticks(x)
    ax.set_xticklabels(bin_labels * n_rows, rotation=90)
    ax.set_xlim(-0.5, len(flat)-0.5)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x):3d}"))

    ymax = ax.get_ylim()[1]
    box_height = ymax * 0.14

    for i, copy_number in enumerate(copy_number_levels * 2):
        rect = patches.Rectangle((i*n_cols-0.5, ymax), n_cols, box_height,
                                 linewidth=1, edgecolor="black", facecolor=colors[i],
                                 transform=ax.transData, clip_on=False)
        ax.add_patch(rect)
        ax.text(i*n_cols + n_cols/2 - 0.5, ymax + box_height/2, copy_number,
                ha="center", va="center", fontsize=9, fontweight="bold")

    for j, chrom in enumerate(chrom_levels):
        rect = patches.Rectangle((j*n_cols*n_rows/2 - 0.5, ymax + box_height),
                                 n_cols*n_rows/2, box_height,
                                 linewidth=1, edgecolor="black", facecolor="white",
                                 transform=ax.transData, clip_on=False)
        ax.add_patch(rect)
        ax.text((2*j + 1)*n_cols*n_rows/4 - 0.5, ymax + 1.5*box_height,
                chrom, ha="center", va="center", fontsize=10, fontweight="bold")
    
    sample = os.path.basename(args.input).replace(".decision.tsv","")
    fig.suptitle(sample, fontsize=11, fontweight="bold", y=0.85)
    plt.tight_layout()
    plt.savefig(summary_path.replace(".tsv", ".pdf"), bbox_inches="tight")
    plt.close()
    print(sample)

if __name__ == "__main__":
    main()

