import argparse, os, glob
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
from scipy.ndimage import uniform_filter
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns


def main():
    parser = argparse.ArgumentParser(description="Cluster genomic rearrangements.")

    parser.add_argument("-i", "--indir", required=True, 
                        help="Input directory containing summary.tsv files.")
    parser.add_argument("-o", "--outdir", default=os.getcwd(), 
                        help="Output directory (default: current directory).")
    parser.add_argument("-w", "--weight", type=float, default=0.4, 
                        help="Weight for the within-sample MS-SSIM component (float between 0 and 1).")
    parser.add_argument("-n", "--n_cluster", type=int, default=5, 
                        help="Number of clusters")

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # -------------------------------
    # Config
    # -------------------------------

    WEIGHT_WITHIN = args.weight 
    WEIGHT_ACROSS = 1 - WEIGHT_WITHIN
    WIN_SIZES = [(1,1), (1,3), (1,5), (1,7), (3,3), (3,5), (3,7)]
    N_CLUSTERS = args.n_cluster
    INPUT_DIR = args.indir + "/*.summary.tsv"
    OUTPUT_DIR = args.outdir

    # -------------------------------
    # Input images
    # -------------------------------
    input_files = sorted(glob.glob(INPUT_DIR))
    samples, images, padded_images = [], [], []

    for f in input_files:
        sample = os.path.basename(f).split('.')[0]
        samples.append(sample)
        
        i = np.loadtxt(f) 
        images.append(i)
        
        pad = np.zeros((4, i.shape[1]), dtype=i.dtype)
        padded_i = np.pad(
            np.vstack((i[0:4, :], pad, i[4:8, :])), 
            pad_width=((2,2),(2,2)), 
            mode="constant", 
            constant_values=0)
        padded_images.append(padded_i)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # -------------------------------
    # Pairwise similarity
    # -------------------------------
    def ssim_rect(a, b, win_size=(3,3), C1=0.01**2, C2=0.03**2):
        """
        Compute the Structural Similarity Index (SSIM) between two 2D arrays `a` and `b`
        using a rectangular window of size `win_size`. Non-square (rectangular) windows
        are supported to handle images with differing height and width scales.
        """
        mu_a = uniform_filter(a, size=win_size)
        mu_b = uniform_filter(b, size=win_size)
        mu_a2, mu_b2, mu_ab = mu_a*mu_a, mu_b*mu_b, mu_a*mu_b
        sigma_a2 = uniform_filter(a*a, size=win_size) - mu_a2
        sigma_b2 = uniform_filter(b*b, size=win_size) - mu_b2
        sigma_ab = uniform_filter(a*b, size=win_size) - mu_ab
        ssim_map = ((2*mu_ab + C1)*(2*sigma_ab + C2)) / ((mu_a2 + mu_b2 + C1)*(sigma_a2 + sigma_b2 + C2))
        return ssim_map.mean()

    def ms_ssim_one_channel(a, b, win_sizes=WIN_SIZES):
        """
        Compute multi-scale SSIM for a single-channel (2D) input using multiple rectangular window sizes.
        Each window is weighted equally by default.
        """
        weights = [1/len(win_sizes)]*len(win_sizes)
        return sum(ssim_rect(a, b, win_size=ws)*wgt for ws, wgt in zip(win_sizes, weights))

    def ms_ssim_two_channel(a, b, weight0=WEIGHT_WITHIN, weight1=WEIGHT_ACROSS, win_sizes=WIN_SIZES):
        """
        Compute multi-scale SSIM combining two different normalization schemes:
            1. Within-sample normalization (proportions) emphasizes relative distribution within each image.
            2. Across-sample normalization (log10 counts) emphasizes absolute magnitude differences between images.
        The final score is a weighted sum of these two components.
        """
        a_within, b_within = a/a.sum(), b/b.sum()
        a_across, b_across = np.log10(a+1), np.log10(b+1)
        s0 = ms_ssim_one_channel(a_within, b_within, win_sizes)
        s1 = ms_ssim_one_channel(a_across, b_across, win_sizes)
        return weight0*s0 + weight1*s1

    n_samples = len(padded_images)
    sim_matrix = np.zeros((n_samples, n_samples))
    for i in range(n_samples):
        for j in range(i, n_samples):
            score = ms_ssim_two_channel(padded_images[i], padded_images[j])
            sim_matrix[i,j] = sim_matrix[j,i] = score

    dist_matrix = 1 - sim_matrix
    np.fill_diagonal(dist_matrix, 0)

    # -------------------------------
    # Dendrogram
    # -------------------------------
    def compute_node_average(current_node, linkage_matrix, leaf_sums, n_samples):
        """
        Recursively compute the average sum of all leaves under a given node:
            1. If the node is a leaf: return the sum of the corresponding sample.
            2. If the node is internal: compute the average of child nodes.
        """
        if current_node < n_samples:
            return leaf_sums[current_node], 1
        else:
            linkage_index = current_node - n_samples
            left_child = int(linkage_matrix[linkage_index, 0])
            right_child = int(linkage_matrix[linkage_index, 1])
            left_sum, left_count = compute_node_average(left_child, linkage_matrix, leaf_sums, n_samples)
            right_sum, right_count = compute_node_average(right_child, linkage_matrix, leaf_sums, n_samples)
            return left_sum + right_sum, left_count + right_count

    def reorder_linkage_matrix(linkage_matrix, leaf_sums, n_samples):
        """
        Reorder a hierarchical clustering linkage matrix so that clusters with larger average sums
        appear on the right side of the dendrogram.
        """
        reordered_matrix = linkage_matrix.copy()
        for i in range(reordered_matrix.shape[0]):
            left_child = int(reordered_matrix[i, 0])
            right_child = int(reordered_matrix[i, 1])
            left_sum, left_count = compute_node_average(left_child, reordered_matrix, leaf_sums, n_samples)
            right_sum, right_count = compute_node_average(right_child, reordered_matrix, leaf_sums, n_samples)
            left_avg = left_sum / left_count
            right_avg = right_sum / right_count
            if left_avg > right_avg:
                reordered_matrix[i, 0], reordered_matrix[i, 1] = reordered_matrix[i, 1], reordered_matrix[i, 0]
        return reordered_matrix

    weights = np.array([1.03, 1.02, 1, 1.01, 1.03, 1.02, 1, 1.01])
    leaf_sums = []
    for img in images: 
        row_sums = img.sum(axis=1)
        weighted_sum = (row_sums * weights).sum()
        leaf_sums.append(weighted_sum)

    linkage_matrix = linkage(squareform(dist_matrix), method="ward")
    reordered_matrix = reorder_linkage_matrix(linkage_matrix, leaf_sums, n_samples)
    pd.DataFrame(reordered_matrix).to_csv(os.path.join(OUTPUT_DIR, "dendrogram_linkage.tsv"), sep="\t", index=False, header=False)

    plt.figure(figsize=(12,5))
    dendro = dendrogram(reordered_matrix, labels=samples, leaf_rotation=90)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "dendrogram.pdf"), bbox_inches="tight")
    plt.close()

    reordered_samples = [samples[i] for i in dendro['leaves']]
    cluster_labels = fcluster(reordered_matrix, t=N_CLUSTERS, criterion='maxclust')
    cluster_labels = [cluster_labels[samples.index(s)] for s in reordered_samples]

    with open(os.path.join(OUTPUT_DIR, "dendrogram_order.tsv"), "w") as f:
        f.write("#sample\tcluster\n")
        for sample, cluster in zip(reordered_samples, cluster_labels):
            f.write(f"{sample}\t{cluster}\n")

    # -------------------------------
    # Heatmap
    # -------------------------------
    heatmap_list = []
    for sample in reordered_samples:
        flatten_image = images[samples.index(sample)].flatten().reshape(-1,1)
        heatmap_list.append(np.log10(flatten_image + 1))
    heatmap_data = np.hstack(heatmap_list)

    n_rows, n_cols = images[0].shape
    fig, ax = plt.subplots(figsize=(len(reordered_samples)*0.3+2, n_cols))

    sns.heatmap(heatmap_data, cmap="mako_r", vmin=0, vmax=2, cbar=True, ax=ax)

    current_label = cluster_labels[0]
    for i, lbl in enumerate(cluster_labels[1:], start=1):
        if lbl != current_label:
            ax.axvline(i, color='black', linewidth=2)
            current_label = lbl

    for h in range(n_cols, n_rows*n_cols, n_cols):
        ax.axhline(h, color='black', linewidth=2)

    ax.set_xticks(np.arange(len(reordered_samples)) + 0.5)
    ax.set_xticklabels(reordered_samples, rotation=90)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "heatmap.pdf"), bbox_inches="tight")
    plt.close()

    heatmap_df = pd.DataFrame(heatmap_data, columns=reordered_samples)
    heatmap_df.to_csv(os.path.join(OUTPUT_DIR, "heatmap_data.tsv"), sep = "\t", index=False) 

if __name__ == "__main__":
    main()