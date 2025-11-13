# SMURF: SSIM-Mediated Unsupervised Rearrangement clustering Framework

**SMURF** quantifies similarities among genomic rearrangements—considering rearrangement type, segment length, and copy number changes—using the Multi-Scale Structural Similarity Index Measure (MS-SSIM), and clusters them accordingly.
Please note that **SMURF** is designed to assess the similarity between rearrangements—either simple or complex—that are derived from a *single* event, such as a simple translocation, chromothripsis, or chromoplexy.

## Application
[figure]
This example demonstrates how SMURF clusters *ALK rearrangements* based on their genomic features.  
  
**SMURF workflow:**  
i) Take structural variations as input.  
ii) Perform segmentation and copy number estimation.  
iii) Generate a simple image for each case, where segments are categorized by size and copy number.  
iv) Quantify similarity between images using windows of different sizes to capture both local and global patterns.  

## Prerequisites
* Python 3.6+
* Required external packages: `numpy`, `pandas`, `scipy`, `matplotlib`, `seaborn`, `pysam`.

## Installation
```bash
git clone https://github.com/p-hansol/smurf.git
cd smurf
pip install -r requirements.txt
```

## Usage 
### 1. Segmentation and Copy Number Estimation
#### * Basic Usage
It is possible for multiple independent events that cause copy number alterations to occur at the same genomic locus. **SMURF** requires only structural variations (SVs) and performs segmentation and copy number estimation based *solely* on these SVs. Therefore, it captures copy number changes that are derived from SVs, but not from other independent events.
```bash
python main.py -i input_sv.tsv -s sample_id -r reference.fasta -v hg38
```
#### * Advanced Usages
If there are missing breakpoints (e.g., due to breakpoints located in repetitive regions), **SMURF** may produce unexpectedly large segments with elevated copy numbers. In such cases, providing additional copy number information is **highly recommended**, as it can help mitigate this issue.
```bash
python main.py -i input_sv.tsv -s sample_id -r reference.fasta -v hg38 \
  -c input_cnv.tsv --search_missing_breakpoint
```
Short copy gains can arise from either *inclusive* or *intersecting* segments. If a BAM/CRAM file is provided, **SMURF** can distinguish between these two cases by examining variant-supporting reads. This analysis is limited to segments smaller than the insert size.
```bash
python main.py -i input_sv.tsv -s sample_id -r reference.fasta -v hg38 \
  -b input.bam --search_duplication_bridge
```
#### * Arguments
| Argument | Type / Default | Description |
|----------|----------------|--------------|
| `-i`, `--input` | **Required** | Path to SV TSV file with columns: `#chrom1`, `pos1`, `strand1`, `chrom2`, `pos2`, `strand2`. |
| `-s`, `--sample` | **Required** | Sample ID (used as output file prefix). |
| `-r`, `--ref_fasta` | **Required** | Path to reference FASTA file. Must have a corresponding `.fai` index. |
| `-v`, `--ref_ver` | `None` <br> *(choices: `hg19`, `hg38`)* | Reference genome version used to retrieve centromere positions. If not provided, centromere information will not be used. |
| `-g`, `--gene` | `None` | Comma-separated gene symbols (e.g., `ALK,EML4`). Used to highlight gene positions in the output plot. `--ref_ver` is required when `--gene` is provided. |
| `-o`, `--outdir` | Current directory | Output directory (default: current directory). |
| `--search_missing_breakpoint` | Flag | Use CNV data to supplement undetected SV breakpoints. `--cnv` is required when `--search_missing_breakpoint` is enabled. |
| `-c`, `--cnv` | `None` | Path to CNV TSV file with required columns: `#chrom`, `start`, `end`, `total_cn`; and optional columns: `major_cn`, `minor_cn`. |
| `--cnv_identical_dist` | `1000` | CNV boundaries within this distance from an SV breakpoint are excluded from segmentation (treated as identical to the breakpoint). |
| `--cnv_association_dist` | `5000000` | CNV boundaries within this distance from an SV breakpoint are included in segmentation (regarded as associated with the breakpoint). |
| `-a`, `--arm` | `None` | Comma-separated chromosome arms (e.g., `2p,5q`). Arm-level segments for the specified chromosome arms will be added. `--ref_ver` is required when `--arm` is provided. |
| `--search_duplication_bridge` | Flag | Examine the BAM/CRAM file to determine whether a short copy gain arises from inclusive or intersecting segments. `--bam` is required when `--search_duplication_bridge` is enabled. |
| `-b`, `--bam` | `None` | Path to BAM or CRAM file. |
| `--min_seg_size` | `20` | Segments smaller than this are automatically considered duplication bridges. |
| `--insert_size` | `500` | Insert size of paired-end reads. Segments smaller than this are checked in BAM/CRAM for duplication bridges. |
#### * Output
- `solution.tsv`: TSV file containing one or more candidate solutions for each chromosome, with columns `#chrom`, `start`, `end`, and `copy_number`.  
- `solution.pdf`: Visual representation of all candidate solutions.  
- `decision.tsv`: TSV file recording the most plausible solution for each chromosome, with columns `#chrom` and `solution`.
  
### 2. Simplification
#### * Usage
```bash
python summary.py -i input_directory/ -c chrN
```
#### * Arguments
| Argument | Type / Default | Description |
|----------|----------------|--------------|
| `-i`, `--input` | **Required** | Path to `decision.tsv` file. |
| `-c`, `--chrA` | **Required** | Chromosome of interest (e.g., `chr2`). |
| `-o`, `--outdir` | Current directory | Output directory (default: current directory). |
#### * Output
- `summary.tsv`: TSV file where rows represent copy numbers and columns represent segment sizes.  
- `summary.pdf`: Visual representation of the data in `summary.tsv`.
  
### 3. Unsupervised Clustering 
#### * Usage
```bash
python cluster.py -i input_direcotry/ 
```
#### * Arguments
| Argument | Type / Default | Description |
|----------|----------------|-------------|
| `-i`, `--indir` | **Required** | Input directory containing `summary.tsv` files. |
| `-o`, `--outdir` | Current directory | Output directory (default: current directory). |
| `-w`, `--weight` | `0.4` | Weight for the within-sample MS-SSIM component (float between 0 and 1). The across-sample component weight is automatically `1 - weight`. |
| `-n`, `--n_cluster` | `5` | Number of clusters. |
#### * Output
- `dendrogram.pdf`: Hierarchical clustering dendrogram of the SVs based on MS-SSIM  
- `heatmap.pdf`: Heatmap of the SVs. Columns represent individual SVs, and rows corresponds to segments defined by size and copy number. Cell colors indicate the number of corresponding segments.  
