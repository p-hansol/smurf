import argparse, os, sys
import pandas as pd
import sv_utils
import plot_utils


def main():
    parser = argparse.ArgumentParser(description="Segment genome and assign copy number based on structural variations.")

    parser.add_argument("-i", "--input", required=True,
                        help="Path to SV TSV with columns: #chrom1, pos1, strand1, chrom2, pos2, strand2.")
    parser.add_argument("-s", "--sample", required=True,
                        help="Sample ID (used as output file prefix).")
    parser.add_argument("-r", "--ref_fasta", required=True,
                        help="Path to reference FASTA. Must have corresponding .fai index.")
    parser.add_argument("-v", "--ref_ver", default=None, choices=["hg19", "hg38"],
                        help="Reference genome version used to retrieve centromere positions (either 'hg19' or 'hg38'). "
                             "If not provided, centromere information will not be used.")
    parser.add_argument("-g", "--gene", default=None, 
                        help="Comma-separated gene symbols (e.g., ALK,EML4). Used to highlight gene positions in the output plot.") 
    parser.add_argument("-o", "--outdir", default=os.getcwd(),
                        help="Output directory (default: current directory).")
    parser.add_argument("--search_missing_breakpoint", action="store_true",
                        help="Use CNV to supplement undetected SV breakpoints.")
    parser.add_argument("-c", "--cnv", default=None,
                        help="Path to CNV TSV with required columns: #chrom, start, end, total_cn; and optional columns: major_cn, minor_cn.")
    parser.add_argument("--cnv_identical_dist", type=int, default=1000,
                        help="CNV boundaries within this distance from an SV breakpoint are excluded from segmentation (i.e., treated as identical to the SV breakpoint).")
    parser.add_argument("--cnv_association_dist", type=int, default=5000000,
                        help="CNV boundaries within this distance from an SV breakpoint are included in segmentation (i.e., regarded as associated with the SV breakpoint).")
    parser.add_argument("-a", "--arm", default=None, 
                        help="Comma-separated chromosome arms (e.g., 2p,5q). Arm-level segments corresponding to the specified chromosome arms will be added.")
    parser.add_argument("--search_duplication_bridge", action="store_true",
                        help="Examine the BAM/CRAM file to determine whether a short copy gain arises from either inclusive or intersecting segments.")
    parser.add_argument("-b", "--bam", default=None,
                        help="Path to BAM or CRAM file.")
    parser.add_argument("--min_seg_size", type=int, default=20,
                        help="Segments smaller than this are automatically considered duplication bridges.")
    parser.add_argument("--insert_size", type=int, default=500,
                        help="Insert size of paired-end reads. Segments smaller than this are checked in BAM/CRAM for duplication bridges.")

    args = parser.parse_args()


    if args.gene and args.ref_ver is None: 
        sys.exit("Error: --ref_ver is required when --gene is provided.")
    if args.arm and args.ref_ver is None: 
        sys.exit("Error: --ref_ver is required when --arm is provided.")
    if args.search_missing_breakpoint and args.cnv is None:
        sys.exit("Error: --cnv is required when --search_missing_breakpoint is enabled.")
    if args.search_duplication_bridge and args.bam is None:
        sys.exit("Error: --bam is required when --search_duplication_bridge is enabled.")


    os.makedirs(args.outdir, exist_ok=True)


    if args.search_missing_breakpoint:
        import cnv_utils
        cnv = cnv_utils.get_cnv_boundary(args.cnv)
        segment_genome = lambda start_pos, end_pos, chrom, chrom_start, chrom_end, cen_start, cen_end: cnv_utils.segment_genome_sv_cnv(
            start_pos, end_pos, chrom, chrom_start, chrom_end, cen_start, cen_end, cnv, args.cnv_identical_dist, args.cnv_association_dist)
    else:
        segment_genome = lambda start_pos, end_pos, chrom, chrom_start, chrom_end, cen_start, cen_end: sv_utils.segment_genome(
            start_pos, end_pos, chrom, chrom_start, chrom_end, cen_start, cen_end)

    if args.search_duplication_bridge:
        import bam_utils
        find_dup_bridge = lambda segment_cn: bam_utils.find_dup_bridge_bam(
            segment_cn, args.bam, args.input, args.ref_fasta, args.insert_size, args.min_seg_size)
    else:
        find_dup_bridge = lambda segment_cn: sv_utils.find_dup_bridge(
            segment_cn, args.min_seg_size)
    

    ref_index = pd.read_csv(args.ref_fasta + ".fai", sep="\t", header=None, usecols=[0, 1], names=["#chrom", "size"], dtype={"size": int})
    cytoband = os.path.join(os.path.dirname(__file__), "resource", args.ref_ver, "cytoBand.txt.gz") if args.ref_ver else None
    centromere = sv_utils.get_centromere(cytoband)
    gene_gff3 = os.path.join(os.path.dirname(__file__), "resource", args.ref_ver, "gencode.v48.basic.annotation.gene.gff3.gz") if args.ref_ver else None
    gene = plot_utils.get_gene_position(args.gene, gene_gff3)
    sv = sv_utils.get_sv_breakpoint(args.input)
    sv = sv_utils.add_centromere_breakpoint(sv, centromere, args.arm)


    for chrom in sv["#chrom"].unique():
        chr_chrom = sv_utils.add_chr_prefix(chrom)

        # Get centomere position
        if centromere is not None:
            cen_start = centromere.loc[centromere[0] == chr_chrom, [1, 2]].values.min()
            cen_end = centromere.loc[centromere[0] == chr_chrom, [1, 2]].values.max()
        else:
            cen_start = cen_end = None

        # Get chromosome size
        try:
            chrom_start = sv_utils.get_chromosome_start(chrom, cen_end)
            chrom_end = ref_index[ref_index["#chrom"] == chrom]["size"].values[0]
        except IndexError:
            sys.exit(f"Error: {chr_chrom} not found in FASTA index.")

        # Segment the genome
        start_pos, end_pos = sv_utils.filter_invalid_position(sv, chrom, chrom_start, chrom_end)
        segment = segment_genome(start_pos, end_pos, chrom, chrom_start, chrom_end, cen_start, cen_end)
        
        # Compute copy number of each segment
        segment_cn, interval_cn = sv_utils.assign_copy_number(segment, chrom, chrom_start, chrom_end)
        segment_alt, interval_alt = sv_utils.alternative_copy_number(segment_cn, interval_cn, chrom, chrom_start, chrom_end, args.insert_size)
        
        # Resolve duplication bridges
        segment_dbg = find_dup_bridge(segment_alt)
        segment_res, interval_res = sv_utils.resolve_dup_bridge(segment_dbg, interval_alt, chrom_start, chrom_end)
        
        # Rank the solutions
        segment_rank, interval_rank = sv_utils.rank_solution(segment_res, interval_res, chrom_start, chrom_end)
        sv_utils.save_solution(segment_rank, chr_chrom, args.outdir, args.sample)
        
        # Plot the solutions
        plot_utils.plot_solution(segment_rank, interval_rank, chr_chrom, gene, cen_start, cen_end, args.outdir, args.sample)

    # Save final decision
    sv_utils.save_decision(sv["#chrom"].unique(), args.outdir, args.sample)


if __name__ == "__main__":
    main()
