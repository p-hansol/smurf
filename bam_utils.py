import pysam
from collections import Counter
import sv_utils


def find_mate_position(sv, chrom, pos):
    match1 = sv[(sv["chrom1"] == chrom) & (sv["pos1"] == pos)].dropna()
    if not match1.empty:
        return tuple(match1[["chrom2", "pos2"]].iloc[0])
    
    match2 = sv[(sv["chrom2"] == chrom) & (sv["pos2"] == pos)].dropna()
    if not match2.empty:
        return tuple(match2[["chrom1", "pos1"]].iloc[0])

    return (None, None) 


def is_dup_bridge(bam, sv, chrom, start, end, insert_size, min_seg_size):
    """
    Evaluate whether the observed short copy gain is due to a duplication bridge 
    involving the overlap of two adjacent segments.
    """
    # Step0: If the copy gain is very short, consider it a duplication bridge.
    dup_bridge = end - start < min_seg_size
    if dup_bridge:
        return True

    # Step1: Get variant reads
    variant_reads1 = []
    variant_reads2 = []

    mate_chrom1, mate_pos1 = find_mate_position(sv, chrom, start)
    mate_chrom2, mate_pos2 = find_mate_position(sv, chrom, end)

    for read in bam.fetch(chrom, start - insert_size, end + insert_size):
        if not read.is_paired or read.is_unmapped or read.mate_is_unmapped:
            continue
        read_mate_chrom = bam.get_reference_name(read.next_reference_id)
        read_mate_pos = read.next_reference_start

        # =: clipped base, -: matched base

        #     [Segement start
        # <===|--- 
        #    <|------    Reverse variant reads
        #     | <------
        if (
            read.is_reverse and 
            read_mate_chrom == mate_chrom1 and 
            abs(read_mate_pos - mate_pos1) < insert_size
        ): 
            variant_reads1.append(read)
        
        #                  Segment end]
        #                          ---|===>
        # Forward variant reads ------|>
        #                    -------> |
        elif (
            not read.is_reverse and 
            read_mate_chrom == mate_chrom2 and 
            abs(read_mate_pos - mate_pos2) < insert_size
        ): 
            variant_reads2.append(read)

    # Step2: Find the exact breakpoint from clipped reads

    #     [Segement start
    # <===|---  Left-clipped read
    left_clipped_pos = [
        read.reference_start
        for read in variant_reads1
        if read.cigartuples and read.cigartuples[0][0] in (4, 5)
    ]
    left_clipped_pos = sorted(left_clipped_pos, key=lambda x: abs(x - start))
    exact_start = Counter(left_clipped_pos).most_common(1)[0][0] if left_clipped_pos else None

    #               Segment end]
    # Right-clipped read    ---|===>
    right_clipped_pos = [
        read.reference_end
        for read in variant_reads2
        if read.cigartuples and read.cigartuples[-1][0] in (4, 5)
    ]
    right_clipped_pos = sorted(right_clipped_pos, key=lambda x: abs(x - end))
    exact_end = Counter(right_clipped_pos).most_common(1)[0][0] if right_clipped_pos else None

    # Step 3: Check if any variant read spans outside the exact breakpoint

    # Duplication segment                          vs.  Duplication brige
    #      [=====]                                  .  [==========]                                 .
    # [===============]                             .       [==========]                            .
    #                                               .                                               .
    #   <==|--                                      .    <==|--                                     .
    #         <--|==    # ends at the exact end     .          <-----    # beyond the exact end     .
    #          --|==>                               .           --|==>                              .
    #    ==|-->         # starts at the exact start .    ----->         # beyond the exact start    .
    
    if exact_end is not None: 
        for read in variant_reads1:
            if read.reference_end > exact_end + 1:
                dup_bridge=True
                break

    if exact_start is not None:
        for read in variant_reads2:
            if read.reference_start < exact_start - 1:
                dup_bridge=True
                break
    
    return dup_bridge


def find_dup_bridge_bam(segment_cn_alt_list, bam_path, sv_path, ref_fasta_path, insert_size=500, min_seg_size=20): 
    """
    Examine the BAM/CRAM file to determine whether a short copy gain is a duplication bridge. 
    """
    bam = pysam.AlignmentFile(bam_path, "rb", reference_filename=ref_fasta_path)
    sv = sv_utils.get_sv(sv_path)

    segment_cn_alt_dbg_list = [df.copy() for df in segment_cn_alt_list]

    for segment_cn_alt_dbg in segment_cn_alt_dbg_list:
        segment_cn_alt_dbg["dup_bridge"] = False  
    
        for i, row in segment_cn_alt_dbg.iterrows():
            if row["size"] <= insert_size and row["cn"] > 1: 
                chrom, start, end = row["#chrom"], row["start"], row["end"]
                if is_dup_bridge(bam, sv, chrom, start, end, insert_size, min_seg_size): 
                    segment_cn_alt_dbg.at[i, "dup_bridge"] = True
    
    bam.close()

    return segment_cn_alt_dbg_list