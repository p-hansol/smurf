import pandas as pd
import sv_utils


def get_cnv_boundary(cnv_path=None):
    """
    Get segment boundaries (start and end positions) from segmental copy number profiles.

    Parameters:
        cnv_path (str or None): 
        Path to a TSV file with required columns: #chrom, start, end, total_cn; 
        and optional columns: major_cn, minor_cn 
        Lines starting with '#' are treated as comments and ignored.
    """
    if cnv_path is None:
        return None

    with open(cnv_path) as f:
        lines = [line.strip().split("\t") for line in f if not line.startswith("#")]

    if len(lines[0]) >= 6:
        cols = ["#chrom", "start", "end", "tot_cn", "maj_cn", "min_cn"]
    elif len(lines[0]) >= 4:
        cols = ["#chrom", "start", "end", "tot_cn"]
    else:
        raise ValueError("CNV file must have at least 4 columns: #chrom, start, end, total_cn.")
    
    cnv_full = pd.DataFrame(lines, columns=cols)
    cnv_full[cols[1:]] = cnv_full[cols[1:]].apply(pd.to_numeric, errors="coerce")

    cnv = []
    for chrom in cnv_full["#chrom"].unique():
        cnv_chrom = cnv_full[cnv_full["#chrom"] == chrom].reset_index(drop=True)

        for cn_col in cols[3:]:
            cnv_cn = cnv_chrom[cols[:3] + [cn_col]] \
                     .dropna(subset=["start", "end", cn_col]) \
                     .sort_values(by=["start"]) \
                     .reset_index(drop=True)

            for i in range(1, len(cnv_cn)):
                cn1 = cnv_cn.loc[i - 1, cn_col]
                cn2 = cnv_cn.loc[i, cn_col]
                
                if cn1 < cn2:
                    cnv.append({
                        "#chrom": chrom, 
                        "pos": cnv_cn.loc[i, "start"], 
                        "strand": "-"
                    })
                elif cn1 > cn2:
                    cnv.append({
                        "#chrom": chrom, 
                        "pos": cnv_cn.loc[i - 1, "end"], 
                        "strand": "+"
                    })
    
    cnv = pd.DataFrame(cnv, columns=["#chrom", "pos", "strand"])
    cnv = cnv.drop_duplicates().reset_index(drop=True)

    return cnv


def filter_cnv_boundary(cnv, sv_start, sv_end, chrom, chrom_start, chrom_end, cen_start, cen_end, identical_dist, association_dist):
    """
    Include CNV boundaries likely associated with SV breakpoints (within association_dist) or near the centromere.
    Exclude CNV boundaries likely representing the same event as SV breakpoints (closest distance within identical_dist).

    Parameters:
        identical_dist: Maximum distance to consider a CNV boundary as identical to an SV breakpoint.
        association_dist: Maximum distance to consider a CNV boundary as associated with an SV breakpoint.
    """
    def exclude_cnv_boundary(cnv_pos_list, sv_pos_list, identical_dist):
        cnv_pos_list = cnv_pos_list.copy()

        for sv_pos in sv_pos_list:
            if cnv_pos_list:
                cnv_pos = min(cnv_pos_list, key = lambda p: abs(p - sv_pos))
                if abs(cnv_pos - sv_pos) <= identical_dist:
                     cnv_pos_list.remove(cnv_pos)
        
        return cnv_pos_list

    cnv_chrom = cnv[cnv["#chrom"] == chrom]
    cnv_start = []
    cnv_end = []
    
    if not cnv_chrom.empty:
        sv_breakpoint = sv_start + sv_end
        if cen_start and cen_end: 
            sv_breakpoint += [cen_start, cen_end]
        
        cnv_chrom = cnv_chrom[cnv_chrom["pos"].apply(
            lambda pos: any(abs(pos - sv_pos) <= association_dist for sv_pos in sv_breakpoint)
        )]
        cnv_start, cnv_end = sv_utils.filter_invalid_position(cnv_chrom, chrom, chrom_start, chrom_end)
        cnv_start = exclude_cnv_boundary(cnv_start, sv_start, identical_dist)
        cnv_end = exclude_cnv_boundary(cnv_end, sv_end, identical_dist)

    return cnv_start, cnv_end


def segment_genome_cnv(start_pos, end_pos, chrom, chrom_start, chrom_end, cen_start, cen_end, cnv, flag="cnv_adj", identical_dist=1000, association_dist=5000000):
    """
    Segment the genome by pairing the closest 5′ and 3′ SV breakpoints.
    Unpaired SV breakpoints are extended to the chromosome start and end or CNV boundaries.
    """
    
    sv_start = sorted(start_pos)
    sv_end = sorted(end_pos)
    unpaired_sv_start = sv_start.copy()
    unpaired_sv_end = sv_end.copy()
    
    sv_cnv_start = sv_start.copy()
    sv_cnv_end = sv_end.copy()
    
    # Filter CNV boundaries
    unpaired_cnv_start, unpaired_cnv_end = filter_cnv_boundary(
        cnv, sv_start, sv_end, chrom, chrom_start, chrom_end, cen_start, cen_end, identical_dist, association_dist
    )
    
    # Supplement CNV boundaries 
    for end in sv_end:
        candidates = [s for s in unpaired_sv_start if s < end]
        if candidates:
            start = max(candidates)
            unpaired_sv_start.remove(start)
            unpaired_sv_end.remove(end)
    
    for end in unpaired_sv_end:
        if end == min(unpaired_sv_end):
            continue
    
        # Find the end of upstream segment
        upstream_ends = [e for e in sv_end if e < end]
        upstream_end = max(upstream_ends) if upstream_ends else 0

        in_btw_cnv_start = [s for s in unpaired_cnv_start if upstream_end < s < end]
        in_btw_cnv_end = [e for e in unpaired_cnv_end if upstream_end < e < end]
        upstream_cnv_start = [s for s in unpaired_cnv_start if s < upstream_end]

        if in_btw_cnv_start:
            # Prefer to minimize the gap with the upstream segment
            start = min(in_btw_cnv_start)
            unpaired_cnv_start.remove(start)
        elif in_btw_cnv_end and upstream_cnv_start:
            # Prefer to create the smallest possible segment
            start = max(upstream_cnv_start)
            unpaired_cnv_start.remove(start)
        else: 
            # If no copy number change, start just after the upstream segment
            start = upstream_end + 1 

        sv_cnv_start.append(start)

    for start in unpaired_sv_start:
        if start == max(unpaired_sv_start):
            continue
    
        # Find the start of downstream segment
        downstream_starts = [s for s in sv_start if s > start]
        downstream_start = min(downstream_starts) if downstream_starts else chrom_end + 1

        in_btw_cnv_end = [e for e in unpaired_cnv_end if start < e < downstream_start]
        in_btw_cnv_start = [s for s in unpaired_cnv_start if start < s < downstream_start]
        downstream_cnv_end = [e for e in unpaired_cnv_end if e > downstream_start]

        if in_btw_cnv_end:
            # Prefer to minimize the gap with the downstream segment
            end = max(in_btw_cnv_end)
            unpaired_cnv_end.remove(end)
        elif in_btw_cnv_start and downstream_cnv_end: 
            # Prefer to create the smallest possible segment
            end = min(downstream_cnv_end)
            unpaired_cnv_end.remove(end)
        else: 
            # If no copy number change, end just before the downstream segment
            end = downstream_start - 1
        
        sv_cnv_end.append(end)
    
    # Segment the genome using SV and CNV
    segment_list = sv_utils.segment_genome(sv_cnv_start, sv_cnv_end, chrom, chrom_start, chrom_end, cen_start, cen_end, flag)
     
    return segment_list


def segment_genome_sv_cnv(start_pos, end_pos, chrom, chrom_start, chrom_end, cen_start, cen_end, cnv, identical_dist=1000, association_dist=5000000):
    segment_list = sv_utils.segment_genome(
        start_pos, end_pos, chrom, chrom_start, chrom_end, cen_start, cen_end, flag="from_sv"
    )
    segment_list_add = segment_genome_cnv(
        start_pos, end_pos, chrom, chrom_start, chrom_end, cen_start, cen_end, cnv, flag="cnv_adj", 
        identical_dist=identical_dist, association_dist=association_dist
    )

    for segment_add in segment_list_add:
        is_duplicate = any(
            segment_add.drop(columns=["flag"]).equals(segment.drop(columns=["flag"]))
            for segment in segment_list
        )
        if not is_duplicate:
            segment_list.append(segment_add)

    return segment_list