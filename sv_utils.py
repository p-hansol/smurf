import os, re, gzip
import warnings
import numpy as np
import pandas as pd
from collections import Counter


def get_sv(sv_path):
    with open(sv_path) as f:
        lines = [line for line in f if not line.startswith("#")]
    
    sv = pd.DataFrame([line.strip().split("\t")[:6] for line in lines],
                      columns=["chrom1", "pos1", "strand1", "chrom2", "pos2", "strand2"])
    sv[["pos1", "pos2"]] = sv[["pos1", "pos2"]].apply(pd.to_numeric, errors="coerce")

    return sv
    

def add_chr_prefix(chrom):
    chrom = str(chrom)
    chr_chrom = chrom if chrom.startswith("chr") else "chr" + chrom
    return chr_chrom


def get_sv_breakpoint(sv_path):
    """
    Collapse paired breakpoints into a unified breakpoint list. 
    
    Parameters:
        sv_path (str):  
            Path to TSV file. The file should have the first six columns in the following order:
            chrom1, pos1, strand1, chrom2, pos2, strand2.
            Lines starting with '#' are treated as comments and ignored.
    """
    sv_ = get_sv(sv_path)
    sv1 = sv_[["chrom1", "pos1", "strand1"]].copy()
    sv1.columns = ["#chrom", "pos", "strand"]
    sv2 = sv_[["chrom2", "pos2", "strand2"]].copy()
    sv2.columns = ["#chrom", "pos", "strand"]
    
    sv = pd.concat([sv1, sv2], ignore_index=True)
    sv = sv.dropna().drop_duplicates()
    sv["pos"] = sv["pos"].astype("int")
    sv = sv[sv["strand"].isin(["+", "-"])]

    return sv


def get_centromere(cytoband_path):
        if cytoband_path is not None: 
            cytoband = pd.read_csv(gzip.open(cytoband_path, "rt"), sep="\t", header=None)
            centromere = cytoband[cytoband[4] == "acen"]
        else:
            centromere = None
        
        return centromere


def add_centromere_breakpoint(sv, centromere, arm_list):
    chr_prefix = sv["#chrom"].iloc[0].startswith("chr")

    if arm_list:
        arm_set = set([a.strip() for a in arm_list.split(",") if a.strip()])
        for arm in arm_set: 
            match = re.match(r"(\d+|X|Y)([pq])$", arm, re.IGNORECASE)
            if not match:
                raise ValueError(f"Invalid arm format: {arm}")
                
            chrom, pq = match.groups()
            chrom = re.sub(r"^chr", "", chrom)
            chr_chrom = "chr" + chrom

            if pq == "p":
                pos = centromere.loc[centromere[0] == chr_chrom, [1, 2]].values.min()
                strand = "+"
            else: 
                pos = centromere.loc[centromere[0] == chr_chrom, [1, 2]].values.max()
                strand = "-"

            arm_sv = pd.DataFrame([{
                "#chrom": chr_chrom if chr_prefix else chrom,
                "pos": pos,
                "strand": strand
            }])

            sv = pd.concat([sv, arm_sv], ignore_index=True)

    return sv.sort_values(by=["#chrom", "pos", "strand"]).reset_index(drop=True)


def get_chromosome_start(chrom, cen_end=None, acrocentric=["13", "14", "15", "21", "22"]):
    chrom = re.sub(r"^chr", "", str(chrom))
    if chrom in acrocentric and cen_end:
        return cen_end
    return 1


def append_arm_segment(segment, segment_list, chrom, start, end):
    flag = segment["flag"].iloc[0]
    segment = segment.drop(columns=["flag"]).copy()
    arm_segment = pd.DataFrame([{
        "#chrom": chrom, 
        "start": start, 
        "end": end
    }])
    segment = pd.concat([segment, arm_segment], ignore_index=True)
    segment["flag"] = flag + ";cen_adj"
    segment_list.append(segment[["#chrom", "start", "end", "flag"]].sort_values(["start", "end"]).reset_index(drop=True))


def filter_invalid_position(df, chrom, chrom_start, chrom_end):
    df = df[df["#chrom"] == chrom].sort_values(by=["pos", "strand"]).reset_index(drop=True)
    invalid_pos = (df["pos"] < chrom_start) | (df["pos"] > chrom_end)
    if invalid_pos.any():
        removed = df.loc[invalid_pos, "pos"].tolist()
        warnings.warn(
            f"The following positions were removed since they are outside of {chrom}:{chrom_start}-{chrom_end}"
            f": {removed}"
        )
    df = df[~invalid_pos].copy()
    start = set(df[df["strand"] == "-"]["pos"].tolist())
    end = set(df[df["strand"] == "+"]["pos"].tolist())

    return start, end


def segment_genome(start_pos, end_pos, chrom, chrom_start, chrom_end, cen_start, cen_end, flag="from_sv"): 
    """
    Segment the genome by pairing the closest 5′ and 3′ SV breakpoints.
    Unpaired SV breakpoints are extended to the chromosome start and end.
    """
    segment_list = []
    
    start_pos = sorted(start_pos)
    end_pos = sorted(end_pos)
    unpaired_start = start_pos.copy()
    unpaired_end = end_pos.copy()

    # Step1: Pair SV breakpoints
    segment = []
    for end in end_pos:
        candidates = [s for s in unpaired_start if s <= end]
        if candidates:
            start = max(candidates)
            unpaired_start.remove(start)
            unpaired_end.remove(end)
            segment.append({
                "#chrom": chrom, 
                "start": start, 
                "end": end
            })

    # Step 2: Handle unpaired breakpoints
    for end in unpaired_end:
        segment.append({
            "#chrom": chrom, 
            "start": chrom_start, 
            "end": end
        })

    for start in unpaired_start:
        segment.append({
            "#chrom": chrom, 
            "start": start, 
            "end": chrom_end
        })

    segment = pd.DataFrame(segment)
    segment["flag"] = flag
    segment_list.append(segment[["#chrom", "start", "end", "flag"]].sort_values(["start", "end"]).reset_index(drop=True))

    # Step 3: Append an arm-level segment if copy loss is suspected 
    if cen_start and cen_end:
        if segment["start"].min() > cen_end and chrom_start < cen_start: # p-arm
            append_arm_segment(segment, segment_list, chrom, chrom_start, cen_start)
        if segment["end"].max() < cen_start: # q-arm
            append_arm_segment(segment, segment_list, chrom, cen_end, chrom_end)
    
    return segment_list


def assign_copy_number(segment_list, chrom, chrom_start, chrom_end):
    """
    Compute copy number across a chromosome using a sweep line algorithm.
    """
    def compute_copy_number(segment, chrom, chrom_start, chrom_end):
        events = Counter()
        events[chrom_start] += 0 
        events[chrom_end + 1] += 0

        for start, end in zip(segment["start"], segment["end"]):
            events[start] += 1
            events[end + 1] -= 1

        interval_cn = []
        boundary = sorted(events)
        cn = 0

        for i in range(len(boundary) - 1):
            cn += events[boundary[i]]
            start = boundary[i]
            end = boundary[i + 1] - 1
            if start <= end:
                interval_cn.append({
                    "#chrom":chrom, 
                    "start": start, 
                    "end": end, 
                    "cn": cn
                })

        interval_cn = pd.DataFrame(interval_cn)[["#chrom", "start", "end", "cn"]]
        interval_cn["cn"] = interval_cn["cn"].astype(int)
 
        return interval_cn

    segment_cn_list = []
    interval_cn_list = []

    segment_list = segment_list.copy()
    for segment in segment_list:
        segment = segment.copy()
        
        # Compute the copy number of each segment
        interval_cn = compute_copy_number(segment, chrom, chrom_start, chrom_end)
        interval_cn_list.append(interval_cn)  
         
        for i, row in segment.iterrows():
            s, e = row["start"], row["end"]
            matched_interval = interval_cn[
                (interval_cn["start"] == s) & 
                (interval_cn["end"] == e)
            ]
            if not matched_interval.empty:
                segment.at[i, "cn"] = matched_interval["cn"].min()
            else:
                sub_intervals = interval_cn[
                    (interval_cn["start"] >= s) & 
                    (interval_cn["end"] <= e)
                ]
                segment.at[i, "cn"] = sub_intervals["cn"].min() if not sub_intervals.empty else np.nan

        # Add lost segment
        lost_segment = interval_cn[interval_cn["cn"] == 0].copy()
        lost_segment["flag"] = segment["flag"].iloc[0]
        segment_cn = pd.concat([segment, lost_segment], ignore_index=True, sort=True)
        segment_cn["cn"] = segment_cn["cn"].astype(int)
        segment_cn = segment_cn[["#chrom", "start", "end", "cn", "flag"]].sort_values(["start", "end"]).reset_index(drop=True)
        segment_cn_list.append(segment_cn)

    return segment_cn_list, interval_cn_list


def alternative_copy_number(segment_list, interval_list, chrom, chrom_start, chrom_end, insert_size=500):
    """
    Propose an alternative segmentation when a small segment (potential duplication bridge)
    is present and not fully covered by other segments.

    The alternative segmentation includes:
    1. Adding a full-chromosome segment with copy number of 1.
    2. Removing all lost segments (whose copy number is 0).
    3. Incrementing the copy number of all non-lost segments by 1.
    """
    segment_alt_list = []
    interval_alt_list = []

    for segment, interval in zip(segment_list, interval_list):
        # Identify small segments with copy number = 1
        segment = segment.copy()
        segment["size"] = segment["end"] - segment["start"] + 1
        segment_alt_list.append(segment[["#chrom", "start", "end", "size", "cn", "flag"]].copy())
        interval_alt_list.append(interval.copy())

        short_segment = segment[(segment["size"] <= insert_size) & (segment["cn"] == 1)]
        if not short_segment.empty:
            # Create a full-length chromosome segment with copy number = 1
            full_chromosome = pd.DataFrame([{
                "#chrom": chrom,
                "start": chrom_start,
                "end": chrom_end,
                "size": chrom_end,
                "cn": 1,
                "flag": ""
            }])

            # Increase copy number by 1 for all segments with copy number > 0
            non_lost_segment = segment[segment["cn"] > 0].copy()
            non_lost_segment["cn"] += 1

            # Generate an alternative solution
            segment_alt = pd.concat([full_chromosome, non_lost_segment], ignore_index=True, sort=True)
            segment_alt = segment_alt.sort_values(["start", "end"]).reset_index(drop=True)
            segment_alt["flag"] = segment["flag"].iloc[0] + ";alt_sol"
            segment_alt_list.append(segment_alt[["#chrom", "start", "end", "size", "cn", "flag"]])

            # Adjust interval copy number
            interval_alt = interval.copy()
            interval_alt["cn"] = interval_alt["cn"] + 1
            interval_alt_list.append(interval_alt)
    
    return segment_alt_list, interval_alt_list


def find_dup_bridge(segment_list, min_seg_size=20):
    """
    A segment is considered a duplication bridge if its size is smaller than `min_seg_size`.
    """
    segment_dbg_list = [df.copy() for df in segment_list]
    
    for segment_dbg in segment_dbg_list:
        segment_dbg["dup_bridge"] = (
            (segment_dbg["size"] < min_seg_size) & 
            (segment_dbg["cn"] > 1)
        )
    
    return segment_dbg_list


def resolve_dup_bridge(segment_list, interval_list, chrom_start, chrom_end): 
    """
    Identify the shortest segment that contains a duplication bridge, 
    and split it into two overlapping adjacent segments.
    """
    segment_res_list = []
    interval_res_list = []

    for segment, interval in zip(segment_list, interval_list): 
        segment = segment.copy()
        dup_bridges = segment[segment["dup_bridge"]]
        
        for _, row in dup_bridges.iterrows():
            dup_start, dup_end, dup_cn = row["start"], row["end"], row["cn"]
            
            # Find the smallest segment that contains the duplication bridge
            containing = segment[
                (segment["start"] < dup_start) & 
                (segment["end"] > dup_end) & 
                (segment["cn"] > 0) &
                (segment["cn"] < dup_cn) &
                (~segment["dup_bridge"])
            ]
            
            if containing.empty:
                continue
            
            smallest_segment = containing.loc[containing["size"].idxmin()]

            # Split the smallest segment into two overlapping segments
            left_segment = {
                "#chrom": smallest_segment["#chrom"],
                "start": smallest_segment["start"],
                "end": dup_end,
                "size": dup_end - smallest_segment["start"],
                "cn": smallest_segment["cn"],
                "flag": smallest_segment["flag"],
                "dup_bridge": False
            }
            right_segment = {
                "#chrom": smallest_segment["#chrom"],
                "start": dup_start,
                "end": smallest_segment["end"],
                "size": smallest_segment["end"] - dup_start,
                "cn": smallest_segment["cn"],
                "flag": smallest_segment["flag"],
                "dup_bridge": False
            }

            # Drop the original segment and add the split segments
            segment = segment.drop(index=smallest_segment.name)
            segment = pd.concat([segment, pd.DataFrame([left_segment, right_segment])], ignore_index=True, sort=True)
            
        # Return None if full-chromosome segment remains after resolving duplication bridges.
        full_chromosome = segment[(segment["start"] == chrom_start) & (segment["end"] == chrom_end)]
        if full_chromosome.empty:
            segment = segment[["#chrom", "start", "end", "size", "cn", "flag", "dup_bridge"]].sort_values(["start", "end"]).reset_index(drop=True)
            segment_res_list.append(segment)
            interval_res_list.append(interval)
 
    return segment_res_list.copy(), interval_res_list.copy()


def rank_solution(segment_list, interval_list, chrom_start, chrom_end):
    """
    Rank segment and interval solutions using 'flag'.

    Ranking Criteria:
    1. Penalize segments with excessive copy number gain.
    2. Prefer segments with alternative solutions.
    3. Prefer SV-based solutions over CNV-adjusted ones.
    4. Penalize segments adjusted near the centromere.
    """

    def has_excessive_gain(segment, chrom_start, chrom_end):
        """
        Determine if a segment exhibits excessive gain.

        Criteria:
        1. Total size of two-copy segments exceeds 90% of the chromosome size.
        2. There are two or more segments starting at the chromosome start.
        3. There are two or more segments ending at the chromosome end.
        """
        two_copy_total_size = segment.loc[segment["cn"] == 2, "size"].sum()
        starts_at_chrom_start = (segment["start"] == chrom_start).sum()
        ends_at_chrom_end = (segment["end"] == chrom_end).sum()

        return (
            two_copy_total_size > chrom_end * 0.9 or
            starts_at_chrom_start >= 2 or
            ends_at_chrom_end >= 2
        )

    def flag_to_binary_score(flag):
        """
        Convert a semicolon-separated flag string to a 4-bit binary score.
        """
        flag_set = set(flag.split(";"))

        b1 = int("ex_gain" in flag_set)
        b2 = int("alt_sol" not in flag_set)
        b3 = int("from_sv" not in flag_set)
        b4 = int("cen_adj" in flag_set)

        return int(f"{b1}{b2}{b3}{b4}", 2)
    
    scores = []
    for i, segment in enumerate(segment_list):
        if has_excessive_gain(segment, chrom_start, chrom_end):
            segment["flag"] = segment["flag"].iloc[0] + ";ex_gain"

        flag = segment["flag"].iloc[0]
        score = flag_to_binary_score(flag)
        scores.append((i, score))

    sorted_indices = [i for i, _ in sorted(scores, key=lambda x: x[1])]
    segment_rank_list = [segment_list[i] for i in sorted_indices]
    interval_rank_list = [interval_list[i] for i in sorted_indices]

    return segment_rank_list, interval_rank_list


def save_solution(segment_list, chr_chrom, outdir, sample):
    os.makedirs(outdir, exist_ok=True)
    for rank, segment in enumerate(segment_list):
        segment = segment[["#chrom", "start", "end", "cn"]].copy()
        segment = segment.rename(columns={"cn": "copy_number"})
        filename = f"{outdir}/{sample}.{chr_chrom}.solution{rank + 1}.tsv"
        segment.to_csv(filename, sep="\t", index=False)


def save_decision(chrom_list, outdir, sample):
    def order_chrom(chrom):
        chrom = chrom.lower().replace("chr", "")
        if chrom.isdigit():
            return (int(chrom), '')
        order = {'x': 23, 'y': 24, 'm': 25, 'mt': 25}
        return (order.get(chrom, 100), chrom)

    os.makedirs(outdir, exist_ok=True)
    sorted_chrom = sorted(chrom_list, key=order_chrom)
    decision = pd.DataFrame({
        "#chrom": sorted_chrom,
        "solution": [1] * len(sorted_chrom)
    })
    filename = f"{outdir}/{sample}.decision.tsv"
    decision.to_csv(filename, sep="\t", index=False)
