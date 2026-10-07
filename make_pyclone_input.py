#!/usr/bin/env python3
"""
Build a PyClone-VI input TSV for one patient from a PURPLE somatic VCF.
v2: adds indel filter, germline-CN-state filter, and mappability filter.
Usage: make_pyclone_input.py <patient_id> <vcf_path> <purity_tsv> <tumor_sample_name> <output_tsv>
"""
import sys, subprocess, csv

patient_id, vcf_path, purity_tsv, tumor_sample, out_path = sys.argv[1:6]

with open(purity_tsv) as f:
    row = next(csv.DictReader(f, delimiter='\t'))
    purity = float(row['purity'])
    gender = row['gender']

cmd = [
    "bcftools", "query",
    "-s", tumor_sample,
    "-f", "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/PURPLE_CN\t%INFO/PURPLE_MACN\t%INFO/TIER\t%INFO/PURPLE_GERMLINE\t%INFO/MAPPABILITY\t[%AD]\n",
    vcf_path,
]
result = subprocess.run(cmd, capture_output=True, text=True, check=True)

KEEP_GERMLINE_STATES = {"DIPLOID", "HET"}
MIN_MAPPABILITY = 0.8

rows_out = []
n_seen = 0
n_dropped_filter = 0
n_dropped_tier = 0
n_dropped_cn_missing = 0
n_dropped_indel = 0
n_dropped_germline = 0
n_dropped_mappability = 0
n_dropped_depth = 0

for line in result.stdout.strip().split("\n"):
    if not line:
        continue
    n_seen += 1
    chrom, pos, ref, alt, filt, cn, macn, tier, germline_state, mappability, ad = line.split("\t")

    if filt != "PASS":
        n_dropped_filter += 1
        continue
    if tier == "LOW_CONFIDENCE":
        n_dropped_tier += 1
        continue
    if cn in (".", "") or macn in (".", ""):
        n_dropped_cn_missing += 1
        continue
    if len(ref) != 1 or len(alt) != 1:
        n_dropped_indel += 1
        continue
    if germline_state not in KEEP_GERMLINE_STATES:
        n_dropped_germline += 1
        continue
    try:
        if float(mappability) < MIN_MAPPABILITY:
            n_dropped_mappability += 1
            continue
    except ValueError:
        n_dropped_mappability += 1
        continue

    ref_count, alt_count = (int(x) for x in ad.split(","))
    depth = ref_count + alt_count
    if depth < 20 or alt_count < 3:
        n_dropped_depth += 1
        continue

    cn_f, macn_f = float(cn), float(macn)
    minor_cn = max(0, round(macn_f))
    major_cn = max(1, round(cn_f - macn_f))

    normal_cn = 2
    if chrom.lstrip("chr") in ("X", "Y") and gender == "MALE":
        normal_cn = 1

    mutation_id = f"{chrom}_{pos}_{ref}_{alt}"
    rows_out.append([mutation_id, patient_id, ref_count, alt_count,
                      major_cn, minor_cn, normal_cn, f"{purity:.4f}", "0.001"])

with open(out_path, "w") as f:
    f.write("mutation_id\tsample_id\tref_counts\talt_counts\tmajor_cn\tminor_cn\tnormal_cn\ttumour_content\terror_rate\n")
    for r in rows_out:
        f.write("\t".join(str(x) for x in r) + "\n")

print(f"{patient_id}: {n_seen} seen | kept {len(rows_out)} | dropped: "
      f"filter={n_dropped_filter} tier={n_dropped_tier} cn_missing={n_dropped_cn_missing} "
      f"indel={n_dropped_indel} germline={n_dropped_germline} mappability={n_dropped_mappability} depth={n_dropped_depth} "
      f"-> {out_path}")
