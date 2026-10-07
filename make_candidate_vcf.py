#!/usr/bin/env python3
"""
Convert a PyClone-VI results TSV into a minimal candidate-SNP VCF for cellsnp-lite -R.
Keeps cluster_id as an INFO field so we can map cells back to bulk clones later.
Usage: python3 make_candidate_vcf.py <results_tsv> <output_vcf>
"""
import sys, pandas as pd

results_tsv, out_vcf = sys.argv[1], sys.argv[2]
df = pd.read_csv(results_tsv, sep="\t")

# one row per mutation per sample_id in pyclone output; dedupe since sample_id is constant here
df = df.drop_duplicates(subset="mutation_id")

with open(out_vcf, "w") as f:
    f.write("##fileformat=VCFv4.2\n")
    f.write('##INFO=<ID=CLUSTER,Number=1,Type=Integer,Description="PyClone-VI bulk cluster ID">\n')
    f.write('##INFO=<ID=CCF,Number=1,Type=Float,Description="Bulk cellular prevalence">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for _, row in df.iterrows():
        chrom, pos, ref, alt = row["mutation_id"].split("_", 3)
        f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\tCLUSTER={row['cluster_id']};CCF={row['cellular_prevalence']:.4f}\n")

print(f"Wrote {len(df)} candidate sites to {out_vcf}")
