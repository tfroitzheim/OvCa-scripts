# OvCa-scripts
Sharing scripts for natgen revisions with Matias 

--> Validate CNA-based single-cell clone calls (ASCENT) against bulk WGS-derived subclonal structure (PyClone-VI, DPClust) by mapping SNPs back and forth between bulk and single-cell data.



Core pipeline:

1. Bulk WGS input (HMF/PURPLE pipeline)
Annoying thing for us: Data was hg19/GRCh37, but the single-cell ASCENT BAMs are hg38 — requires CrossMap liftover.
Maybe similar for you? Worth double checking before doing analyses. We used CrossMap liftover.
Filter somatic VCF: FILTER=PASS, drop TIER=LOW_CONFIDENCE, drop indels (cellsnp-lite is SNV-only), require PURPLE_GERMLINE state in {DIPLOID, HET}, MAPPABILITY ≥0.9, depth ≥20, alt≥3


2. Clonal deconvolution on bulk mutations
PyClone-VI: build input (ref/alt counts, major/minor CN rounded from PURPLE_CN/MACN, purity), fit with beta-binomial, sweep --mix-weight-prior to check cluster stability
DPClust: cross-check tool; critical fix — must properly estimate no.chrs.bearing.mut (multiplicity) rather than assuming 1, or spurious high-CN clusters appear


3. Single-cell side (cellsnp-lite)
Build per-patient candidate SNP whitelist from PyClone-VI results (CLUSTER/CCF in VCF INFO)
Run cellsnp-lite in targeted mode (-R) against ASCENT's single-cell DNA+RNA BAMs, --minCOUNT 1 --minMAF 0 (shallow coverage means standard thresholds filter out everything)
Per-patient whitelist only — not one whitelist for all patients. 


4. Joining back to ASCENT clones
ASCENT clone labels come from .phased_on_manual.Rds ($cells$clone_final); I rewrote the join key based on dna_bam/rna_bam path but surely theres easier  --> I will also still send you the .Rds files. 
Key finding: naive "best cluster per cell" assignment is badly biased by cluster size imbalance (trunk clusters can be 30-70x larger than subclonal ones)

Even more important finding: single-cell coverage is too shallow to reliably detect subclonal (non-trunk) mutations at the individual-cell level — strict test (≥2 independent alt-supporting SNPs within one subclonal cluster) returns zero or near-zero hits in most patients, because subclonal clusters are a small minority of the SNP whitelist and per-cell coverage of that minority is just incredibly thin. 
