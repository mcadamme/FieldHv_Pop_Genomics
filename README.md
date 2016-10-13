# FieldHv_Pop_Genomics
# The following scripts were used for Hv popgen data analysis

pyrethroidResistance_AnalysisforAlexD.R - for the candidate gene analysis, test for statistically significant associations between candidate gene (voltage-gated sodium channel), and ddRAD-seq marker.
bowtieScript - to align reads to my reference genome using Bowtie2.
stacks_refmap_script - to cluster my reference genome aligned reads.
populations.sh - to get summary data from Stacks, enabling me to see which ddRAD-seq markers aligned to Contig 4600.
shell_script_loc11322_fasta.sh - invokes stacks_to_fasta_matches.py (written by Ryan Waples) to generate a fasta file for each field-collected individual with their Stacks-generated ddRAD-seq marker genotypes at locus 11322. These fasta files were then edited by hand to include indels, remove genotypes where marker counts were low, etc.
samtools_script.sh - convert sam files to bam files that could be fed into samtools mpileup and bcftools for SNP calling.
vcftools_script - looking at contig-specific and genome-wide SNP changes over time.
