################################################################################
# Calculate GC content and read length of each contig in the coverage-filtered abot439
# manual bins to see if the ones with too much length have 2 different distributions
#
# Irina Velsko, 13/07/2023
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/RIII_simple/04-analysis/gc_rl"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/RIII_simple/03-preprocessing/by_sample/eager2_out/samtools/filter/*.fastq"):
	SAMPLES[os.path.basename(sample).split(".f")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("{sample}.gc_rl.tsv.gz", sample=SAMPLES.keys())

rule gc_rl:
    output:
        "{sample}.gc_rl.tsv.gz"
    message: "Run emboss infoseq on {wildcards.sample}"
    params: 
        fasta = lambda wildcards: SAMPLES[wildcards.sample]
    shell:
        """
        infoseq -auto -outfile {wildcards.sample}.gc_rl.tsv -only -name -length -pgc {params.fasta}
        pigz -p 8 {wildcards.sample}.gc_rl.tsv
        """
        
