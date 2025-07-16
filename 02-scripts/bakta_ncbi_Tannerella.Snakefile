################################################################################
# Run Bakta on NCBI Tannerella forsythia/serpentiformis genomes to get 
# annotations in the same format as for the ancient MAGs
#
# Irina Velsko, 04/12/2023
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/RIII_simple/04-analysis/bakta/"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/RIII_simple/04-analysis/bakta/input/*.fna"):
	SAMPLES[os.path.basename(sample).split(".b")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("{sample}.done.txt", sample=SAMPLES.keys())

rule damage:
    output:
        touch("{sample}.done.txt")
    message: "Run {wildcards.sample} through Bakta"
    conda: "ENVS_bakta.yaml"
    params: 
        infile = lambda wildcards: SAMPLES[wildcards.sample],
        db = "/mnt/archgen/microbiome_paleobiotech/calcBGCecoevo/03-data/refdbs/bakta/db/"
    threads: 8
    shell:
        """
        bakta --db {params.db} \
        --keep-contig-headers \
        --meta \
        --skip-crispr \
        --threads {threads} \
        {params.infile}
#         {wildcards.sample} 
        """

