################################################################################
# Align ancient calculus against Tannerella serpentiformis representative genome 
#
# Irina Velsko 30/11/2023
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/RIII_simple/04-analysis/tannerella_mapping"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/RIII_simple/04-analysis/tannerella_mapping/input/G12*.gz"):
	SAMPLES[os.path.basename(sample).split(".u")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("map_out/{sample}.Ts.mapped_q25l30_rmdup.bam.bai", sample=SAMPLES.keys()),
        expand("map_out/{sample}.Ts.mapped_q25l30_rmdup.cov", sample=SAMPLES.keys())

rule bwa_aln:
    output:
        temp("map_out/{sample}.sai")
    message: "Align sample {wildcards.sample} against the T. serpentiformis rep genome using BWA aln"
    params: 
        reffa = "/mnt/archgen/microbiome_calculus/RIII_simple/01-data/T_serpentiformis/GCA_003033925.1.fna",
        fastq = lambda wildcards: SAMPLES[wildcards.sample]
    threads: 4
    shell:
        """
        bwa aln -n 0.01 -o 2 -l 16500 -t {threads} \
            {params.reffa} \
            {params.fastq} > {output}
        """

rule bwa_samse:
    input:
        "map_out/{sample}.sai"
    output:
        "map_out/{sample}.Ts.mapped_q25l30.bam"
    message: "Generate alignment file for sample {wildcards.sample} against T. serpentiformis rep genome"
    params:
        reffa = "/mnt/archgen/microbiome_calculus/RIII_simple/01-data/T_serpentiformis/GCA_003033925.1.fna",
        fastq = lambda wildcards: SAMPLES[wildcards.sample]
    shell:
        """
        bwa samse \
            {params.reffa} \
            {input} \
            {params.fastq} | \
        samtools view -Sb -q 25 -F 4 -e 'qlen >= 30' - | samtools sort - -o {output} 
        """

rule samtools_rmdup:
    input:
        "map_out/{sample}.Ts.mapped_q25l30.bam"
    output:
        "map_out/{sample}.Ts.mapped_q25l30_rmdup.bam"
    message: "Remove duplicate mapped reads for sample {wildcards.sample} against the T. serpentiformis rep genome"
    params:
    shell:
        """
        samtools rmdup -s {input} {output}
        """

rule samtools_index:
    input:
        "map_out/{sample}.Ts.mapped_q25l30_rmdup.bam"
    output:
        "map_out/{sample}.Ts.mapped_q25l30_rmdup.bam.bai"
    message: "Index bam file for sample {wildcards.sample} mapped against the T. serpentiformis rep genome"
    params:
    shell:
        """
        samtools index {input}
        """

rule bedtools_coverage:
    input:
        "map_out/{sample}.Ts.mapped_q25l30_rmdup.bam"
    output:
        "map_out/{sample}.Ts.mapped_q25l30_rmdup.cov"
    message: "Calculate depth/breadth of coverage of {wildcards.sample} mapped against each genome"
    params:
        refbed = "/mnt/archgen/microbiome_calculus/RIII_simple/01-data/T_serpentiformis/GCA_003033925.1.bed",
    shell:
        """
        bedtools coverage -a {params.refbed} -b {input} -hist > {output}
        """
