################################################################################
# Run DamageProfiler on RIII, etc libraries mapped against the T. forsythia genome
#
# Irina Velsko, 30/11/2023
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/RIII_simple/04-analysis/DamageProfiler/"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/RIII_simple/04-analysis/tannerella_mapping/map_out/*rmdup.bam"):
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
    message: "Run {wildcards.sample} through DamageProfiler"
    params: 
        infile = lambda wildcards: SAMPLES[wildcards.sample],
    threads: 32
    shell:
        """
        set +u
        source $HOME/miniconda3/etc/profile.d/conda.sh
        conda activate damageprofiler
        set -u

        java -jar /home/irina_marie_velsko/miniconda3/envs/damageprofiler/share/damageprofiler-1.1-2/DamageProfiler-1.1-java11.jar \
        -i {params.infile}\
        -o out/{wildcards.sample}
        """

