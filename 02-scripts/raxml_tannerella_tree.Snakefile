################################################################################
# Run RaxML on the Tannerella forsythia marker gene protein alignment from  Phylophlan 
# to get bootstrap values
#
# Irina Velsko 11/01/2025
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/RIII_simple/04-analysis/phylogenies"

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")

rule all:
    input: 
        "phylophlan_aln.raxml.done",
        "phylophlan_aln.raxml.bs.done",
        "phylophlan_aln.raxml.support.done"

rule raxml:
    output:
        touch("phylophlan_aln.raxml.done")
    message: "Run raxml on phylophlan marker gene protein alignment"
    conda: "ENVS_raxml.yaml"
    params: 
        infa = "../phylophlan/tf_mags_refs_concatenated.aln"
    threads: 8
    shell:
        """
         raxml-ng --msa {params.infa} \
         --model LG \
         --prefix phylophlan_aln \
         --threads {threads} \
         --seed 1989
        """

rule bootstrap:
    input:
        "phylophlan_aln.raxml.done"
    output:
        touch("phylophlan_aln.raxml.bs.done")
    message: "Run raxml on phylophlan marker gene protein alignment"
    conda: "ENVS_raxml.yaml"
    params: 
        infa = "../phylophlan/tf_mags_refs_concatenated.aln"
    threads: 24
    shell:
        """
         raxml-ng --msa {params.infa} \
         --model LG \
         --bootstrap \
         --bs-trees 200 \
         --prefix phylophlan_aln \
         --threads {threads} \
         --seed 1989
        """


# next step for branch support
rule branch_support:
    input: 
        "phylophlan_aln.raxml.bs.done"
    output:
        touch("phylophlan_aln.raxml.support.done")
    message: "Get branch support for RaxML tree"
    conda: "ENVS_raxml.yaml"
    params: 
        tree = "phylophlan_aln.raxml.bestTree",
        boots = "phylophlan_aln.raxml.bootstraps"
    threads: 2
    shell:
        """
        raxml-ng --support \
        --tree {params.tree} \
        --bs-trees {params.boots} \
        --prefix phylophlan_aln \
        --threads {threads} 
        """
