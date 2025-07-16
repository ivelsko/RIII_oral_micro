################################################################################
# Run dRep on Tannerella forsythia genomes
#
# Irina Velsko, 28/11/2023
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/RIII_simple/04-analysis/dRep"


if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        "dRep_40_Tf_refs_out/data_tables/Ndb.csv"

rule dRep:
    output:
        "dRep_40_Tf_refs_out/data_tables/Ndb.csv"
    message: "Run dRep to cluster Tannerella forsythia genomes NCBI"
    params: 
        genomes = "/mnt/archgen/microbiome_calculus/RIII_simple/04-analysis/dRep/tf_refs_list.tsv"
    threads: 12
    shell:
        """
        dRep dereplicate -p {threads} dRep_40_Tf_refs_out/ -g {params.genomes} -comp 40 --S_algorithm ANImf -pa 0.95 -sa 0.99
        """
     
