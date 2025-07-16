################################################################################
# Run SourceTracker on a rarefied species table for the RIII data
#
# Irina Velsko 02/01/2023
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/RIII_simple/04-analysis/sourcetracker/by_library"

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input:
        "sourcetracker_by_library.done"
        
rule sourcetracker:
    output:
        "sourcetracker_by_library.done"
    message: "Run SourceTracker for the RIII samples at species level"
    params:
        otu = "/mnt/archgen/microbiome_calculus/RIII_simple/05-results/sourcetracker_input_table_w_blanks_rarefied_by_library.tsv",
        map = "/mnt/archgen/microbiome_calculus/RIII_simple/00-documentation/sourcetracker_mapping_file_w_blanks_riii_by_library.tsv"
    shell:
        """
        Rscript \
        /projects1/users/velsko/bin/sourcetracker-1.0.1/sourcetracker_for_qiime.r \
        -i {params.otu} \
        -m {params.map} \
        -o shotgun_sourcetracker_species \
        -r 10000 \
        --train_rarefaction 5000 \
        -v
        """

