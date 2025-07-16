################################################################################
# Get a  list of Tannerella genomes from NCBI
#
# Irina Velsko 24/11/2023
################################################################################

from glob import glob
import os
import shutil

import pandas as pd


if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")
    

rule all:
    input:
        "01-data/T_forsythia/ncbi/ncbi_assembly_t_forsythia_taxids.txt"

rule extract_Eubacteriales_XIII:
    output:
        "01-data/T_forsythia/ncbi/ncbi_assembly_t_forsythia_taxids.txt"
    message: "Extract the information on the AssemblyAccession, Taxid, and SpeciesName from NCBI Assembly"
    conda: "ENVS_entrez.yaml"
    shell:
        """
        esearch -db assembly -query '"Tannerella forsythia"[Organism] AND (latest[filter] AND all[filter] NOT anomalous[filter])' | \
        esummary | xtract -pattern DocumentSummary -element AssemblyAccession,SpeciesName,Taxid > {output}
        """

