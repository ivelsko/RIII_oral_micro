################################################################################
# Generate a phylogenetic tree of the RIII + other MAGs and reference genomes with Phylophlan
#
# Irina Velsko 24/11/2023
################################################################################

from glob import glob
import os
import shutil

import pandas as pd

# workdir: "/mnt/archgen/microbiome_calculus/RIII_simple/04-analysis/phylophlan"

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")
    
#### SAMPLES ###################################################################
SAMPLES, = glob_wildcards("05-results/genomes/04-analysis/phylophlan/{sample}.fa.gz")
################################################################################

#### Auxilliary functions ######################################################

def list_Tf_accids(wildcards):
    filepaths = [line.rstrip()
                 for line in open(checkpoints.extract_genus_Tf_urls.get(**wildcards).output[0], "rt")]
    accid_list = ["_".join(os.path.basename(fn).split("_")[:2]) for fn in filepaths]
    return [f"03-data/refgenomes/04-analysis/phylophlan/{accid}.fna.gz" for accid in accid_list]


def determine_assembly_url(wildcards):
    filepaths = {"_".join(os.path.basename(line).split("_")[:2]): line.rstrip().replace("_assembly_report.txt", "_genomic.fna.gz")
                 for line in open(checkpoints.extract_genus_Tf_urls.get(**wildcards).output[0], "rt")}
    return filepaths[wildcards.accid]


def return_fasta_fn(wildcards):
    mag_fasta_fns = []
    with open(checkpoints.fastani_samplelist.get(**wildcards).output[0], "rt") as infile:
        for line in infile:
            mag_fasta_fns.append(line.rstrip())
    return mag_fasta_fns


################################################################################

rule all:
    input:
        "05-results/PHYL_T_forsythia_proteintree_RAxML.tre",
        "05-results/PHYL_T_forsythia_proteintree_taxa.tsv"

rule extract_T_forsythia_taxids:
    output:
        "04-analysis/phylophlan/ncbi_assembly_tf_taxids.txt"
    message: "Extract the information on the AssemblyAccession, Taxid, and SpeciesName from NCBI Assembly"
    conda: "ENVS_entrez.yaml"
    shell:
        """
        esearch -db assembly -query '"Tannerella forsythia"[Organism] AND (latest[filter] AND all[filter] NOT anomalous[filter])' | \
        esummary | xtract -pattern DocumentSummary -element AssemblyAccession,SpeciesName,Taxid > {output}
        """

rule phylophlan3_write_config:
    output:
        "04-analysis/phylophlan/config.cfg"
    message: "Write the config file for PhyloPhlAn3"
    conda: "ENVS_PhyloPhlAn3.yaml"
    shell:
        """
        phylophlan_write_config_file \
            -d a \
            -o {output} \
            --db_aa diamond \
            --map_dna diamond \
            --map_aa diamond \
            --msa mafft \
            --trim trimal \
            --tree1 fasttree \
            --tree2 raxml
        """

rule setup_T_forsythia_coregenes:
    output:
        "tmp/Tf_tree/core_genes/s__Tannerella_forsythia/s__Tannerella_forsythia.faa"
    message: "Download the core genes of Tannerella forsythia and generate a database"
    conda: "ENVS_PhyloPhlAn3.yaml"
    params:
        dir = "tmp/Tf_tree/core_genes"
    shell:
        """
        mkdir -p {params.dir} && \
        phylophlan_setup_database \
            -g s__Tannerella_forsythia --database_update \
            -o {params.dir} \
            --verbose
        """

rule phylophlan3_tree:
    input:
        database = "tmp/Tf_tree/core_genes/s__Tannerella_forsythia/s__Tannerella_forsythia.faa",
        config = "04-analysis/phylophlan/config.cfg"
    output:
        "04-analysis/phylophlan/RAxML_bestTree.tf_mags_refs_refined.tre"
    message: "Run PhyloPhlAn3 to generate a tree of T. forsythia, with T. serpentiformis for an outgoup"
    conda: "ENVS_PhyloPhlAn3.yaml"
    params:
        genomes = "01-data/tf_mags_refs",
        db_folder = "tmp/Tf_tree/core_genes",
        outdir = "04-analysis/phylophlan"
    log: "04-analysis/phylophlan/phylophlan.log"
    threads: 36
    shell:
        """
        mkdir -p {params.genomes}/phylophlan
        phylophlan \
            -i {params.genomes} \
            -d s__Tannerella_forsythia \
            --databases_folder {params.db_folder} \
            --diversity medium \
            --accurate \
            --genome_extension .fa \
            -f {input.config} \
            -o T_forsythia \
            --output_folder {params.outdir} \
            --nproc {threads} \
            --verbose 2>&1 | tee {log}
        for fn in tf_mags_refs_concatenated.aln tf_mags_refs.tre tf_mags_refs_resolved.tre \
                  RAxML_bestTree.tf_mags_refs_refined.tre RAxML_info.tf_mags_refs_refined.tre RAxML_result.tf_mags_refs_refined.tre; do
            mv 04-analysis/phylophlan/T_forsythia/${{fn}} {params.outdir}/
        done
        """

rule summarise_results:
    input:
        taxids = "04-analysis/phylophlan/ncbi_assembly_tf_taxids.txt",
        drep = "04-analysis/phylophlan/fake_dRep/Wdb.csv",
        tree = "04-analysis/phylophlan/RAxML_bestTree.tf_mags_refs_refined.tre"
    output:
        taxa = "05-results/PHYL_T_forsythia_proteintree_taxa.tsv",
        tree = "05-results/PHYL_T_forsythia_proteintree_RAxML.tre"
    message: "Summarise the dereplication of the taxa and provide the RAxML tree"
    run:
        # List of downloaded genomes
        taxids = pd.concat([pd.read_csv(input.taxids, sep="\t", header=None,
                                        names=['accession Id', 'name', 'NCBI taxonomy Id']),
                            pd.DataFrame.from_dict({'accession Id': "GCA_003033925.1",
                                                    'name': "Tannerella serpentiformis",
                                                    'NCBI taxonomy Id': 712710}, orient="index") \
                            .transpose()])

        # Dereplication results
        drep_report = pd.read_csv(input.drep, sep=",")
        drep_report['genome'] = drep_report['genome'].str.replace(".fna", "")
        taxids = taxids.loc[(taxids['accession Id'].isin(drep_report['genome'].tolist())) |
                            (taxids['accession Id'] == "GCA_003033925.1")]

        taxids.to_csv(output.taxa, sep="\t", index=False)
        shutil.copy2(input.tree, output.tree)
