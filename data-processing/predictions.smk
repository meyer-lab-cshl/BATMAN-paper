# how to run from ./:
# conda activate snakemake
# snakemake -s predictions.smk -p --profile uge

import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("7.30.0")


##### functions #####

def construct_combine_input(wildcards):
    fastadir = checkpoints.split_fasta_files.get(**wildcards).output[0]
    print(fastadir)
    mhc=glob_wildcards(os.path.join(fastadir,
                                           "peptides_{mhc}_{plength}.fasta")).mhc
    plength=glob_wildcards(os.path.join(fastadir,
                                           "peptides_{mhc}_{plength}.fasta")).plength
    return expand("data/predictions/netMHCpan_{mhc}_{plength}_formated.csv",
                zip,
                mhc=mhc,
                plength=plength
                )

##### wildcards #####
wildcard_constraints:
    plength="\d+"

##### rules #####

rule all:
    input:
        "data/predictions/netMHCpan_formated.csv",
        "data/sanitycheck_peptides_predictions.txt"

checkpoint split_fasta_files:
    input:
        "data/netmhcpan_inputs_mhc_i.csv"
    output:
        directory("data/fastafiles")
    resources:
        mem_mb = 10000
    conda:
        "envs/peptides.yaml"
    script:
        "scripts/01_generate-fasta.R"

rule netmhcpan:
    input:
        peptides="data/fastafiles/peptides_{mhc}_{plength}.fasta",
    output:
        "data/predictions/netMHCpan_{mhc}_{plength}.txt",
    resources:
        mem_mb = 10000
    params:
        plength=9
    shell:
        """
        netMHCpan -f {input.peptides} \
            -l {wildcards.plength} \
            -a  {wildcards.mhc} \
            -BA \
            > {output}
        """

rule format:
    input:
        "data/predictions/netMHCpan_{mhc}_{plength}.txt",
    output:
        "data/predictions/netMHCpan_{mhc}_{plength}_formated.csv",
    conda:
        "envs/format.yaml"
    resources:
        mem_mb = 5000
    shell:
        """
        Rscript scripts/02_format-predictions.R  \
            --dirpred data/predictions \
            --dirseq data/fastafiles \
            --sampleid {wildcards.mhc}_{wildcards.plength}
        """

rule combine:
    input:
        construct_combine_input
    output:
        "data/predictions/netMHCpan_formated.csv",
    resources:
        mem_mb = 1000
    shell:
        """
        echo "Pos,MHC,Peptide,Core,Of,Gp,Gl,Ip,Il,Icore,Identity,Score_EL,Rank_EL,Score_BA,Rank_BA,Aff_nM,BindLevel"> {output}
        for file in {input}; do
            tail -n +2 $file >> {output}
        done
        """

rule sanitycheck_io:
    input:
        pred="data/netmhcpan_inputs_mhc_i.csv",
        pep="data/predictions/netMHCpan_formated.csv",
    output:
        "data/sanitycheck_peptides_predictions.txt"
    shell:
        """
        pep=`wc -l {input.pep}`
        pred=`wc -l {input.pred}`
        echo "$pep\t$pred" >> {output}.tmp
        awk  '{{ print $3, $1, $2, ($2-1) - ($1-1)*2}}' {output}.tmp > {output}
        rm {output}.tmp
        """

