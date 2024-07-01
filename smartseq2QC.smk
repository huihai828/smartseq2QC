#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Snakemake pipelines
Smart-seq2 quantification and quality control pipeline, developed for single-cell service at the Earlham Institue.
"""

import sys
import os
import re
import pandas as pd

configfile: "config.yaml"

FASTQ_DIR = config["fastq_dir"]
QUANT_DIR = config["quant_dir"]
QCOUT_DIR = config["qcout_dir"]
PLATE_IDS = config["plate_ids"]
count_type = config["count_type"]

samplesheet = pd.read_csv(config["samplesheet"])
SAMPLES = samplesheet[["Sample_Plate", "Sample_Name", "Sample_ID",]].agg("_".join, axis=1) \
          + "_" + samplesheet[["index", "index2"]].agg("-".join, axis=1) \
          + "_L00" + samplesheet["Lane"].astype(str)



doc_complete_check = expand("Finished_{plate_id}.txt", plate_id = PLATE_IDS)

if len(PLATE_IDS) > 1:
    PLATE_IDS += ['all']

TARGETS = []
TARGETS.append(expand(os.path.join(f"{QUANT_DIR}", "{sample}", "run_info.json"), sample=SAMPLES))
TARGETS.append(expand(f"qc_{count_type}{{plate_id}}_matrix.txt", plate_id = PLATE_IDS))
TARGETS.append(expand("QC_meanexp_vs_freq_{plate_id}.png", plate_id = PLATE_IDS))
TARGETS.append(expand("{plate_id}_QC_report.pdf", plate_id = PLATE_IDS))
TARGETS.append(os.path.join(f"{QUANT_DIR}", f"{count_type}all_matrix.tsv"))
TARGETS.append("tx2gene_finished.txt")


rule all:
    input: TARGETS


rule kallisto_quant:
    input:
        read1 = os.path.join(FASTQ_DIR, "{sample}_R1.fastq.gz"),
        read2 = os.path.join(FASTQ_DIR, "{sample}_R2.fastq.gz")
    output:
        os.path.join(f"{QUANT_DIR}", "{sample}", "run_info.json")
    params:
        idx = config["idx"],
        outdir = os.path.join(f"{QUANT_DIR}", "{sample}")
    threads: 2
    log:
        "logs/kallisto_quant_{sample}.log"
    benchmark:
        "benchmarks/kallisto_quant_{sample}.tsv"
    singularity:
        "/ei/software/cb/singlecellQC/1.1/x86_64/rscater.img"
    shell:
         "kallisto quant -t {threads} -i {params.idx} -o {params.outdir} -b 100 {input.read1} {input.read2} &> {log}"


rule kallisto_mapping_scrape:
    output:
        mapping_stats = os.path.join(f"{QCOUT_DIR}", "percent_pseudoaligned.txt")
    params:
        quant_dir = QUANT_DIR,
        metric = "runinfo_metric"
    threads: 1
    log:
        "logs/kallisto_mapping_scrape.log"
    benchmark:
        "benchmarks/kallisto_mapping_scrape.tsv"
    singularity:
        "/ei/software/cb/singlecellQC/1.1/x86_64/rscater.img"
    shell:
         "Rscript scripts/kallisto_mapping_scrape.R {params.quant_dir} {output.mapping_stats} {params.metric} &> {log}"


rule merge_kallisto_quant:
    output:
        os.path.join(f"{QCOUT_DIR}", f"{count_type}{{plate_id}}_matrix.tsv")
    params:
        quant_dir = QUANT_DIR,
        count_type = count_type
    threads: 1
    log:
        "logs/merge_kallisto_quant_{plate_id}.log"
    benchmark:
        "benchmarks/merge_kallisto_quant_{plate_id}.tsv"
    singularity:
        "/ei/software/cb/singlecellQC/1.1/x86_64/rscater.img"
    shell:
         """
            Rscript scripts/merge_kallisto_quant.R {params.quant_dir} {params.count_type} {wildcards.plate_id} &> {log};
            ln -s -v -f {params.quant_dir}/{params.count_type}{wildcards.plate_id}_matrix.tsv {params.count_type}{wildcards.plate_id}_matrix.tsv
         """


rule qc:
    input:
        mapping_stats = os.path.join(f"{QCOUT_DIR}", "percent_pseudoaligned.txt"),
        plate_matrix = os.path.join(f"{QCOUT_DIR}", f"{count_type}{{plate_id}}_matrix.tsv")
    output:
        os.path.join(f"{QCOUT_DIR}", "QC_meanexp_vs_freq_{plate_id}.pdf"),
        os.path.join(f"{QCOUT_DIR}", "{plate_id}qc_for_doc.Rdata"),
        done = f"qc_{count_type}{{plate_id}}_matrix.txt"
    params:
        qcout_dir = QCOUT_DIR,
        count_type = count_type,
        samplesheet = samplesheet,
        mtnamefile = config["mtnamefile"]
    threads: 1
    log:
        "logs/scqc_from_matrix_{plate_id}.log"
    benchmark:
        "benchmarks/scqc_from_matrix_{plate_id}.tsv"
    singularity:
        "/ei/software/cb/singlecellQC/1.1/x86_64/rscater.img"
    shell:
         """
            Rscript scripts/scqc_from_matrix.meta.R {input.plate_matrix} {input.mapping_stats} {params.qcout_dir} \
            {params.samplesheet} {params.mtnamefile} &> {log};
            echo 'QC of {wildcards.plate_id} is done. ' > {output.done}
         """


rule gs:
    input:
        pdf_file = os.path.join(f"{QCOUT_DIR}", "QC_meanexp_vs_freq_{plate_id}.pdf"),
    output:
        png_file = "QC_meanexp_vs_freq_{plate_id}.png",
        pdf_rename = os.path.join(f"{QCOUT_DIR}", "QC_meanexp_vs_freq_{plate_id}_changename.pdf")
    params:
        qcout_dir = QCOUT_DIR,
    threads: 1
    log:
        "logs/gs_png_{plate_id}.log"
    benchmark:
        "benchmarks/gs_png_{plate_id}.tsv"
    singularity:
        "/ei/software/cb/singlecellQC/1.1/x86_64/rscater.img"
    shell:
         """
            gs -dNOPAUSE -dQUIET -dBATCH -sDEVICE=png16m -sOutputFile={output.png_file} -r256 {input.pdf_file} &> {log};
            mv {input.pdf_file} {output.pdf_rename}
         """


rule doc:
    input:
        mapping_stats = os.path.join(f"{QCOUT_DIR}", "percent_pseudoaligned.txt"),
        qc_rdata = os.path.join(f"{QCOUT_DIR}", "{plate_id}qc_for_doc.Rdata"),
    output:
        done = "Finished_{plate_id}.txt",
        report = "{plate_id}_QC_report.pdf"
    params:
        qcout_dir = QCOUT_DIR,
        out_dir = os.path.join(f"{QCOUT_DIR}", "{plate_id}"),
    threads: 1
    log:
        "logs/doc_qc_report_{plate_id}.log"
    benchmark:
        "benchmarks/doc_qc_report_{plate_id}.tsv"
    singularity:
        "/ei/software/cb/singlecellQC/1.1/x86_64/rknitr.img"
    shell:
         """
            Rscript -e \"options(warn=-1);objects<-\'{input.qc_rdata}\'; \
            mapping_file <- read.table(\'{input.mapping_stats}\'); \
            rmarkdown::render(\'scripts/QCreport.Rmd\', 'pdf_document', output_file=\'{output.report}\', output_dir=\'{params.out_dir}\')\" &> {log};
            echo '{wildcards.plate_id}' > {output.done}
         """


rule plate_merge:
    output:
        all_tsv = os.path.join(f"{QUANT_DIR}", f"{count_type}all_matrix.tsv")
    params:
        quant_dir = QUANT_DIR,
        count_type = count_type,
    threads: 1
    log:
        "logs/plate_merge.log"
    benchmark:
        "benchmarks/plate_merge.tsv"
    singularity:
        "/ei/software/cb/singlecellQC/1.1/x86_64/rscater.img"
    shell:
         """
            Rscript scripts/plate_merge.R {params.quant_dir} {params.count_type} \'{doc_complete_check}\' &> {log};
            ln -s -v -f {output.all_tsv} all_plates.tsv
         """


rule tx2gene:
    output:
        done = "tx2gene_finished.txt",
        quant_gene = os.path.join(f"{QUANT_DIR}", "plates_as_genelevel.tsv")
    params:
        qcout_dir = QCOUT_DIR,
        quant_dir = QUANT_DIR,
        species = config["species"],
        trans2gen_tsv = config["trans2gen_tsv"]
    threads: 1
    log:
        "logs/est_counts_tx2gene.log"
    benchmark:
        "benchmarks/est_counts_tx2gene.tsv"
    singularity:
        "/ei/software/cb/singlecellQC/1.1/x86_64/rscater.img"
    shell:
         """
            Rscript scripts/est_counts_tx2gene.R {params.qcout_dir} ${params.quant_dir} {params.species} {params.trans2gen_tsv} &> {log};
            echo 'tx2g finished' > {output.done}
         """

















