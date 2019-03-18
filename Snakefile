__author__ = "Joey Estabrook"
__email__ = "estabroj@ohsu.edu"
__license__ = "MIT"

"""Computation Hub omic data processing pipeline for SE Bulk RNA-sequencing"""

import datetime
import sys
import os
import pandas as pd
import json

timestamp = ('{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now()))

configfile:"omic_config.yaml"
project_id = config["project_id"]

SAMPLES, = glob_wildcards("samples/raw/{sample}.fastq.gz")

ext = ['r','pdf','xls']
fastqscreen_ext = ['html','png','txt']
insertion_and_clipping_prof_ext = ['r','pdf','xls']
inner_distance_ext = ['_freq.txt','_plot.pdf','_plot.r','.txt']
read_dist_ext = ['txt']
read_gc_ext = ['.xls','_plot.r','_plot.pdf']

with open('cluster.json') as json_file:
    json_dict = json.load(json_file)

rule_dirs = list(json_dict.keys())
rule_dirs.pop(rule_dirs.index('__default__'))

for rule in rule_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'logs',rule)):
        log_out = os.path.join(os.getcwd(), 'logs', rule)
        os.makedirs(log_out)
        print(log_out)

result_dirs = ['diffexp','tables']
for rule in result_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'results',rule)):
        log_out = os.path.join(os.getcwd(), 'results', rule)
        os.makedirs(log_out)
        print(log_out)

def message(mes):
    sys.stderr.write("|--- " + mes + "\n")

def format_plot_columns():
    factors = config['meta_columns_to_plot'].keys()
    reformat_factors = '"' + '","'.join(factors) + '"'
    return 'c({})'.format(reformat_factors)

def get_deseq2_threads(wildcards=None):
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(config["omic_meta_data"]) < 100 or few_coeffs else 6

def get_contrast(wildcards):
    """Return each contrast provided in the configuration file"""
    return config["diffexp"]["contrasts"][wildcards.contrast]

for sample in SAMPLES:
    message("Sample " + sample + " will be processed")

rule all:
    input:
        expand("results/tables/{project_id}_STAR_mapping_statistics.txt", project_id = config['project_id']),
        expand("samples/fastqc/{sample}/{sample}_t_fastqc.zip", sample = SAMPLES),
        expand("samples/fastqscreen/{sample}/{sample}_t_screen.{fastqscreen_ext}", sample=SAMPLES, fastqscreen_ext=fastqscreen_ext),
        expand("rseqc/insertion_profile/{sample}/{sample}.insertion_profile.{ext}",sample=SAMPLES, ext=insertion_and_clipping_prof_ext),
        expand("rseqc/clipping_profile/{sample}/{sample}.clipping_profile.{ext}", sample = SAMPLES, ext = insertion_and_clipping_prof_ext),
        expand("rseqc/read_distribution/{sample}/{sample}.read_distribution.{ext}", sample = SAMPLES, ext = read_dist_ext),
        expand("rseqc/read_GC/{sample}/{sample}.GC{ext}", sample = SAMPLES, ext = read_gc_ext),
        "results/tables/{}_Normed_with_Ratio_and_Abundance.txt".format(config['project_id']),
        "results/diffexp/pca.pdf",
        expand(["results/diffexp/GOterms/{contrast}.diffexp.downFC.2.adjp.0.01_BP_GO.txt", "results/diffexp/GOterms/{contrast}.diffexp.upFC.2.adjp.0.01_BP_GO.txt", "results/diffexp/GOterms/{contrast}.diffexp.downFC.2.adjp.0.01.BP.pdf", "results/diffexp/GOterms/{contrast}.diffexp.upFC.2.adjp.0.01.BP.pdf","results/diffexp/GOterms/{contrast}.diffexp.downFC.2.adjp.0.01_BP_classic_5_all.pdf","results/diffexp/GOterms/{contrast}.diffexp.upFC.2.adjp.0.01_BP_classic_5_all.pdf"], contrast = config["diffexp"]["contrasts"]),
        expand("results/diffexp/{contrast}.diffexp.01.VolcanoPlot.pdf", contrast = config["diffexp"]["contrasts"]),
        expand("results/diffexp/permutationTest/Histogram.{contrast}.Permutation.Test.pdf", contrast = config["diffexp"]["contrasts"])

include: "rules/align_rmdp.smk"
include: "rules/omic_qc.smk"
include: "rules/deseq.smk"
