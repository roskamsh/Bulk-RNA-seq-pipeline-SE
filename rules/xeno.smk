rule genecount:
    input:
        "samples/xeno/{sample}/Filtered_bams/Aligned.sortedByCoord.out_Filtered.bam"
    output:
        "samples/htseq/{sample}_htseq_gene_count.txt"
    log:
        "logs/genecount/{sample}_genecount.log"
    params:
        name = "genecount_{sample}",
        gtf = config["gtf_file_graft"]
    conda:
        "../envs/omic_qc_wf.yaml"
    shell:
        """htseq-count \
                -f bam \
                -r name \
                -s reverse \
                -m intersection-strict \
                {input} \
                {params.gtf} > {output}"""

rule compile_counts:
    input:
        expand("samples/htseq/{sample}_htseq_gene_count.txt",sample=SAMPLES)
    output:
        "data/{project_id}_counts.txt".format(project_id=config["project_id"])
    script:
        "../scripts/compile_counts_table.py"

rule filter_counts:
    input:
        countsFile="data/{project_id}_counts.txt".format(project_id=config["project_id"])
    output:
        "data/{project_id}_counts.filt.txt".format(project_id=config["project_id"])
    params:
        anno=config["filter_anno"],
        biotypes=config["biotypes"],
        mito=config['mito']
    script:
        "../scripts/RNAseq_filterCounts.R"

