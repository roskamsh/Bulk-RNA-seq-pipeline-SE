rule trimming:
    input:
        "samples/raw/{sample}.fastq.gz"
    output:
        "samples/trimmed/{sample}_t.fastq"
    params:
        adapter = config["adapter"]
    conda:
        "../envs/trim.yaml"
    message:
        """--- Trimming."""
    shell:
        """trimmomatic SE -phred33 {input} {output} ILLUMINACLIP:{params.adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50"""


rule fastqc:
    input:
        "samples/trimmed/{sample}_t.fastq"
    output:
        "samples/fastqc/{sample}/{sample}_t_fastqc.zip"

    conda:
        "../envs/fastqc.yaml"
    message:
        """--- Quality check of raw data with Fastqc."""
    shell:
        """fastqc --outdir  samples/fastqc/{wildcards.sample} --extract  -f fastq {input}"""

rule fastqscreen:
    input:
        "samples/trimmed/{sample}_t.fastq"
    output:
        "samples/fastqscreen/{sample}/{sample}_t_screen.html",
        "samples/fastqscreen/{sample}/{sample}_t_screen.png",
        "samples/fastqscreen/{sample}/{sample}_t_screen.txt"
    params:
        conf = config["conf"]
    conda:
        "../envs/fastqscreen.yaml"
    shell:
        """fastq_screen --aligner bowtie2 --conf {params.conf} --outdir samples/fastqscreen/{wildcards.sample} {input}"""

rule STAR:
    input:
        "samples/trimmed/{sample}_t.fastq"
    output:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam",
        "samples/star/{sample}_bam/ReadsPerGene.out.tab",
        "samples/star/{sample}_bam/Log.final.out"
    threads: 12
    params:
        gtf=config["gtf_file"]

    run:
         STAR=config["star_tool"],
         pathToGenomeIndex = config["star_index"]

         shell("""
                {STAR} --runThreadN {threads} --runMode alignReads --genomeDir {pathToGenomeIndex} \
                --readFilesIn {input} \
                --outFileNamePrefix samples/star/{wildcards.sample}_bam/ \
                --sjdbGTFfile {params.gtf} --quantMode GeneCounts \
                --sjdbGTFtagExonParentGene gene_name \
                --outSAMtype BAM SortedByCoordinate \
                #--readFilesCommand zcat \
                --twopassMode Basic
                """)

rule index:
    input:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam"
    output:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam.bai"
    conda:
        "../envs/omic_qc_wf.yaml"
    shell:
        """samtools index {input} {output}"""

rule star_statistics:
    input:
        expand("samples/star/{sample}_bam/Log.final.out",sample=SAMPLES)
    output:
        "results/tables/{project_id}_STAR_mapping_statistics.txt".format(project_id = config["project_id"])
    script:
        "../scripts/compile_star_log.py"

rule compile_star_counts:
    input:
        expand("samples/star/{sample}_bam/ReadsPerGene.out.tab",sample=SAMPLES)
    params:
        samples=SAMPLES
    output:
        "data/{project_id}_counts.txt".format(project_id=config["project_id"])
    script:
        "../scripts/compile_star_counts.py"

rule filter_counts:
    input:
        "data/{project_id}_counts.txt".format(project_id=config["project_id"])
    output:
        "data/{project_id}_counts.filt.txt".format(project_id=config["project_id"])
    params:
        anno=config["filter_anno"],
        biotypes=config["biotypes"]
    script:
        "..scripts/RNAseq_filterCounts.R"
