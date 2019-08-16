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


rule star_graft:
    input:
        "samples/trimmed/{sample}_t.fastq"
    output:
        "samples/star_graft/{sample}_bam/Aligned.sortedByCoord.out.bam",
        "samples/star_graft/{sample}_bam/ReadsPerGene.out.tab",
        "samples/star_graft/{sample}_bam/Log.final.out"
    threads: 12
    params:
        gtf=config["gtf_file_graft"]

    run:
         STAR=config["star_tool"],
         pathToGenomeIndex = config["star_index_graft"]

         shell("""
                {STAR} --runThreadN {threads} --runMode alignReads --genomeDir {pathToGenomeIndex} \
                --readFilesIn {input} \
                --outFileNamePrefix samples/star_graft/{wildcards.sample}_bam/ \
                --sjdbGTFfile {params.gtf} --quantMode GeneCounts \
                --sjdbGTFtagExonParentGene gene_name \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMattributes NM \
                --twopassMode Basic
                """)


rule star_host:
    input:
        "samples/trimmed/{sample}_t.fastq"
    output:
        "samples/star_host/{sample}_bam/Aligned.sortedByCoord.out.bam",
        "samples/star_host/{sample}_bam/ReadsPerGene.out.tab",
        "samples/star_host/{sample}_bam/Log.final.out"
    threads: 12
    params:
        gtf=config["gtf_file_host"]
    run:
         STAR=config["star_tool"],
         pathToGenomeIndex = config["star_index_host"]

         shell("""
                {STAR} --runThreadN {threads} --runMode alignReads --genomeDir {pathToGenomeIndex} \
                --readFilesIn {input} \
                --outFileNamePrefix samples/star_host/{wildcards.sample}_bam/ \
                --sjdbGTFfile {params.gtf} --quantMode GeneCounts \
                --sjdbGTFtagExonParentGene gene_name \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMattributes NM \
                --twopassMode Basic
                """)


rule star_statistics_graft:
    input:
        expand("samples/star_graft/{sample}_bam/Log.final.out",sample=SAMPLES)
    output:
        "results/tables/{project_id}_graft_STAR_mapping_statistics.txt".format(project_id = config["project_id"])
    script:
        "../scripts/compile_star_log.py"

rule star_statistics_host:
    input:
        expand("samples/star_host/{sample}_bam/Log.final.out",sample=SAMPLES)
    output:
        "results/tables/{project_id}_host_STAR_mapping_statistics.txt".format(project_id = config["project_id"])
    script:
        "../scripts/compile_star_log.py"
