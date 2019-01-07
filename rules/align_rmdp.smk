rule trimming:
    input:
        "samples/raw/{sample}.fq"
    output:
        "samples/trimmed/{sample}_t.fq"
    log:
        "logs/trimming/{sample}_trimming.log"
    params:
        adapter = config["adapter-SE"]
    conda:
        "../envs/trim.yaml"
    message:
        """--- Trimming."""
    shell:
        """trimmomatic SE -phred33 {input} {output} ILLUMINACLIP:{params.adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"""

rule fastqc:
    input:
        "samples/trimmed/{sample}_t.fq"
    output:
        "samples/fastqc/{sample}/{sample}_t_fastqc.zip"
    log:
        "logs/fastqc/{sample}_fastqc.log"
    conda:
        "../envs/omic_qc_wf.yaml"
    message:
        """--- Quality check of raw data with Fastqc."""
    shell:
        """fastqc --outdir  samples/fastqc/{wildcards.sample} --extract  -f fastq {input}"""

rule fastqscreen:
    input:
        fwd = "samples/trimmed/{sample}_R1_t.fq",
        rev = "samples/trimmed/{sample}_R2_t.fq"
    output:
        "samples/fastqscreen/{sample}/{sample}_R1_t_screen.html",
        "samples/fastqscreen/{sample}/{sample}_R1_t_screen.png",
        "samples/fastqscreen/{sample}/{sample}_R1_t_screen.txt",
        "samples/fastqscreen/{sample}/{sample}_R2_t_screen.html",
        "samples/fastqscreen/{sample}/{sample}_R2_t_screen.png",
        "samples/fastqscreen/{sample}/{sample}_R2_t_screen.txt"
    params:
        conf = config["conf"]
    conda:
        "../envs/fastqscreen.yaml"
    shell:
        """fastq_screen --aligner bowtie2 --conf {params.conf} --outdir samples/fastqscreen/{wildcards.sample} {input.fwd} {input.rev}"""

rule STAR:
    input:
        "samples/trimmed/{sample}_t.fq"
    output:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam",
        "samples/star/{sample}_bam/ReadsPerGene.out.tab",
        "samples/star/{sample}_bam/Log.final.out"
    threads: 12
    params:
        gtf=config["gtf_file"]
    log:
        "logs/star/{sample}_star.log"
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

rule bam_statistics:
    input:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam"
    output:
        "samples/bamstats/{sample}/genome_coverage.json"
    run:
        bamstats=config["bamstats_tool"]
        gtf = config["gtf_file"]

        shell("{bamstats} -a {gtf} -i {input} -o {output} -u")

rule get_bam_coverage:
    input:
        expand("samples/bamstats/{sample}/genome_coverage.json", sample=SAMPLES)
    output:
        "data/{project_id}_coverage.txt".format(project_id=config["project_id"])
    script:
        "../scripts/get_coverage.py"

rule genecount:
    input:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam"
    output:
        "samples/htseq_count/{sample}_htseq_gene_count.txt",
    log:
        "logs/genecount/{sample}_genecount.log"
    params:
        name = "genecount_{sample}",
        gtf = config["gtf_file"]
    conda:
        "../envs/omic_qc_wf.yaml"
    threads: 1
    shell:
        """
          htseq-count \
                -f bam \
                -r name \
                -s reverse \
                -m union \
                {input} \
                {params.gtf} > {output}"""

rule count_exons:
    input:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam"
    output:
        "samples/htseq_exon_count/{sample}_htseq_exon_count.txt"
    params:
        exon_gtf = config["exon_gtf"]
    conda:
        "../envs/omic_qc_wf.yaml"
    shell:
        """htseq-count \
                -f bam \
                -m intersection-nonempty \
                -i exon_id \
                --additional-attr=gene_name \
                {input} \
                {params.exon_gtf} > {output}"""

rule compile_counts:
    input:
        expand("samples/htseq_count/{sample}_htseq_gene_count.txt",sample=SAMPLES)
    output:
        "data/{project_id}_counts.txt".format(project_id=config["project_id"])
    script:
        "../scripts/compile_counts_table.py"


rule compile_counts_and_stats:
    input:
        expand("samples/htseq_count/{sample}_htseq_gene_count.txt",sample=SAMPLES)
    output:
        "data/{project_id}_counts_w_stats.txt".format(project_id=config["project_id"])
    script:
        "../scripts/compile_counts_table_w_stats.py"


rule compile_exon_counts:
    input:
        expand("samples/htseq_exon_count/{sample}_htseq_exon_count.txt", sample=SAMPLES)
    output:
        "data/{project_id}_exon_counts.txt".format(project_id = config["project_id"])
    conda:
        "../envs/junction_counts.yaml"
    script:
        "../scripts/compile_exon_counts.R"

