
rule insertion_profile:
    input:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam"
    output:
        "rseqc/insertion_profile/{sample}/{sample}.insertion_profile.r",
        "rseqc/insertion_profile/{sample}/{sample}.insertion_profile.pdf",
        "rseqc/insertion_profile/{sample}/{sample}.insertion_profile.xls",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "insertion_profile.py -s SE -i {input} -o rseqc/insertion_profile/{wildcards.sample}/{wildcards.sample}"


rule clipping_profile:
    input:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam"
    output:
        "rseqc/clipping_profile/{sample}/{sample}.clipping_profile.r",
        "rseqc/clipping_profile/{sample}/{sample}.clipping_profile.pdf",
        "rseqc/clipping_profile/{sample}/{sample}.clipping_profile.xls",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "clipping_profile.py -i {input} -s SE -o rseqc/clipping_profile/{wildcards.sample}/{wildcards.sample}"


rule read_distribution:
    input:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam"
    params:
        bed=config['bed_file']
    output:
        "rseqc/read_distribution/{sample}/{sample}.read_distribution.txt",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_distribution.py -i {input} -r {params.bed} > {output}"


rule read_GC:
    input:
        "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam"
    output:
        "rseqc/read_GC/{sample}/{sample}.GC.xls",
        "rseqc/read_GC/{sample}/{sample}.GC_plot.r",
        "rseqc/read_GC/{sample}/{sample}.GC_plot.pdf",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_GC.py -i {input} -o rseqc/read_GC/{wildcards.sample}/{wildcards.sample}"


rule generate_qc_qa:
 input:
    counts = "data/{project_id}_counts.filt.txt".format(project_id=config['project_id'])
 params:
    project_id = config["project_id"],
    read_dir = config['base_dir'],
    meta = config["omic_meta_data"],
    baseline = config["baseline"],
    linear_model = config["linear_model"],
    sample_id = config["sample_id"],
    gtf_file = config["gtf_file"],
    meta_viz = format_plot_columns(),
 output:
    "analysis_code/{project_id}_analysis.R".format(project_id=project_id)
 script:
    "../scripts/GenerateAbundanceFile.py"


rule run_qc_qa:
    input:
        rules.generate_qc_qa.output
    output:
        "results/tables/{project_id}_Normed_with_Ratio_and_Abundance.txt".format(project_id=config['project_id'])
    conda:
        "../envs/omic_qc_wf.yaml"
    script:
        "../analysis_code/{project_id}_analysis.R".format(project_id=config['project_id'])
