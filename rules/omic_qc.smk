
rule insertion_profile:
    input:
        "samples/xeno/{sample}/Filtered_bams/Aligned.sortedByCoord.out_Filtered.bam"
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
        "samples/xeno/{sample}/Filtered_bams/Aligned.sortedByCoord.out_Filtered.bam"
    output:
        "rseqc/clipping_profile/{sample}/{sample}.clipping_profile.r",
        "rseqc/clipping_profile/{sample}/{sample}.clipping_profile.pdf",
        "rseqc/clipping_profile/{sample}/{sample}.clipping_profile.xls",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "clipping_profile.py -i {input} -s SE -o rseqc/clipping_profile/{wildcards.sample}/{wildcards.sample}"


#rule read_distribution:
#    input:
#        "samples/xeno/{sample}/Filtered_bams/Aligned.sortedByCoord.out_Filtered.bam"
#    params:
#        bed=config['bed_file']
#    output:
#        "rseqc/read_distribution/{sample}/{sample}.read_distribution.txt",
#    conda:
#        "../envs/rseqc.yaml"
#    shell:
#        "read_distribution.py -i {input} -r {params.bed} > {output}"


rule read_GC:
    input:
        "samples/xeno/{sample}/Filtered_bams/Aligned.sortedByCoord.out_Filtered.bam"
    output:
        "rseqc/read_GC/{sample}/{sample}.GC.xls",
        "rseqc/read_GC/{sample}/{sample}.GC_plot.r",
        "rseqc/read_GC/{sample}/{sample}.GC_plot.pdf",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_GC.py -i {input} -o rseqc/read_GC/{wildcards.sample}/{wildcards.sample}"
