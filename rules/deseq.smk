contrast = get_contrast

rule deseq2_init:
    input:
        counts = "data/{project_id}_counts.filt.txt".format(project_id=config["project_id"])
    output:
        rds="results/diffexp/pairwise/{contrast}_all.rds",
        rld_out = "results/diffexp/pairwise/{contrast}_rlog_dds.rds"
    params:
        samples=config["omic_meta_data"],
        sample_id = config["sample_id"],
        linear_model = config["linear_model"],
        contrast = get_contrast
    conda:
        "../envs/permutation.yaml"
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2-init.R"

rule deseq2_pairwise:
    input:
        rds="results/diffexp/pairwise/{contrast}_all.rds",
        rld = "results/diffexp/pairwise/{contrast}_rlog_dds.rds"
    output:
        table="results/diffexp/pairwise/{contrast}.diffexp.tsv",
        ma_plot="results/diffexp/pairwise/{contrast}.ma_plot.pdf",
        p_hist="results/diffexp/pairwise/{contrast}.phist_plot.pdf",
        heatmap_plot = "results/diffexp/pairwise/{contrast}.heatmap_plot.pdf",
        panel_ma = "results/diffexp/pairwise/{contrast}.panel_ma.pdf",
        var_heat = "results/diffexp/pairwise/{contrast}.variance_heatmap.pdf",
        pca_plot = "results/diffexp/pairwise/{contrast}.pca_plot.pdf"
    params:
        contrast=get_contrast,
        linear_model = config["linear_model"],
        pca_labels = config["pca"]["labels"],
        sample_id = config["linear_model"]
    conda:
        "../envs/deseq2.yaml"
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2_pairwise.R"

rule deseq2_group:
    input:
        counts = "data/{project_id}_counts.filt.txt".format(project_id = project_id)
    output:
        pca="results/diffexp/group/LRT_pca.pdf",
        sd_mean_plot="results/diffexp/group/LRT_sd_mean_plot.pdf",
        distance_plot = "results/diffexp/group/LRT_distance_plot.pdf",
        heatmap_plot = "results/diffexp/group/LRT_heatmap_plot.pdf",
        rds="results/diffexp/group/LRT_all.rds",
        rld_out = "results/diffexp/group/LRT_rlog_dds.rds"
    params:
        pca_labels=config["pca"]["labels"],
        samples=config["omic_meta_data"],
        sample_id = config["sample_id"],
        linear_model = config["linear_model"],
        LRT = config["diffexp"]["LRT"]
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/deseq2_group.R"

rule GO:
    input:
        "results/diffexp/pairwise/{contrast}.diffexp.tsv"
    output:
        "results/diffexp/pairwise/GOterms/{contrast}.diffexp.downFC.2.adjp.0.01_BP_GO.txt",
        "results/diffexp/pairwise/GOterms/{contrast}.diffexp.upFC.2.adjp.0.01_BP_GO.txt",
        "results/diffexp/pairwise/GOterms/{contrast}.diffexp.downFC.2.adjp.0.01_BP_classic_5_all.pdf",
        "results/diffexp/pairwise/GOterms/{contrast}.diffexp.upFC.2.adjp.0.01_BP_classic_5_all.pdf"
    params:
        contrast = get_contrast,
        assembly = config["assembly"]
    conda:
        "../envs/runGO.yaml"
    shell: """Rscript scripts/runGOforDESeq2.R --degFile={input} --assembly={params.assembly} --FC=2 --adjp=0.01 --printTree=1"""

rule volcano:
    input:
        "results/diffexp/pairwise/{contrast}.diffexp.tsv"
    output:
        "results/diffexp/pairwise/{contrast}.diffexp.01.VolcanoPlot.pdf"
    params:
       contrast = get_contrast
    shell: """Rscript scripts/RNAseq_makeVolcano.R --degFile={input} --adjp=0.01 --FC=2"""

rule permutation:
    input:
        counts = "data/{project_id}_counts.filt.txt".format(project_id=config["project_id"])
    output:
        numGenes = "results/diffexp/pairwise/permutationTest/{contrast}.number.diff.genes.csv",
        permList = "results/diffexp/pairwise/permutationTest/{contrast}.permutation.list.csv",
        histogram = "results/diffexp/pairwise/permutationTest/Histogram.{contrast}.Permutation.Test.pdf"
    params:
        contrast = get_contrast,
        samples = config["omic_meta_data"],
        sample_id = config["sample_id"],
        linear_model = config["linear_model"]
    conda:
        "../envs/permutation.yaml"
    script:
        "../scripts/permutation_test.R"

