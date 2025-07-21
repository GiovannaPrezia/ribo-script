configfile: "config.yaml"

include: "rules/download.smk"
include: "rules/preprocessing.smk"
include: "rules/alignment.smk"
include: "rules/quantification.smk"
include: "rules/ribotricer.smk"
include: "rules/analysis.smk"

rule all:
    input:
        expand("results/alignment/{sample}_Aligned.sortedByCoord.out.bam", sample=config["samples"].values()),
        "results/deseq2_MCF7_vs_MCF10A_with_counts.csv",
        "results/pca_plot.png",
        "results/volcano_plot.png",
        "results/heatmap_top30.png",
        "results/heatmap_genes.csv",
        "results/ribotricer/ORFs_results_stats.tsv"
