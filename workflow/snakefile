############################################################
# RiboLongShort
# Main Snakemake workflow
############################################################

configfile: "config/config.yaml"

SAMPLES = {
    sample["sample_name"]: sample["run_id"]
    for sample in config["samples"]
}

SIZE_MODES = ["all_lengths", "28_36"]

include: "rules/download.smk"
include: "rules/qc.smk"
include: "rules/trimming.smk"
include: "rules/contaminants.smk"
include: "rules/alignment.smk"
include: "rules/quantification.smk"
include: "rules/reporting.smk"

rule all:
    input:
        expand("02_fastq/ribo_seq/{sample}.fastq.gz", sample=SAMPLES.keys()),
        expand("03_trimmed/ribo_seq/{sample}.trim.fastq.gz", sample=SAMPLES.keys()),
        expand("04_cleaned/ribo_seq/{mode}/{sample}.{mode}.clean.fastq.gz", sample=SAMPLES.keys(), mode=SIZE_MODES),
        expand("05_alignment/ribo_seq/{mode}/{sample}.{mode}_Aligned.sortedByCoord.out.bam", sample=SAMPLES.keys(), mode=SIZE_MODES),
        expand("07_counts/ribo_seq/{mode}/{sample}.{mode}.CDS_counts.txt", sample=SAMPLES.keys(), mode=SIZE_MODES),
        "11_MultiQC/" + config["multiqc_report_name"],
        "11_MultiQC/" + config["qc_summary_name"]
