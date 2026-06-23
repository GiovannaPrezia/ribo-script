#!/bin/bash

set -eu

PROJECT_ROOT="$1"
PROJECT_NAME="$2"

RAW_DIR="$PROJECT_ROOT/02_fastq/ribo_seq"
TRIM_DIR="$PROJECT_ROOT/03_trimmed/ribo_seq"
CLEAN_DIR="$PROJECT_ROOT/04_cleaned/ribo_seq"
STAR_QC_DIR="$PROJECT_ROOT/06_star_qc/ribo_seq"
LOG_DIR="$PROJECT_ROOT/logs"
OUTDIR="$PROJECT_ROOT/13_Report/tables"

mkdir -p "$OUTDIR"

OUTFILE="$OUTDIR/${PROJECT_NAME}_QC_master_table.tsv"

echo -e "Sample\tRun\tMode\tRaw_Reads\tRaw_Sequence_Length\tRaw_%GC\tTrimmed_Reads\t%Removed\tTrimmed_Sequence_Length\tTrimmed_%GC\tTrim_Info\tPost_N_PolyG_Reads\tPost_Contaminants_Reads\tRetention_After_Contaminants_%\tClean_%GC\tClean_Length\tClean_Info\tReads_Post_Contaminants\tUnique_Reads_STAR\tUnique_Mapping_%" > "$OUTFILE"

for RAW_FASTQ in "$RAW_DIR"/*.fastq.gz; do

    SAMPLE=$(basename "$RAW_FASTQ" .fastq.gz)

    RUN=$(grep -R "$SAMPLE" "$PROJECT_ROOT/logs" 2>/dev/null | head -n 1 | awk '{print $2}' || true)
    [[ -z "$RUN" ]] && RUN="NA"

    RAW_STATS=$(seqkit stats -T "$RAW_FASTQ" | awk 'NR==2')
    RAW_READS=$(echo "$RAW_STATS" | awk '{print $4}')
    RAW_LEN=$(echo "$RAW_STATS" | awk '{print $7}')
    RAW_GC=$(echo "$RAW_STATS" | awk '{print $8}')

    TRIM_FASTQ="$TRIM_DIR/${SAMPLE}.trim.fastq.gz"
    TRIM_STATS=$(seqkit stats -T "$TRIM_FASTQ" | awk 'NR==2')
    TRIM_READS=$(echo "$TRIM_STATS" | awk '{print $4}')
    TRIM_LEN=$(echo "$TRIM_STATS" | awk '{print $7}')
    TRIM_GC=$(echo "$TRIM_STATS" | awk '{print $8}')

    CUTADAPT_LOG="$LOG_DIR/${SAMPLE}.cutadapt.log"
    REMOVED=$(grep "Reads written" "$CUTADAPT_LOG" | awk -F'[()]' '{print $2}' | sed 's/%//' || echo "NA")
    TRIM_INFO="Cutadapt_-u3_polyA_A10_min17"

    POST_N_POLYG=$(zcat "$CLEAN_DIR/${SAMPLE}.trim.noN.noPolyG.all_lengths.fastq.gz" | awk 'END{print NR/4}')

    for MODE in all_lengths 28_36; do

        CLEAN_FASTQ="$CLEAN_DIR/$MODE/${SAMPLE}.${MODE}.clean.fastq.gz"
        BOWTIE_LOG="$LOG_DIR/${SAMPLE}.${MODE}.bowtie.log"
        STAR_LOG="$STAR_QC_DIR/${SAMPLE}.${MODE}_Log.final.out"

        if [[ ! -f "$CLEAN_FASTQ" || ! -f "$BOWTIE_LOG" || ! -f "$STAR_LOG" ]]; then
            continue
        fi

        CLEAN_STATS=$(seqkit stats -T "$CLEAN_FASTQ" | awk 'NR==2')
        CLEAN_READS=$(echo "$CLEAN_STATS" | awk '{print $4}')
        CLEAN_LEN=$(echo "$CLEAN_STATS" | awk '{print $7}')
        CLEAN_GC=$(echo "$CLEAN_STATS" | awk '{print $8}')

        PROCESSED=$(grep "# reads processed:" "$BOWTIE_LOG" | awk '{print $4}')
        POST_CONTAM=$(grep "# reads that failed to align:" "$BOWTIE_LOG" | awk '{print $7}')
        RETENTION=$(awk -v clean="$POST_CONTAM" -v total="$PROCESSED" 'BEGIN{printf "%.2f", clean/total*100}')

        INPUT_STAR=$(grep "Number of input reads" "$STAR_LOG" | awk -F'|' '{gsub(/ /,"",$2); print $2}')
        UNIQUE_READS=$(grep "Uniquely mapped reads number" "$STAR_LOG" | awk -F'|' '{gsub(/ /,"",$2); print $2}')
        UNIQUE_PCT=$(grep "Uniquely mapped reads %" "$STAR_LOG" | awk -F'|' '{gsub(/ /,"",$2); print $2}')

        CLEAN_INFO="noN_noPolyG_Bowtie_contaminant_filter"

        echo -e "${SAMPLE}\t${RUN}\t${MODE}\t${RAW_READS}\t${RAW_LEN}\t${RAW_GC}\t${TRIM_READS}\t${REMOVED}\t${TRIM_LEN}\t${TRIM_GC}\t${TRIM_INFO}\t${POST_N_POLYG}\t${POST_CONTAM}\t${RETENTION}\t${CLEAN_GC}\t${CLEAN_LEN}\t${CLEAN_INFO}\t${INPUT_STAR}\t${UNIQUE_READS}\t${UNIQUE_PCT}" >> "$OUTFILE"

    done
done

echo "QC master table saved:"
echo "$OUTFILE"
