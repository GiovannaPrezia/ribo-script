#!/bin/bash

set -eu

# ======================================================
# CONFIGURAÇÕES — PRJEB29208 / iPSC-CM Day 21
# ======================================================

BASE_DIR="$HOME/Diretório/data/bioprojects/Germany_PRJEB29208"
RESOURCE_DIR="$HOME/Diretório/data/bioprojects/Resources"

SRA_DIR="$BASE_DIR/01_SRAs/ribo_seq"
RAW_DIR="$BASE_DIR/02_fastq/ribo_seq"
FASTQC_RAW_DIR="$RAW_DIR/fastqc_raw"

TRIM_DIR="$BASE_DIR/03_trimmed/ribo_seq"
FASTQC_TRIM_DIR="$TRIM_DIR/fastqc_trimmed"

CLEAN_DIR="$BASE_DIR/04_cleaned/ribo_seq"
FASTQC_CLEAN_DIR="$CLEAN_DIR/fastqc_cleaned"

ALIGN_DIR="$BASE_DIR/05_alignment/ribo_seq"
STAR_QC_DIR="$BASE_DIR/06_star_qc/ribo_seq"
COUNT_DIR="$BASE_DIR/07_counts/ribo_seq"

RIBOTRICER_DIR="$BASE_DIR/10_Ribotricer"
MULTIQC_DIR="$BASE_DIR/11_MultiQC"
LOG_DIR="$BASE_DIR/logs"
QC_DIR="$BASE_DIR/QC_tables"

RNA_DICT="$RESOURCE_DIR/rnas_dictionary/indexes/bowtie1/rnas_dictionary_human"
STAR_INDEX="$RESOURCE_DIR/hg38_star_index"
GTF="$RESOURCE_DIR/gencode.v45.annotation.gtf"
GENOME_FA="$RESOURCE_DIR/GRCh38.primary_assembly.genome.fa"

THREADS=20

SIZE_MODES=("all_lengths" "28_36")

SAMPLES=(
"iPSC_DAY21_Ger_rep1"
"iPSC_DAY21_Ger_rep2"
)

RUNS=(
"ERR3367797"
"ERR3367798"
)

declare -A RUN_MAP
RUN_MAP["ERR3367797"]="iPSC_DAY21_Ger_rep1"
RUN_MAP["ERR3367798"]="iPSC_DAY21_Ger_rep2"

mkdir -p \
"$RAW_DIR" "$FASTQC_RAW_DIR" \
"$TRIM_DIR" "$FASTQC_TRIM_DIR" \
"$CLEAN_DIR" "$FASTQC_CLEAN_DIR" \
"$ALIGN_DIR" "$STAR_QC_DIR" \
"$COUNT_DIR" \
"$RIBOTRICER_DIR" "$MULTIQC_DIR" \
"$LOG_DIR" "$QC_DIR"

[[ -f "$GTF" ]] || { echo "ERRO: GTF não encontrado: $GTF"; exit 1; }
[[ -f "$GENOME_FA" ]] || { echo "ERRO: genoma não encontrado: $GENOME_FA"; exit 1; }
[[ -f "$STAR_INDEX/Genome" ]] || { echo "ERRO: STAR index não encontrado: $STAR_INDEX"; exit 1; }
[[ -f "${RNA_DICT}.1.ebwt" ]] || { echo "ERRO: Bowtie index não encontrado: ${RNA_DICT}.1.ebwt"; exit 1; }

# ======================================================
# MENU
# ======================================================

echo "======================================"
echo "Ribo-seq Pipeline — PRJEB29208 D21"
echo "======================================"
echo ""
echo "1 - Pipeline contínua completa"
echo "2 - Pipeline modular/interativa"
echo ""

read -p "Opção: " MODE

if [[ "$MODE" == "1" ]]; then
    PIPELINE_MODE="continuous"
    MODULE="10"
elif [[ "$MODE" == "2" ]]; then
    PIPELINE_MODE="interactive"

    echo ""
    echo "Módulos disponíveis:"
    echo "0 - Converter SRA para FASTQ"
    echo "1 - FastQC bruto"
    echo "2 - Cutadapt"
    echo "3 - QC pós-trimming"
    echo "4 - Filtro noN/noPolyG + all_lengths/28_36"
    echo "5 - Bowtie contaminantes"
    echo "6 - FastQC dados limpos"
    echo "7 - STAR"
    echo "8 - featureCounts"
    echo "9 - MultiQC + relatório"
    echo "10 - Pipeline completa"
    echo ""

    read -p "Escolha o módulo: " MODULE
else
    echo "Opção inválida."
    exit 1
fi

MASTER_LOG="$LOG_DIR/PRJEB29208_D21_pipeline_master.log"

echo "==================================================" | tee "$MASTER_LOG"
echo "Pipeline iniciada em: $(date)" | tee -a "$MASTER_LOG"
echo "Projeto: PRJEB29208" | tee -a "$MASTER_LOG"
echo "Condição: iPSC-derived cardiomyocytes Day 21" | tee -a "$MASTER_LOG"
echo "Size modes: all_lengths e 28_36" | tee -a "$MASTER_LOG"
echo "Modo: $PIPELINE_MODE" | tee -a "$MASTER_LOG"
echo "Módulo: $MODULE" | tee -a "$MASTER_LOG"
echo "==================================================" | tee -a "$MASTER_LOG"

pause_step () {
    if [[ "$PIPELINE_MODE" == "interactive" ]]; then
        echo ""
        while true; do
            read -p "Continuar (c), modo contínuo (s), pausar para ajuste (r) ou sair (q)? " ANSWER
            case "$ANSWER" in
                c|C) break ;;
                s|S) PIPELINE_MODE="continuous"; echo "Modo contínuo ativado."; break ;;
                r|R) echo "Pipeline pausada para ajuste manual."; exit 0 ;;
                q|Q) echo "Pipeline encerrada."; exit 0 ;;
                *) echo "Resposta inválida. Use c, s, r ou q." ;;
            esac
        done
    fi
}

# ======================================================
# FUNÇÕES
# ======================================================

run_sra_conversion () {
    echo "Convertendo SRAs para FASTQ..." | tee -a "$MASTER_LOG"

    for RUN in "${RUNS[@]}"; do
        SAMPLE="${RUN_MAP[$RUN]}"
        SRA_FILE="$SRA_DIR/$RUN/$RUN.sra"
        OUT_FASTQ="$RAW_DIR/${SAMPLE}.fastq.gz"

        echo "Run: $RUN → $SAMPLE" | tee -a "$MASTER_LOG"

        [[ -f "$SRA_FILE" ]] || { echo "ERRO: SRA não encontrado: $SRA_FILE"; exit 1; }

        if [[ -f "$OUT_FASTQ" ]]; then
            echo "FASTQ já existe, pulando: $OUT_FASTQ" | tee -a "$MASTER_LOG"
            continue
        fi

        fasterq-dump "$SRA_FILE" \
            -O "$RAW_DIR" \
            -e "$THREADS" \
            > "$LOG_DIR/${RUN}.fasterq_dump.log" 2>&1

        if [[ -f "$RAW_DIR/${RUN}.fastq" ]]; then
            mv "$RAW_DIR/${RUN}.fastq" "$RAW_DIR/${SAMPLE}.fastq"
            pigz -p "$THREADS" "$RAW_DIR/${SAMPLE}.fastq"
        else
            echo "ERRO: FASTQ não foi gerado para $RUN" | tee -a "$MASTER_LOG"
            exit 1
        fi

        seqkit stats "$OUT_FASTQ" | tee "$QC_DIR/${SAMPLE}.seqkit_raw_stats.txt"
    done
}

run_fastqc_raw () {
    SAMPLE="$1"
    RAW_FASTQ="$RAW_DIR/${SAMPLE}.fastq.gz"

    echo "FastQC bruto: $SAMPLE" | tee -a "$MASTER_LOG"

    fastqc "$RAW_FASTQ" \
        -o "$FASTQC_RAW_DIR" \
        -t "$THREADS" \
        > "$LOG_DIR/${SAMPLE}.fastqc_raw.log" 2>&1
}

run_cutadapt () {
    SAMPLE="$1"
    RAW_FASTQ="$RAW_DIR/${SAMPLE}.fastq.gz"
    TRIM_FASTQ="$TRIM_DIR/${SAMPLE}.trim.fastq.gz"

    echo "Cutadapt: $SAMPLE" | tee -a "$MASTER_LOG"

    cutadapt \
        -u 3 \
        -a "A{10}" \
        -m 17 \
        -o "$TRIM_FASTQ" \
        "$RAW_FASTQ" \
        > "$LOG_DIR/${SAMPLE}.cutadapt.log"

    tail -n 20 "$LOG_DIR/${SAMPLE}.cutadapt.log"
}

run_qc_trimmed () {
    SAMPLE="$1"
    TRIM_FASTQ="$TRIM_DIR/${SAMPLE}.trim.fastq.gz"

    echo "QC pós-trimming: $SAMPLE" | tee -a "$MASTER_LOG"

    zcat "$TRIM_FASTQ" \
        | head -n 4000000 \
        | awk 'NR%4==2 {print length($0)}' \
        | sort -n \
        | uniq -c \
        > "$QC_DIR/${SAMPLE}.length_distribution_trimmed.txt"

    zcat "$TRIM_FASTQ" \
        | head -n 4000000 \
        | awk 'NR%4==2 {count[$0]++} END {for(seq in count) print count[seq], seq}' \
        | sort -nr \
        | head -n 30 \
        > "$QC_DIR/${SAMPLE}.top_sequences_trimmed.txt"

    cat "$QC_DIR/${SAMPLE}.length_distribution_trimmed.txt"
    cat "$QC_DIR/${SAMPLE}.top_sequences_trimmed.txt"

    fastqc "$TRIM_FASTQ" \
        -o "$FASTQC_TRIM_DIR" \
        -t "$THREADS" \
        > "$LOG_DIR/${SAMPLE}.fastqc_trimmed.log" 2>&1
}

run_filter () {
    SAMPLE="$1"
    TRIM_FASTQ="$TRIM_DIR/${SAMPLE}.trim.fastq.gz"

    FILTERED_ALL="$CLEAN_DIR/${SAMPLE}.trim.noN.noPolyG.all_lengths.fastq.gz"
    FILTERED_RPF="$CLEAN_DIR/${SAMPLE}.trim.noN.noPolyG.28_36.fastq.gz"

    echo "Filtro noN/noPolyG + geração all_lengths e 28_36: $SAMPLE" | tee -a "$MASTER_LOG"

    zcat "$TRIM_FASTQ" \
        | awk 'NR%4==1{h=$0} NR%4==2{s=$0} NR%4==3{p=$0} NR%4==0{q=$0; if(s !~ /N/ && s !~ /^G+$/) print h"\n"s"\n"p"\n"q}' \
        | gzip > "$FILTERED_ALL"

    zcat "$TRIM_FASTQ" \
        | awk 'NR%4==1{h=$0} NR%4==2{s=$0} NR%4==3{p=$0} NR%4==0{q=$0; if(s !~ /N/ && s !~ /^G+$/ && length(s)>=28 && length(s)<=36) print h"\n"s"\n"p"\n"q}' \
        | gzip > "$FILTERED_RPF"

    zcat "$FILTERED_ALL" \
        | head -n 4000000 \
        | awk 'NR%4==2 {print length($0)}' \
        | sort -n \
        | uniq -c \
        > "$QC_DIR/${SAMPLE}.length_distribution_noN_noPolyG_all_lengths.txt"

    zcat "$FILTERED_RPF" \
        | head -n 4000000 \
        | awk 'NR%4==2 {print length($0)}' \
        | sort -n \
        | uniq -c \
        > "$QC_DIR/${SAMPLE}.length_distribution_noN_noPolyG_28_36.txt"

    cat "$QC_DIR/${SAMPLE}.length_distribution_noN_noPolyG_all_lengths.txt"
    cat "$QC_DIR/${SAMPLE}.length_distribution_noN_noPolyG_28_36.txt"
}

run_bowtie () {
    SAMPLE="$1"

    for SIZE_MODE in "${SIZE_MODES[@]}"; do
        INPUT_FASTQ="$CLEAN_DIR/${SAMPLE}.trim.noN.noPolyG.${SIZE_MODE}.fastq.gz"
        MODE_CLEAN_DIR="$CLEAN_DIR/$SIZE_MODE"
        mkdir -p "$MODE_CLEAN_DIR"

        CLEAN_FASTQ="$MODE_CLEAN_DIR/${SAMPLE}.${SIZE_MODE}.clean.fastq"

        echo "Bowtie contaminantes: $SAMPLE - $SIZE_MODE" | tee -a "$MASTER_LOG"

        bowtie \
            -p "$THREADS" \
            -v 1 \
            --un "$CLEAN_FASTQ" \
            "$RNA_DICT" \
            "$INPUT_FASTQ" \
            /dev/null \
            2> "$LOG_DIR/${SAMPLE}.${SIZE_MODE}.bowtie.log"

        gzip -f "$CLEAN_FASTQ"

        cat "$LOG_DIR/${SAMPLE}.${SIZE_MODE}.bowtie.log"
    done
}

run_fastqc_clean () {
    SAMPLE="$1"

    for SIZE_MODE in "${SIZE_MODES[@]}"; do
        CLEAN_FASTQ="$CLEAN_DIR/$SIZE_MODE/${SAMPLE}.${SIZE_MODE}.clean.fastq.gz"

        echo "FastQC clean: $SAMPLE - $SIZE_MODE" | tee -a "$MASTER_LOG"

        fastqc "$CLEAN_FASTQ" \
            -o "$FASTQC_CLEAN_DIR" \
            -t "$THREADS" \
            > "$LOG_DIR/${SAMPLE}.${SIZE_MODE}.fastqc_clean.log" 2>&1
    done
}

run_star () {
    SAMPLE="$1"

    for SIZE_MODE in "${SIZE_MODES[@]}"; do
        MODE_ALIGN_DIR="$ALIGN_DIR/$SIZE_MODE"
        mkdir -p "$MODE_ALIGN_DIR"

        CLEAN_FASTQ="$CLEAN_DIR/$SIZE_MODE/${SAMPLE}.${SIZE_MODE}.clean.fastq.gz"

        echo "STAR: $SAMPLE - $SIZE_MODE" | tee -a "$MASTER_LOG"

        STAR \
            --runThreadN "$THREADS" \
            --genomeDir "$STAR_INDEX" \
            --readFilesIn "$CLEAN_FASTQ" \
            --readFilesCommand zcat \
            --outFileNamePrefix "$MODE_ALIGN_DIR/${SAMPLE}.${SIZE_MODE}_" \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS nM MD \
            --seedSearchStartLmaxOverLread 0.5 \
            --outFilterMismatchNmax 2 \
            --outMultimapperOrder Random \
            --outFilterMultimapNmax 20 \
            --outSAMmultNmax 1

        BAM="$MODE_ALIGN_DIR/${SAMPLE}.${SIZE_MODE}_Aligned.sortedByCoord.out.bam"

        samtools index "$BAM"

        samtools idxstats "$BAM" \
            > "$STAR_QC_DIR/${SAMPLE}.${SIZE_MODE}.idxstats.txt"

        cp "$MODE_ALIGN_DIR/${SAMPLE}.${SIZE_MODE}_Log.final.out" \
           "$STAR_QC_DIR/${SAMPLE}.${SIZE_MODE}_Log.final.out"

        cat "$MODE_ALIGN_DIR/${SAMPLE}.${SIZE_MODE}_Log.final.out"
    done
}

run_featurecounts () {
    SAMPLE="$1"

    for SIZE_MODE in "${SIZE_MODES[@]}"; do
        MODE_ALIGN_DIR="$ALIGN_DIR/$SIZE_MODE"
        MODE_COUNT_DIR="$COUNT_DIR/$SIZE_MODE"
        mkdir -p "$MODE_COUNT_DIR"

        BAM="$MODE_ALIGN_DIR/${SAMPLE}.${SIZE_MODE}_Aligned.sortedByCoord.out.bam"

        echo "featureCounts: $SAMPLE - $SIZE_MODE" | tee -a "$MASTER_LOG"

        featureCounts \
            -T "$THREADS" \
            -t CDS \
            -g gene_id \
            -a "$GTF" \
            -o "$MODE_COUNT_DIR/${SAMPLE}.${SIZE_MODE}.CDS_counts.txt" \
            "$BAM" \
            2> "$LOG_DIR/${SAMPLE}.${SIZE_MODE}.featureCounts.log"

        cat "$MODE_COUNT_DIR/${SAMPLE}.${SIZE_MODE}.CDS_counts.txt.summary"
    done
}

run_report () {
    echo "Gerando MultiQC e relatório final..." | tee -a "$MASTER_LOG"

    multiqc "$BASE_DIR" \
        -o "$MULTIQC_DIR" \
        -n PRJEB29208_iPSC_DAY21_Ger_MultiQC.html \
        -f

    REPORT_TSV="$MULTIQC_DIR/PRJEB29208_iPSC_DAY21_Ger_QC_report.tsv"

    echo -e "Sample\tRun\tMode\tInput_STAR\tUnique_Reads\tUnique_%\tAvg_Mapped_Length\tToo_Many_Loci_%\tToo_Short_%\tCDS_Assigned\tNoFeatures\tAmbiguity" > "$REPORT_TSV"

    for RUN in "${RUNS[@]}"; do
        SAMPLE="${RUN_MAP[$RUN]}"

        for SIZE_MODE in "${SIZE_MODES[@]}"; do
            STAR_LOG="$ALIGN_DIR/$SIZE_MODE/${SAMPLE}.${SIZE_MODE}_Log.final.out"
            FC_SUMMARY="$COUNT_DIR/$SIZE_MODE/${SAMPLE}.${SIZE_MODE}.CDS_counts.txt.summary"

            [[ -f "$STAR_LOG" && -f "$FC_SUMMARY" ]] || continue

            INPUT_STAR=$(grep "Number of input reads" "$STAR_LOG" | awk -F'|' '{gsub(/ /,"",$2); print $2}')
            UNIQUE_READS=$(grep "Uniquely mapped reads number" "$STAR_LOG" | awk -F'|' '{gsub(/ /,"",$2); print $2}')
            UNIQUE_PCT=$(grep "Uniquely mapped reads %" "$STAR_LOG" | awk -F'|' '{gsub(/ /,"",$2); print $2}')
            AVG_LEN=$(grep "Average mapped length" "$STAR_LOG" | awk -F'|' '{gsub(/ /,"",$2); print $2}')
            TOO_MANY=$(grep "% of reads mapped to too many loci" "$STAR_LOG" | awk -F'|' '{gsub(/ /,"",$2); print $2}')
            TOO_SHORT=$(grep "% of reads unmapped: too short" "$STAR_LOG" | awk -F'|' '{gsub(/ /,"",$2); print $2}')

            ASSIGNED=$(grep "^Assigned" "$FC_SUMMARY" | awk '{print $2}')
            NOFEATURES=$(grep "^Unassigned_NoFeatures" "$FC_SUMMARY" | awk '{print $2}')
            AMBIGUITY=$(grep "^Unassigned_Ambiguity" "$FC_SUMMARY" | awk '{print $2}')

            echo -e "${SAMPLE}\t${RUN}\t${SIZE_MODE}\t${INPUT_STAR}\t${UNIQUE_READS}\t${UNIQUE_PCT}\t${AVG_LEN}\t${TOO_MANY}\t${TOO_SHORT}\t${ASSIGNED}\t${NOFEATURES}\t${AMBIGUITY}" >> "$REPORT_TSV"
        done
    done

    echo "Relatório final gerado em: $MULTIQC_DIR"
}

# ======================================================
# EXECUÇÃO
# ======================================================

if [[ "$MODULE" == "0" || "$MODULE" == "10" ]]; then
    run_sra_conversion
    pause_step
fi

for SAMPLE in "${SAMPLES[@]}"; do
    echo "======================================"
    echo "Amostra: $SAMPLE"
    echo "======================================" | tee -a "$MASTER_LOG"

    case "$MODULE" in
        1) run_fastqc_raw "$SAMPLE"; pause_step ;;
        2) run_cutadapt "$SAMPLE"; pause_step ;;
        3) run_qc_trimmed "$SAMPLE"; pause_step ;;
        4) run_filter "$SAMPLE"; pause_step ;;
        5) run_bowtie "$SAMPLE"; pause_step ;;
        6) run_fastqc_clean "$SAMPLE"; pause_step ;;
        7) run_star "$SAMPLE"; pause_step ;;
        8) run_featurecounts "$SAMPLE"; pause_step ;;
        10)
            run_fastqc_raw "$SAMPLE"; pause_step
            run_cutadapt "$SAMPLE"; pause_step
            run_qc_trimmed "$SAMPLE"; pause_step
            run_filter "$SAMPLE"; pause_step
            run_bowtie "$SAMPLE"; pause_step
            run_fastqc_clean "$SAMPLE"; pause_step
            run_star "$SAMPLE"; pause_step
            run_featurecounts "$SAMPLE"; pause_step
            ;;
        0|9) ;;
        *) echo "Módulo inválido."; exit 1 ;;
    esac
done

if [[ "$MODULE" == "9" || "$MODULE" == "10" ]]; then
    run_report
fi

echo "Finalizado em: $(date)" | tee -a "$MASTER_LOG"
