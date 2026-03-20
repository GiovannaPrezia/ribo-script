#!/bin/bash

# ======================================
# Configurações
# ======================================

BASE_DIR="/home/giovanna.prezia/Diretório/data/bioprojects/Columbia_Data"

GTF="/home/giovanna.prezia/Diretório/data/bioprojects/genome/hg38.ensGene.gtf"

THREADS=20

BAM_DIR="$BASE_DIR/05_alignment/ribo_seq"

OUTPUT_DIR="$BASE_DIR/07_counts/ribo_seq"
OUTPUT_FILE="$OUTPUT_DIR/riboseq_counts.txt"

# ======================================
# Criar diretórios
# ======================================

mkdir -p "$OUTPUT_DIR"

# ======================================
# Buscar BAMs
# ======================================

echo "🔹 Procurando arquivos BAM..."

mapfile -t BAMS < <(find "$BAM_DIR" -type f -name "*_Aligned.sortedByCoord.out.bam")

echo "🔹 Arquivos encontrados: ${#BAMS[@]}"

if [ ${#BAMS[@]} -eq 0 ]; then
    echo "❌ Nenhum BAM encontrado em $BAM_DIR"
    exit 1
fi

echo "🔹 Lista de BAMs:"
printf "%s\n" "${BAMS[@]}"

# ======================================
# Rodar featureCounts
# ======================================

echo "🔹 Rodando featureCounts para Ribo-seq..."

featureCounts \
-T "$THREADS" \
-t CDS \
-g gene_id \
-a "$GTF" \
-o "$OUTPUT_FILE" \
"${BAMS[@]}"

# ======================================
# Verificar se rodou corretamente
# ======================================

if [ ! -f "$OUTPUT_FILE" ]; then
    echo "❌ featureCounts falhou — arquivo de saída não foi criado."
    exit 1
fi

# ======================================
# Converter para matriz limpa
# ======================================

echo "🔹 Convertendo saída para CSV..."

tail -n +2 "$OUTPUT_FILE" | cut -f1,7- > "$OUTPUT_DIR/riboseq_counts.csv"

echo "✅ Ribo-seq counts gerados com sucesso!"
echo "📁 Arquivo final:"
echo "$OUTPUT_DIR/riboseq_counts.csv"
