# 🧬 RiboScript Pipeline: Detecção de smORFs em lncRNAs

Este repositório contém um pipeline automatizado em Bash para processar dados de Ribo-seq, com foco na detecção de smORFs presentes em lncRNAs. O objetivo é gerar arquivos alinhados e quantificados prontos para downstream de análise de tradução e expressão.




## 🔁 Etapas do Workflow

O pipeline executa as seguintes etapas:

1. Download dos dados SRA com `prefetch`
2. Conversão para FASTQ com `fasterq-dump`
3. Controle de qualidade com `FastQC`
4. Remoção de adaptadores com `Cutadapt`
5. Filtragem de contaminantes com `Bowtie1`
6. Alinhamento com `STAR`
7. Quantificação com `featureCounts`
8. QC do STAR
9. Relatório integrativo com `MultiQC`

## 🧪 Requisitos

- `sra-tools`
- `fastqc`
- `cutadapt`
- `bowtie`
- `STAR`
- `samtools`
- `subread` (featureCounts)
- `multiqc`

## ⚙️ Como usar

1. Clone o repositório:
```bash
git clone https://github.com/seu-usuario/ribo-seq-smorfs-lncrnas.git
cd ribo-seq-smorfs-lncrnas
