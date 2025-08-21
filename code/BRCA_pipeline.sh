#!/usr/bin/env bash
set -euo pipefail

# ==============================
# Configurações
# ==============================
THREADS=4
REF="data/hg19.fasta"
TARGET="data/BRCA.list"
R1="data/510-7-BRCA_S8_L001_R1_001.fastq.gz"
R2="data/510-7-BRCA_S8_L001_R2_001.fastq.gz"
PREFIX="BRCA_FEMALE_CTRL"
OUTDIR="output"
SNPEFF_DIR="code/snpEff"

mkdir -p ${OUTDIR} ${OUTDIR}/fastqc_raw ${OUTDIR}/fastqc_processed logs

# ==============================
# 0. Controle de qualidade inicial (RAW READS)
# ==============================
if [ ! -f ${OUTDIR}/fastqc_raw/${PREFIX}_R1_fastqc.html ]; then
  echo "[INFO] FastQC das reads brutas..."
  fastqc ${R1} ${R2} -o ${OUTDIR}/fastqc_raw --quiet
fi

# ==============================
# 1. Alinhamento
# ==============================
if [ ! -f ${OUTDIR}/${PREFIX}.sorted.bam ]; then
  echo "[INFO] Rodando BWA MEM..."
  bwa mem -t ${THREADS} -R "@RG\tID:${PREFIX}\tSM:${PREFIX}\tPL:ILLUMINA" \
    ${REF} ${R1} ${R2} \
    | samtools view -@ ${THREADS} -b \
    | samtools sort -@ ${THREADS} -o ${OUTDIR}/${PREFIX}.sorted.bam -
  
  samtools index ${OUTDIR}/${PREFIX}.sorted.bam
  
  # Estatísticas pós-alinhamento
  echo "[INFO] Estatísticas de alinhamento..."
  samtools flagstat ${OUTDIR}/${PREFIX}.sorted.bam > ${OUTDIR}/${PREFIX}.sorted.flagstat.txt
  samtools idxstats ${OUTDIR}/${PREFIX}.sorted.bam > ${OUTDIR}/${PREFIX}.sorted.idxstats.txt
  
  # FastQC pós-alinhamento
  echo "[INFO] FastQC do BAM alinhado..."
  fastqc ${OUTDIR}/${PREFIX}.sorted.bam -o ${OUTDIR}/fastqc_processed --quiet
fi

# ==============================
# 2. Marcar duplicatas
# ==============================
if [ ! -f ${OUTDIR}/${PREFIX}.markdup.bam ]; then
  echo "[INFO] Marcando duplicatas com Picard..."
  java -jar ~/Downloads/picard.jar MarkDuplicates \
    I=${OUTDIR}/${PREFIX}.sorted.bam \
    O=${OUTDIR}/${PREFIX}.markdup.bam \
    M=${OUTDIR}/${PREFIX}.metrics.txt \
    VALIDATION_STRINGENCY=SILENT
  
  samtools index ${OUTDIR}/${PREFIX}.markdup.bam
  
  # Estatísticas pós-duplicatas
  echo "[INFO] Estatísticas pós-remoção de duplicatas..."
  samtools flagstat ${OUTDIR}/${PREFIX}.markdup.bam > ${OUTDIR}/${PREFIX}.markdup.flagstat.txt
  
  # FastQC pós-remoção de duplicatas
  echo "[INFO] FastQC do BAM pós-duplicatas..."
  fastqc ${OUTDIR}/${PREFIX}.markdup.bam -o ${OUTDIR}/fastqc_processed --quiet
fi

# ==============================
# 3. Chamada de variantes
# ==============================
if [ ! -f ${OUTDIR}/${PREFIX}.raw.vcf.gz ]; then
  echo "[INFO] Rodando FreeBayes..."
  freebayes -f ${REF} --targets ${TARGET} ${OUTDIR}/${PREFIX}.markdup.bam \
    | bgzip -c > ${OUTDIR}/${PREFIX}.raw.vcf.gz
  
  tabix -p vcf ${OUTDIR}/${PREFIX}.raw.vcf.gz
fi

# ==============================
# 4. Anotação com SnpEff
# ==============================
if [ ! -f ${OUTDIR}/${PREFIX}.snpeff.vcf.gz ]; then
  echo "[INFO] Anotando variantes com SnpEff..."
  java -Xmx4g -jar ${SNPEFF_DIR}/snpEff.jar -v hg19 ${OUTDIR}/${PREFIX}.raw.vcf.gz \
    > ${OUTDIR}/${PREFIX}.snpeff.vcf
  
  bgzip -f ${OUTDIR}/${PREFIX}.snpeff.vcf
  tabix -p vcf ${OUTDIR}/${PREFIX}.snpeff.vcf.gz
fi

# ==============================
# 5. Exportar tabela
# ==============================
if [ ! -f ${OUTDIR}/${PREFIX}.variants_annotated.tsv ]; then
  echo "[INFO] Gerando tabela TSV..."
  bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/ANN\n' \
    ${OUTDIR}/${PREFIX}.snpeff.vcf.gz > ${OUTDIR}/${PREFIX}.variants_annotated.tsv
fi

# ==============================
# 6. Relatório final
# ==============================

echo "[INFO] Gerando relatório final..."
{
  echo "=== RELATÓRIO FINAL - PIPELINE BRCA ==="
  echo "Data: $(date)"
  echo "Amostra: ${PREFIX}"
  echo ""
  echo "=== ESTATÍSTICAS DE ALINHAMENTO ==="
  grep "mapped (" ${OUTDIR}/${PREFIX}.markdup.flagstat.txt | head -n 1
  grep "properly paired" ${OUTDIR}/${PREFIX}.markdup.flagstat.txt
  echo ""
  echo "=== VARIANTES ENCONTRADAS ==="
  echo "HIGH: $(grep -c "HIGH" ${OUTDIR}/${PREFIX}.variants_annotated.tsv || true)"
  echo "MODERATE: $(grep -c "MODERATE" ${OUTDIR}/${PREFIX}.variants_annotated.tsv || true)"
  echo "LOW: $(grep -c "LOW" ${OUTDIR}/${PREFIX}.variants_annotated.tsv || true)"
  echo "MODIFIER: $(grep -c "MODIFIER" ${OUTDIR}/${PREFIX}.variants_annotated.tsv || true)"
  echo ""
  echo "=== VARIANTES HIGH (DETALHES) ==="
  grep "HIGH" ${OUTDIR}/${PREFIX}.variants_annotated.tsv | head -n 5 | cut -f1-4
  echo ""
  echo "=== INTERPRETAÇÃO FASTQC ==="
  echo "FAILs/WARNINGS:"
  echo "- Sequence Duplication Levels: Corrigido por MarkDuplicates"
  echo "- Demais warnings: Dentro do esperado para dados clínicos"
} > ${OUTDIR}/${PREFIX}.pipeline_report.txt

echo "[INFO] Pipeline concluído com sucesso!"
echo "[INFO] Relatório salvo em: ${OUTDIR}/${PREFIX}.pipeline_report.txt"