#!/bin/bash

SAMPLEID='93ff786e-0165-4b02-8d27-806d422e93fc'

REFGENOME='/srv/shared/vanloo/pipeline-files/human/references/alignment/GRCh37-1000G/genome.fa'
NORMALBAM='/srv/shared/vanloo/ICGC/COAD-US/WGS/PCAWG.7e39e79a-7e20-4c7f-87f3-577806f42a04.bam'
TUMORBAM='/srv/shared/vanloo/ICGC/COAD-US/WGS/PCAWG.e72c727c-557e-424d-b045-3562a5a0f9c6.bam'
KNOWNINDELS='/srv/shared/vanloo/pipeline-files/human/references/b37indels_and_snps/Mills_and_1000G_gold_standard.indels.b37.vcf'
KNOWNSNPS='/srv/shared/vanloo/pipeline-files/human/references/b37indels_and_snps/af-only-gnomad.raw.sites.b37.vcf.gz'
OUTBASE=/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/data/pcawg_recall/${SAMPLEID}/

module purge
ml GATK

mkdir ${OUTBASE}

gatk --java-options "-Xmx64G" BaseRecalibrator \
      -I ${NORMALBAM} \
      -R ${REFGENOME} \
      --known-sites ${KNOWNSNPS} \
      -O ${OUTBASE}${SAMPLEID}_normal_recal_data.table

gatk --java-options "-Xmx64G" BaseRecalibrator \
      -I ${TUMORBAM} \
      -R ${REFGENOME} \
      --known-sites ${KNOWNSNPS} \
      -O ${OUTBASE}${SAMPLEID}_tumor_recal_data.table

gatk --java-options "-Xmx64G" ApplyBQSR \
      -R ${REFGENOME} \
      -I ${NORMALBAM} \
      --bqsr-recal-file ${OUTBASE}${SAMPLEID}_normal_recal_data.table \
      -O ${OUTBASE}${SAMPLEID}_normal.aln.recal.bam

gatk --java-options "-Xmx64G" ApplyBQSR \
      -R ${REFGENOME} \
      -I ${TUMORBAM} \
      --bqsr-recal-file ${OUTBASE}${SAMPLEID}_tumor_recal_data.table \
      -O ${OUTBASE}${SAMPLEID}_tumor.aln.recal.bam


echo "Running Mutect2"
gatk --java-options "-Xmx64G" Mutect2 \
      -R ${REFGENOME} \
      -I ${OUTBASE}${SAMPLEID}_tumor.aln.recal.bam \
      -I ${OUTBASE}${SAMPLEID}_normal.aln.recal.bam \
      -tumor 93ff786e-0165-4b02-8d27-806d422e93fc \
      -normal 71807b57-89d9-40ab-bbbb-d63e2c86040c \
      --germline-resource ${KNOWNSNPS} \
      --af-of-alleles-not-in-resource 0.000025 \
      -O ${OUTBASE}${SAMPLEID}_mutect2_snvs_indels.vcf.gz

gatk --java-options "-Xmx64G" FilterMutectCalls \
      -V ${OUTBASE}${SAMPLEID}_mutect2_snvs_indels.vcf.gz \
      --max-alt-allele-count 2 \
      -O ${OUTBASE}${SAMPLEID}_mutect2_snvs_indels_oncefiltered.vcf.gz

gatk --java-options "-Xmx64G" CollectSequencingArtifactMetrics \
      -I ${OUTBASE}${SAMPLEID}_tumor.aln.recal.bam \
      -O ${SAMPLEID}_ \
      -R ${REFGENOME} \
      --DB_SNP ${KNOWNSNPS}

# and FilterByOrientationBias / Oxo-G etc
gatk --java-options "-Xmx64G" FilterByOrientationBias \
      --artifact-modes 'G/T' \
      -V ${OUTBASE}${SAMPLEID}_mutect2_snvs_indels_oncefiltered.vcf.gz \
      -P ${SAMPLEID}_pre_adapter_detail_metrics \
      -O ${OUTBASE}${SAMPLEID}_mutect2_snvs_indels_twicefiltered.vcf.gz
