#!/bin/bash

# arg parsing
NAID=$1
TAID=$2
NBAM=$3
TBAM=$4
PROJ=$5

# statics
REFGENOME='/srv/shared/vanloo/pipeline-files/human/references/alignment/GRCh37-1000G/genome.fa'
KNOWNINDELS='/srv/shared/vanloo/pipeline-files/human/references/b37indels_and_snps/Mills_and_1000G_gold_standard.indels.b37.vcf'
KNOWNSNPS='/srv/shared/vanloo/pipeline-files/human/references/b37indels_and_snps/af-only-gnomad.raw.sites.b37.vcf.gz'
PON="/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/data/pcawg_recall/ISA_PCAWG_Po205N.vcf.gz"

# variables
NORMALBAM="/srv/shared/vanloo/ICGC/${PROJ}/WGS/${NBAM}"
TUMORBAM="/srv/shared/vanloo/ICGC/${PROJ}/WGS/${TBAM}"

OUTBASE="/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/data/pcawg_recall/${PROJ}/${TAID}/"


echo "Running sample ${TAID}, project ${PROJ}"

if [ ! -f "${OUTBASE}${TAID}_tumor_mutect2_snvs_indels.vcf.gz" ]
then
      echo "Running Mutect2 on tumour"
      gatk --java-options "-Xmx16G" Mutect2 \
            -R ${REFGENOME} \
            -I ${OUTBASE}${TAID}_tumor.aln.recal.bam \
            -I ${OUTBASE}${TAID}_normal.aln.recal.bam \
            -tumor ${TAID} \
            -normal ${NAID} \
            -pon ${PON} \
            --germline-resource ${KNOWNSNPS} \
            --af-of-alleles-not-in-resource 0.000025 \
            -O ${OUTBASE}${TAID}_tumor_mutect2_snvs_indels.vcf.gz
fi


if [ ! -f "${OUTBASE}${TAID}_tumor_mutect2_snvs_indels_oncefiltered.vcf.gz" ]
then
      echo "Filtering Mutect2 tumour calls allowing biallelics"
      gatk --java-options "-Xmx16G" FilterMutectCalls \
            -V ${OUTBASE}${TAID}_tumor_mutect2_snvs_indels.vcf.gz \
            --max-alt-allele-count 2 \
            -O ${OUTBASE}${TAID}_tumor_mutect2_snvs_indels_oncefiltered.vcf.gz
fi


if [ ! -f "${OUTBASE}${TAID}.pre_adapter_detail_metrics" ]
then
      echo "Collecting additional Sequencing artefact metrics"
      gatk --java-options "-Xmx16G" CollectSequencingArtifactMetrics \
            -I ${OUTBASE}${TAID}_tumor.aln.recal.bam \
            -O ${OUTBASE}${TAID} \
            -R ${REFGENOME} \
            --DB_SNP ${KNOWNSNPS}
fi


if [ ! -f "${OUTBASE}${TAID}_tumor_mutect2_snvs_indels_twicefiltered.vcf.gz" ]
then
      echo "Filtering Mutect2 tumour calls second time for orientation biases"
      gatk --java-options "-Xmx16G" FilterByOrientationBias \
            --artifact-modes 'G/T' \
            -V ${OUTBASE}${TAID}_tumor_mutect2_snvs_indels_oncefiltered.vcf.gz \
            -P ${OUTBASE}${TAID}.pre_adapter_detail_metrics \
            -O ${OUTBASE}${TAID}_tumor_mutect2_snvs_indels_twicefiltered.vcf.gz
fi

echo "Finished running sample ${TAID}, project ${PROJ}"
