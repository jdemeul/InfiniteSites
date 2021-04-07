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
VAR="/srv/shared/vanloo/home/jdemeul/sw/VariantBam/src/variant"

# variables
NORMALBAM="/srv/shared/vanloo/ICGC/${PROJ}/WGS/${NBAM}"
TUMORBAM="/srv/shared/vanloo/ICGC/${PROJ}/WGS/${TBAM}"

OUTBASE="/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/data/pcawg_recall/${PROJ}/${TAID}/"


echo "Running sample ${TAID}, project ${PROJ}"


if [ ! -f "${OUTBASE}${TAID}_normal.aln.recal.mini.bam" ]
then
      echo "Generating mini normal bam"
      ${VAR} ${OUTBASE}${TAID}_normal.aln.recal.bam \
            -l ${OUTBASE}${TAID}_tumor_mutect2_snvs_indels.vcf.gz \
            -o ${OUTBASE}${TAID}_normal.aln.recal.mini.bam \
            -b \
            -v \
            -t 3
fi


if [ ! -f "${OUTBASE}${TAID}_tumor.aln.recal.mini.bam" ]
then
      echo "Generating mini tumour bam"
      ${VAR} ${OUTBASE}${TAID}_tumor.aln.recal.bam \
            -l ${OUTBASE}${TAID}_tumor_mutect2_snvs_indels.vcf.gz \
            -o ${OUTBASE}${TAID}_tumor.aln.recal.mini.bam \
            -b \
            -v \
            -t 3
fi

echo "Finished running sample ${TAID}, project ${PROJ}"
