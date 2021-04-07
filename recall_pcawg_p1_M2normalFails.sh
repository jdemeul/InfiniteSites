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

# variables
NORMALBAM="/srv/shared/vanloo/ICGC/${PROJ}/WGS/${NBAM}"
TUMORBAM="/srv/shared/vanloo/ICGC/${PROJ}/WGS/${TBAM}"

OUTBASE="/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/data/pcawg_recall/${PROJ}/${TAID}/"


echo "Running sample ${TAID}, project ${PROJ}"
mkdir -p ${OUTBASE}

if [ ! -f "${OUTBASE}${TAID}_normal4PoN_mutect2_snvs_indels.vcf.gz" ]
then
echo "Running Mutect2 tumour-only mode PoN"
gatk --java-options "-Xmx6G" Mutect2 \
      -R ${REFGENOME} \
      -I ${OUTBASE}${TAID}_normal.aln.recal.bam \
      -tumor ${NAID} \
      --germline-resource ${KNOWNSNPS} \
      --max-population-af 0.01 \
      -O ${OUTBASE}${TAID}_normal4PoN_mutect2_snvs_indels.vcf.gz
fi

echo "Finished running sample ${TAID}, project ${PROJ}"
