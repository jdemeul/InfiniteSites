#!/bin/bash

# statics
PON="/srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/data/pcawg_recall/ISA_PCAWG_PoN.vcf.gz"

echo "Generating PoN"
find /srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/data/pcawg_recall/ -name "*_normal4PoN_mutect2_snvs_indels.vcf.gz" > /srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/data/pcawg_recall/pon.samplelist.args

gatk --java-options "-Xmx6G" CreateSomaticPanelOfNormals \
      --vcfs /srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/data/pcawg_recall/pon.samplelist.args \
      -O ${PON}

echo "Done generating PoN"
