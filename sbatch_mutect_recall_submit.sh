#!/bin/bash
#
#SBATCH --job-name=isa_m2_pcawg4
#SBATCH --output=isa_m2_pcawg4%A_%a.out
#SBATCH --error=isa_m2_pcawg4%A_%a.err
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#
#SBATCH --array=2-209
#SBATCH --exclude=compute001

declare -a naidar
declare -a taidar
declare -a nbamar
declare -a tbamar
declare -a projar

while read naid taid nbam tbam proj;
do
    naidar+=($naid)
    taidar+=($taid)
    nbamar+=($nbam)
    tbamar+=($tbam)
    projar+=($proj)
done < /srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/data/samples_to_recall.txt

ml GATK

srun /srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/code/recall_pcawg_p1_BQSR_M2normal.sh ${naidar[$SLURM_ARRAY_TASK_ID]} ${taidar[$SLURM_ARRAY_TASK_ID]} ${nbamar[$SLURM_ARRAY_TASK_ID]} ${tbamar[$SLURM_ARRAY_TASK_ID]} ${projar[$SLURM_ARRAY_TASK_ID]}

module purge
