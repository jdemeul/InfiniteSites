#!/bin/bash
#
#SBATCH --job-name=isa_m2_callssa
#SBATCH --output=isa_m2_callssa%A_%a.out
#SBATCH --error=isa_m2_callssa%A_%a.err
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#
#SBATCH --array=1
##SBATCH --exclude=thin-worker0[0-1][0-9]
##SBATCH --array=2-209

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
srun /srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/code/recall_pcawg_p3_M2tumours.sh ${naidar[$SLURM_ARRAY_TASK_ID]} ${taidar[$SLURM_ARRAY_TASK_ID]} ${nbamar[$SLURM_ARRAY_TASK_ID]} ${tbamar[$SLURM_ARRAY_TASK_ID]} ${projar[$SLURM_ARRAY_TASK_ID]}
srun /srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/code/recall_pcawg_p4_minibams.sh ${naidar[$SLURM_ARRAY_TASK_ID]} ${taidar[$SLURM_ARRAY_TASK_ID]} ${nbamar[$SLURM_ARRAY_TASK_ID]} ${tbamar[$SLURM_ARRAY_TASK_ID]} ${projar[$SLURM_ARRAY_TASK_ID]}

module purge
