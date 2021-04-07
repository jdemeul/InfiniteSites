#!/bin/bash
#
#SBATCH --job-name=isa_minibams
#SBATCH --output=isa_minibams%A_%a.out
#SBATCH --error=isa_minibams%A_%a.err
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=10G
#
#SBATCH --array=2-209
#SBATCH --exclude=fat-worker00[1-3]

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

ml foss

srun /srv/shared/vanloo/home/jdemeul/projects/2016-17_ICGC/infinite_sites/code/recall_pcawg_p4_minibams.sh ${naidar[$SLURM_ARRAY_TASK_ID]} ${taidar[$SLURM_ARRAY_TASK_ID]} ${nbamar[$SLURM_ARRAY_TASK_ID]} ${tbamar[$SLURM_ARRAY_TASK_ID]} ${projar[$SLURM_ARRAY_TASK_ID]}

module purge
