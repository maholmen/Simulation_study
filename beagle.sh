#!/bin/bash

############SLURM arguments#############

#SBATCH --job-name=beagle
#SBATCH --time=03:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=8
#SBATCH --partition=orion
#SBATCH --output=slurm-%x_%A.out
#SBATCH -e slurm-%x_%A.err
#########################################

###############Main SCRIPT####################

##Variables###
RSYNC='rsync -aPLhv --no-perms --no-owner --no-group'
dir=/mnt/project/CWD_reindeer/maholmen_phd/simulations/VCF
name=MERGED_FOR_PRNP_NORWEGIAN_FILTERED.2
LIST=$1
#number=$3

chr=$(head -n $SLURM_ARRAY_TASK_ID $LIST | tail -n 1)

## activate conda environment
#module load gzip
cd $dir

#gunzip $name.vcf.gz

module purge
echo "Activating Miniconda3 module for $USER"
module load Miniconda3
eval "$(conda shell.bash hook)"
echo "conda is running. Please type conda activate to load the basic conda functions..."
conda activate /mnt/project/CWD_reindeer/condaenvironments/BEAGLE
echo "I'm working with this CONDAENV"
echo $CONDA_PREFIX

export _JAVA_OPTIONS="-Xmx50g"

beagle gt=$name.vcf chrom=$chr nthreads=$SLURM_CPUS_ON_NODE out=$name.$chr.beagle.out
