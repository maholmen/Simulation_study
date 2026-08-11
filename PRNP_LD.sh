#!/bin/bash

############SLURM arguments#############

#SBATCH --job-name=PRNP_LD
#SBATCH --time=06:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=4
#SBATCH --partition=orion
#SBATCH --output=slurm-%x_%A.out
#SBATCH -e slurm-%x_%A.err
#########################################

###############Main SCRIPT####################
VCF_IN=/mnt/project/CWD_reindeer/maholmen_phd/PRNP/VCF/MERGED_FOR_PRNP_NORWEGIAN_FILTERED.2.CHR11.recode.vcf
SNP_list=/mnt/project/CWD_reindeer/maholmen_phd/PRNP/LD/ORF_SNP_list.txt
NAME=LD_PRNP

cd /mnt/project/CWD_reindeer/maholmen_phd/PRNP/26.02.2026/LD

module load PLINK/1.9b_6.17-x86_64

plink --vcf $VCF_IN --allow-extra-chr --set-missing-var-ids @:# \
    --vcf-half-call missing --ld-snp-list $SNP_list --ld-window 999999 \
    --r2 --ld-window-r2 0.2 --out $NAME
