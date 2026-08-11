#!/bin/bash
############ SLURM arguments #############
#SBATCH --job-name=PRNP_matrix
#SBATCH --time=01:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=4
#SBATCH --partition=orion
#SBATCH --output=slurm-%x_%A.out
#SBATCH -e slurm-%x_%A.err
#########################################

set -euo pipefail

#################### Variables ##################
LIST=$1
chr=$(head -n $SLURM_ARRAY_TASK_ID $LIST | tail -n 1)

VCF=/mnt/project/CWD_reindeer/maholmen_phd/simulations/VCF/MERGED_FOR_PRNP_NORWEGIAN_FILTERED.2.$chr.beagle.out.vcf.gz
CHR_NUM=$(echo "$chr" | awk -F'_' '{print $2}')
DIR=/mnt/project/CWD_reindeer/maholmen_phd/simulations/selection/matrixes/CHR${CHR_NUM}
NAME=haplotype.$CHR_NUM

cd $DIR
#################### Extract haplotypes ####################
echo "Processing haplotypes..."

module load gzip
gunzip -c $VCF > $NAME.vcf

module load BCFtools

# Extract genotypes and sample names
bcftools query -f '%CHROM\t%POS[\t%GT]\n' $NAME.vcf > all_genotypes.txt
bcftools query -l $NAME.vcf > samples.txt

# Separate haplotypes
awk '{
    hap1 = $1 "\t" $2;
    hap2 = $1 "\t" $2;
    for (i = 3; i <= NF; i++) {
        split($i, gt, "|");
        hap1 = hap1 "\t" (length(gt[1]) > 0 ? gt[1] : ".");
        hap2 = hap2 "\t" (length(gt[2]) > 0 ? gt[2] : ".");
    }
    print hap1 >> "haplotype1_body.txt";
    print hap2 >> "haplotype2_body.txt";
}' all_genotypes.txt

# Add headers
echo -e "CHROM\tPOS\t$(paste -sd '\t' samples.txt)" > header.txt
cat header.txt haplotype1_body.txt > haplotype1.txt
cat header.txt haplotype2_body.txt > haplotype2.txt

# Remove duplicate positions
awk 'NR==1 {print; next} !seen[$2]++ {print}' haplotype1.txt > haplotype1.unique.txt
awk 'NR==1 {print; next} !seen[$2]++ {print}' haplotype2.txt > haplotype2.unique.txt

mv haplotype1.unique.txt haplotype1.txt
mv haplotype2.unique.txt haplotype2.txt

############################################################
# MATRIX CONSTRUCTION
############################################################

# Rename hap1 headers
(head -n1 haplotype1.txt | awk 'BEGIN{OFS="\t"} {for (i=3; i<=NF; i++) $i=$i"_hap1"; print}' > hap1.renamed.txt; tail -n +2 haplotype1.txt >> hap1.renamed.txt)

# Rename hap2 headers
(head -n1 haplotype2.txt | awk 'BEGIN{OFS="\t"} {for (i=3; i<=NF; i++) $i=$i"_hap2"; print}' > hap2.renamed.txt; tail -n +2 haplotype2.txt >> hap2.renamed.txt)

# Sort both
sort -k1,1 -k2,2 hap1.renamed.txt > hap1_sorted.txt
sort -k1,1 -k2,2 hap2.renamed.txt > hap2_sorted.txt

# Merge haplotypes by position
awk -F'\t' '
FNR==1 && NR==FNR {
    for (i=3; i<=NF; i++) hap1_header[i] = $i
    next
}
FNR!=1 && NR==FNR {
    key = $1 FS $2
    hap1[key] = $0
    hap1_cols = NF
    next
}
FNR==1 && NR!=FNR {
    for (i=3; i<=NF; i++) hap2_header[i] = $i
    print "CHROM\tPOS", (length(hap1_header) ? "\t" : "") join(hap1_header), (length(hap2_header) ? "\t" : "") join(hap2_header)
    next
}
{
    key = $1 FS $2
    seen[key] = 1
    hap2[key] = $0
    hap2_cols = NF
}
END {
    PROCINFO["sorted_in"] = "@ind_str_asc"
    for (key in seen) {
        split(key, fields, FS)
        out = fields[1] FS fields[2]

        if (key in hap1) {
            split(hap1[key], f1, FS)
            for (i=3; i<=hap1_cols; i++) out = out FS f1[i]
        } else {
            for (i=3; i<=hap1_cols; i++) out = out FS "."
        }

        if (key in hap2) {
            split(hap2[key], f2, FS)
            for (i=3; i<=hap2_cols; i++) out = out FS f2[i]
        } else {
            for (i=3; i<=hap2_cols; i++) out = out FS "."
        }

        print out
    }
}

function join(arr,  i, s) {
    s = arr[3]
    for (i=4; i in arr; i++) s = s FS arr[i]
    return s
}
' hap1_sorted.txt hap2_sorted.txt > merged_haplotype_genotypes.txt

# Transpose into matrix (haplotypes as rows)
awk '
NR==1 {
    for(i=3;i<=NF;i++) id[i]=$i
    next
}
{
    pos[NR]=$2
    for(i=3;i<=NF;i++){
        hap[i,NR]=$i
    }
    max_row=NR
}
END {
    printf "ID"
    for(r=2;r<=max_row;r++){
        printf "\t" pos[r]
    }
    printf "\n"

    for(i=3;i<=NF;i++){
        printf id[i]
        for(r=2;r<=max_row;r++){
            printf "\t" hap[i,r]
        }
        printf "\n"
    }
}
' merged_haplotype_genotypes.txt > haplo_matrix.$CHR_NUM.txt

echo "✅ Matrix created: haplo_matrix.$CHR_NUM.txt"
