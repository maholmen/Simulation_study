#!/bin/bash

############ SLURM arguments #############
#SBATCH --job-name=PRNP_haplotypes
#SBATCH --time=01:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=4
#SBATCH --partition=orion
#SBATCH --output=slurm-%x_%A.out
#SBATCH -e slurm-%x_%A.err
#########################################

set -euo pipefail

#################### Variables ####################
DIR=/mnt/project/CWD_reindeer/maholmen_phd/PRNP/26.02.2026/chr11_phased
VCF=/mnt/project/CWD_reindeer/maholmen_phd/PRNP/26.02.2026/BEAGLE_PHASED/MERGED_FOR_PRNP_NORWEGIAN_FILTERED.2.CHR11_phased.beagle.out.vcf.gz
NAME=genotyping
ALLELES=("E" "A" "B" "D")

cd $DIR

#################### Haplotype processing (ONLY ONCE) ####################
echo "Processing haplotypes and labeling alleles..."
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

# making sure there is no duplication
awk 'NR==1 {print; next} !seen[$2]++ {print}' haplotype1.txt > haplotype1.unique.txt
awk 'NR==1 {print; next} !seen[$2]++ {print}' haplotype2.txt > haplotype2.unique.txt

mv haplotype1.unique.txt haplotype1.txt
mv haplotype2.unique.txt haplotype2.txt

############################## GENOTYPING OF PRNP ##################################################

# Subset PRNP region
awk 'NR==1 || !seen[$2]++ && ($2 >= 44485593 && $2 <= 44486363)' haplotype1.txt > haplotype1_PRNP.txt
awk 'NR==1 || !seen[$2]++ && ($2 >= 44485593 && $2 <= 44486363)' haplotype2.txt > haplotype2_PRNP.txt

# Collect haplotype strings
for HAP in 1 2; do
    awk '
    NR==1 {for (i=3;i<=NF;i++) {samples[i]=$i; hap_str[i]=""}; next}
    {for (i=3;i<=NF;i++) hap_str[i]=hap_str[i] $i}
    END {print "Sample\tHaplotype"; for(i=3;i in samples;i++) print samples[i] "\t" hap_str[i]}
    ' haplotype${HAP}_PRNP.txt > haplotype${HAP}_PRNP_strings.txt
done

# Label haplotypes with alleles
for HAP in 1 2; do
    awk -v hap="$HAP" '
    BEGIN {map["00000"]="B"; map["00001"]="A"; map["00011"]="D"; map["11101"]="E"; map["100101"]="E2"}
    NR==1 {print "Sample\tHaplotype\tLabel"; next}
    {label=map[$2]; if(label=="") label="?"; print $1"_hap"hap "\t" $2 "\t" label}' haplotype${HAP}_PRNP_strings.txt \
    > haplotype${HAP}_PRNP_labeled.txt
done



---------------------------
#making a matrix of all variants
# Rename hap1 headers
(head -n1 haplotype1.txt | awk 'BEGIN{OFS="\t"} {for (i=3; i<=NF; i++) $i=$i"_hap1"; print}' > hap1.renamed.txt; tail -n +2 haplotype1.txt >> hap1.renamed.txt)

# Rename hap2 headers
(head -n1 haplotype2.txt | awk 'BEGIN{OFS="\t"} {for (i=3; i<=NF; i++) $i=$i"_hap2"; print}' > hap2.renamed.txt; tail -n +2 haplotype2.txt >> hap2.renamed.txt)

sort -k1,1 -k2,2 hap1.renamed.txt > hap1_sorted.txt
sort -k1,1 -k2,2 hap2.renamed.txt > hap2_sorted.txt


# merge by position
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
        chrom_pos = key
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
