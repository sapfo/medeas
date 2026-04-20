#!/bin/bash

#Defining input file
vcffile="CEU_YRI.chr22.vcf.gz"
label="label.lab"

#Defining output file
output_dir="result"
output_log=${output_dir}"/result.log"


#Name for temp file
plink_dir="plink"
plink_name="all_individual"
label_haploid="label_haploid.lab"
med_dir="med"
med_sorted="all_haplotype.med"
med_random="all_haplotype_randomized.med"

#From diploid to haploid label
awk '{print $1}{print $1}' $label > $label_haploid

#Using plink to get a file with 1 and 2
mkdir $plink_dir
plink --vcf  $vcffile  --allow-extra-chr --recode12 --transpose --out ${plink_dir}/${plink_name}


# Removing the first useless column
mkdir $med_dir
cut -d" " -f 5- ${plink_dir}/${plink_name}.tped > ${med_dir}/${med_sorted}

# Mixing the two haplotypes of each individuals
awk 'BEGIN{srand()}{for (i=0;i<NF/2;i++) if (rand() > 0.5) printf "%s %s ", $(2*i+1) ,$(2*i+2); else printf "%s %s ", $(2*i+2), $(2*i+1);print ""}' ${med_dir}/${med_sorted} > ${med_dir}/${med_random}

#Running the actual software
mkdir $output_dir
python ../main.py -sf ${med_dir}/${med_random} -lf $label_haploid -of $output_dir -bws 1000 > $output_log
