#! /bin/bash

samf="../data/SamplesID.txt"

while read SAMPLE_ID

./ExmFilt.CustomGenotype.py -v JOINTVCF.vcf.gz -o SAMPLE_ID -n SAMPLE_ID -V –S
 
##-n = output any non-reference variant for SAMPLE_ID
##-V = output a vcf file as well as a tsv
##-S = output synonymous variants as well
done < "$samf"
