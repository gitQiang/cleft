#! /bin/bash

samf="../data/SamplesID.txt"
JOINTVCF="../data/CLCP_June2015.vcf.gz"
OUTPUT="../variantcalling/"

while read SAMPLE_ID
do
echo ./ExmFilt.CustomGenotype.py -v $JOINTVCF -o $OUTPUT$SAMPLE_ID -n $SAMPLE_ID -V –S
./ExmFilt.CustomGenotype.py -v $JOINTVCF -o $OUTPUT$SAMPLE_ID -n $SAMPLE_ID -V –S
 
##-n = output any non-reference variant for SAMPLE_ID
##-V = output a vcf file as well as a tsv
##-S = output synonymous variants as well
done < "$samf"
