#!/bin/bash


QUERY=$1
DATABASE=$2
TRANS=$3

ls | grep $TRANS | parallel cut -f1,2,3 {} ">" {.}.bed 

samtools faidx $QUERY


samtools faidx $DATABASE

ls | grep \.bed | parallel bedtools getfasta -fi $DATABASE -bed {} -fo {.}_pre.fa -name

ls | grep -Po ".*(?=_msa_${3}_pre\.fa)" | parallel samtools faidx $QUERY {} ">>" {}_msa_${3}_pre.fa
