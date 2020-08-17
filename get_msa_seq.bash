#!/bin/bash


QUERY=$1
DATABASE=$2
TRANS=$3

PARALLEL=/path/to/parallel
BEDTOOLS=/path/to/bedtools
SAMTOOLS=/path/to/samtools

ls | grep $TRANS | $PARALLEL cut -f1,2,3 {} ">" {.}.bed 

samtools faidx $QUERY


samtools faidx $DATABASE

ls | grep \.bed | $PARALLEL $BEDTOOLS getfasta -fi $DATABASE -bed {} -fo {.}_pre.fa -name

ls | grep -Po ".*(?=_msa_${3}_pre\.fa)" | $PARALLEL $SAMTOOLS faidx $QUERY {} ">>" {}_msa_${3}_pre.fa
