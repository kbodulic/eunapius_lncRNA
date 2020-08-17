#!/bin/bash

QUERY_SEQUENCE=$1
DATABASE=$2
WINDOW=$3
LEVEL=$4
DUSTMASKER=path/to/dustmasker
PERL_SCRIPT_RENAME=/path/to/rename.pl
MAKEBLASTDB=/path/to/makeblastdb
BLASTN=/path/to/blastn

DUSTMASKER_OUTPUT=${DATABASE%.*}_w{$WINDOW}_l{$LEVEL}.fa 
$DUSTMASKER -in $DATABASE -infmt fasta -outfmt fasta -window $WINDOW -level $LEVEL -out $DUSTMASKER_OUTPUT 
TR_OUTPUT=${DUSTMASKER_OUTPUT%.*}_tr.fa

cat $DUSTMASKER_OUTPUT  | tr "a" "N" | tr "t" "N" | tr "g" "N" | tr "c" "N" | tr "n" "N" > $TR_OUTPUT 

GREP_OUTPUT=${DATABASE%.*}_headers.txt

cat $DATABASE | grep ">" > $GREP_OUTPUT 


PERL_OUTPUT=${TR_OUTPUT%.*}_replaced_headers.fa

perl $PERL_SCRIPT_RENAME $GREP_OUTPUT $TR_OUTPUT > $PERL_OUTPUT 


$MAKEBLASTDB -in $PERL_OUTPUT -input_type fasta -dbtype nucl -title ${TR_OUTPUT%.*}  

QUERY_SEQUENCE_BLAST_VAR=$(basename $QUERY_SEQUENCE)
DATABASE_BLAST_VAR=$(basename $DATABASE)
BLAST_OUTPUT=sponges_blast.outfmt6

$BLASTN -task blastn -query $QUERY_SEQUENCE -db $PERL_OUTPUT -num_threads 10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore" -out $BLAST_OUTPUT 


