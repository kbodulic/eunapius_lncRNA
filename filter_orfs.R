args <- commandArgs(trailingOnly = TRUE)  
library(data.table)
library(stringr)
setwd(".")

#ORF filtering - remove transcripts which have no ORFS and transcripts whose longest ORF is shorter than 150 nucleotides (50 amino acids)
#Arguments: 1- rRNA and length filtered transcript table, 2 - Transdecoder.cds table 
all_transcripts<-fread(args[1],header=F)
orf_table_transdecoder<-fread(args[2])

setmaes(all_trancripts,"transcript")
setnames(orf_table_transdecoder,c("protein","feature","feature_status","feature_length","strand","score"))

orf_table_transdecoder[,transcript:=str_extract(protein,".*(?=\\.p.+)")]

orf_table_transdecoder_maxlen<-orf_table_transdecoder[,.(max_length=max(feature_length)),by=transcript]

less_than_150_orfs<-unique(orf_table_transdecoder_maxlen[max_length<=150]$transcript)

no_orfs<-setdiff(all_transcripts$transcript,orf_table_transdecoder$transcript)

less_than_150_or_no_orfs<-c(no_orfs,less_than_150_orfs)

write.table(less_than_150_or_no_orfs,file=paste(args[1],"rrna_length_orfs_filtered",sep="_"),row.names = F,col.names = F,quote = F,sep = "\t")