args <- commandArgs(trailingOnly = TRUE)  
library(data.table)
setwd(".")


#rRNA filtering
all_transcripts<-fread(args[1])
setmaes(all_trancripts,c("transcript","length"))


blastn_rrna<-fread(args[2])
setnames(blastn_rrna,c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","qlen","sstart","send","slen","evalue","bitscore"))

rrnas<-unique(blastn_rrna[evalue<=1e-20]$qseqid)
transcripts_without_rrna<-all_transcripts[transcript%in%rrnas==F]


#Length filtering - 200 

all_transcripts_200<-transcripts_without_rrna[length>200]
write.table(all_transcripts_200$transcript,file=paste(args[1],"without_rrna_more_than_200",sep="_"),row.names = F,col.names = F,quote = F,sep = "\t")

#getting the information about the number of filtered transcripts
write.table(cbind(nrow(transcripts_without_rrna),nrow(all_transcripts_200)),file=paste(args[1],"info_table_rrna_length",sep="_"),row.names = F,col.names = F,quote = F,sep = "\t")