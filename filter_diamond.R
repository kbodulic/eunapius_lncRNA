args <- commandArgs(trailingOnly = TRUE)  
library(data.table)
setwd(".")
#Dimond filtering - Remove transcripts which show a significant similarity to known proteeins (NR database)
#Arguments: 1 - List of transcripts fitlered by rRNA, length and ORF length, 2 - DIAMOND results
all_transcripts<-fread(args[1],header=F)
setmaes(all_trancripts,"transcript")

diamond_trancripts<-fread(args[2])

setnames(diamond_transcripts,c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"))



diamond_transcripts_besthit<-diamond_transcripts[,.(eval_min=min(evalue)),by=qseqid]



less_than_150_ORF_without_blast<-setdiff(all_transcripts,diamond_transcripts_besthit[eval_min<=1e-5]$qseqid)


write.table(less_than_150_ORF_without_blast,file=paste(args[1],"less_than_150_ORF_without_diamond",sep="_"),row.names = F,col.names = F,quote = F,sep = "\t")