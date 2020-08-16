args <- commandArgs(trailingOnly = TRUE)  
library(data.table)
library(rhmmer)
setwd(".")


#HMMER filtering script - remove translated transcripts with a significant similarity to conserved protein domains
#Arguments: 1 - table of transcripts filtered by rRNA, length, ORF length and DIAMOND hits, 2 - hmmscan results 
all_trancripts<-fread(args[1],header=F)
setmaes(all_trancripts,"transcript")
hmmer_results<-as.data.table(read_domtblout(args[2]))


hmmer_results[,transcript:=str_extract(query_accession,".*(?=\\.p.+)")]


hmmer_results_besthit<-hmmer_results[,.(eval_min=min(sequence_bias)),by=transcript]



no_hmmer<-setdiff(all_trancripts,hmmer_results_besthit[eval_min<=1e-5]$transcript)


write.table(no_hmmer,file=paste(args[1],"no_hmmer",sep="_"),row.names = F,col.names = F,quote = F,sep = "\t")
