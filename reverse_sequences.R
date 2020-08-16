args <- commandArgs(trailingOnly = TRUE)  
library(Biostrings)
library(IRanges)
library(GenomicRanges)
library(BSgenome)
setwd(".")
#reversing the rev_comp hits
#Arguments: 1 - concatanated fasta files of transcripts with similarities 
sequences_list<-list()
for(i in list.files()) {
  sequences_list<-c(sequences_list,readDNAStringSet(i,format = "fasta"))
}

for (i in 1:length(sequences_list)){
  for(j in 1:length(sequences_list[[i]])) {
    if(score(pairwiseAlignment(sequences_list[[i]][j],sequences_list[[i]][length(sequences_list[[i]])])) < score(pairwiseAlignment(sequences_list[[i]][j],reverseComplement(sequences_list[[i]][length(sequences_list[[i]])])))) {
      sequences_list[[i]][length(sequences_list[[i]])]<-reverseComplement(sequences_list[[i]][length(sequences_list[[i]])])
    }
      
  }
  writeXStringSet(x = sequences_list[[i]],filepath = paste(names(sequences_list[[i]][length(sequences_list[[i]])]),"R.fa",sep="_"),)
}



