args <- commandArgs(trailingOnly = TRUE)  
library(data.table)
library(IRanges)
library(GenomicRanges)
library(DescTools)
library(dplyr)
library(ggplot2)
setwd(".")
# finding the lncRNA consensus - R script which reduces transcripts to genes, finds unique genes and between-transcriptome gene overlaps, takes only unique genes, longer genes from overlaps longer than 20% of the length of the longer gene in the pair and longer + shorter genes from overlaps shorter than 20% of the length of the longer gene in the pair
#Arguments: 1 - filtered rnaspades mapping (paf), 2 - filtered trinity mapping (paf)

#rnaspades and trinity union

rnaspades_paf_filtered<-fread(args[1])
trinity_paf_filtered<-fread(args[2])

rnaspades_paf_filtered[,gene_name:=str_extract(query_name,"NODE_\\d*")]
trinity_paf_filtered[,gene_name:=str_extract(query_name,".*(?=_i\\d*)")]
rnaspades_paf_filtered_granges<-makeGRangesFromDataFrame(rnaspades_paf_filtered,seqnames.field = "target_name",start.field = "target_start",end.field = "target_end",ignore.strand=T,keep.extra.columns = T)
trinity_paf_filtered_granges<-makeGRangesFromDataFrame(trinity_paf_filtered,seqnames.field = "target_name",start.field = "target_start",end.field = "target_end",ignore.strand=T,keep.extra.columns = T)


#finding genes
rnaspades_paf_filtered_granges_split<-split(rnaspades_paf_filtered_granges,rnaspades_paf_filtered_granges$gene_name)
trinity_paf_filtered_granges_split<-split(trinity_paf_filtered_granges,trinity_paf_filtered_granges$gene_name)

rnaspades_paf_filtered_granges_reduced<-reduce(rnaspades_paf_filtered_granges_split,ignore.strand=T)
trinity_paf_filtered_granges_reduced<-reduce(trinity_paf_filtered_granges_split,ignore.strand=T)

rnaspades_paf_filtered_granges_reduced<-unlist(rnaspades_paf_filtered_granges_reduced,use.names = T)
trinity_paf_filtered_granges_reduced<-unlist(trinity_paf_filtered_granges_reduced,use.names = T)
rnaspades_paf_filtered_granges_reduced$gene<-names(rnaspades_paf_filtered_granges_reduced)
trinity_paf_filtered_granges_reduced$gene<-names(trinity_paf_filtered_granges_reduced)

#removing the nested overlaps of genes which don't overlap between the transcriptomes
rnaspades_paf_filtered_granges_reduced$ID<-1:length(rnaspades_paf_filtered_granges_reduced)
trinity_paf_filtered_granges_reduced$ID<-1:length(trinity_paf_filtered_granges_reduced)


nested_overlaps_rnaspades<-findOverlaps(rnaspades_paf_filtered_granges_reduced,rnaspades_paf_filtered_granges_reduced,type = "within")
nested_overlap_table_rnaspades<-cbind(as.data.table(rnaspades_paf_filtered_granges_reduced[queryHits(nested_overlaps_rnaspades)]),as.data.table(rnaspades_paf_filtered_granges_reduced[subjectHits(nested_overlaps_rnaspades)]))

setnames(nested_overlap_table_rnaspades,c("seqnames1","start1","end1","width1","strand1","gene1","ID1","seqnames2","start2","end2","width2","strand2","gene2","ID2"))

inter_width_vector<-c()
for (i in 1:nrow(nested_overlap_table_rnaspades)) {
  inter_width_vector<-c(inter_width_vector,1+Overlap(c(nested_overlap_table_rnaspades$start1[i],nested_overlap_table_rnaspades$end1[i]),c(nested_overlap_table_rnaspades$start2[i],nested_overlap_table_rnaspades$end2[i])))
}
nested_overlap_table_rnaspades[,intersect_widths:=inter_width_vector]
nested_overlap_table_rnaspades[width1<=width2,max_width:=width2]
nested_overlap_table_rnaspades[width1>width2,max_width:=width1]
nested_overlap_table_rnaspades[,length_frac:=100*(intersect_widths/max_width)]
nested_overlap_table_rnaspades<-nested_overlap_table_rnaspades[length_frac!=100]
nested_overlap_table_rnaspades_genes<-unique(nested_overlap_table_rnaspades[,.SD,.SDcols=c("seqnames1","start1","end1","ID1")])

nested_overlaps_trinity<-findOverlaps(trinity_paf_filtered_granges_reduced,trinity_paf_filtered_granges_reduced,type = "within")
nested_overlap_table_trinity<-cbind(as.data.table(trinity_paf_filtered_granges_reduced[queryHits(nested_overlaps_trinity)]),as.data.table(trinity_paf_filtered_granges_reduced[subjectHits(nested_overlaps_trinity)]))

setnames(nested_overlap_table_trinity,c("seqnames1","start1","end1","width1","strand1","gene1","ID1","seqnames2","start2","end2","width2","strand2","gene2","ID2"))

inter_width_vector<-c()
for (i in 1:nrow(nested_overlap_table_trinity)) {
  inter_width_vector<-c(inter_width_vector,1+Overlap(c(nested_overlap_table_trinity$start1[i],nested_overlap_table_trinity$end1[i]),c(nested_overlap_table_trinity$start2[i],nested_overlap_table_trinity$end2[i])))
}
nested_overlap_table_trinity[,intersect_widths:=inter_width_vector]
nested_overlap_table_trinity[width1<=width2,max_width:=width2]
nested_overlap_table_trinity[width1>width2,max_width:=width1]
nested_overlap_table_trinity[,length_frac:=100*(intersect_widths/max_width)]
nested_overlap_table_trinity<-nested_overlap_table_trinity[length_frac!=100]
nested_overlap_table_trinity_genes<-unique(nested_overlap_table_trinity[,.SD,.SDcols=c("seqnames1","start1","end1","ID1")])



#finding between-transcriptome gene overlaps

reduced_overlaps<-findOverlaps(rnaspades_paf_filtered_granges_reduced,trinity_paf_filtered_granges_reduced,ignore.strand=T)
overlap_table<-cbind(as.data.table(rnaspades_paf_filtered_granges_reduced[queryHits(reduced_overlaps)]),as.data.table(trinity_paf_filtered_granges_reduced[subjectHits(reduced_overlaps)]))

setnames(overlap_table,c("seqnames1","start1","end1","width1","strand1","gene1","IDR","seqnames2","start2","end2","width2","strand2","gene2","IDT"))

overlap_table<-overlap_table[IDR %in% nested_overlap_table_rnaspades_genes$ID==F]
overlap_table<-overlap_table[IDT %in% nested_overlap_table_trinity_genes$ID1==F]



overlap_table$IDR<-NULL
overlap_table$IDT<-NULL

inter_width_vector<-c()
for (i in 1:nrow(overlap_table)) {
  inter_width_vector<-c(inter_width_vector,Overlap(c(overlap_table$start1[i],overlap_table$end1[i]),c(overlap_table$start2[i],overlap_table$end2[i])))
}
overlap_table[,intersect_widths:=inter_width_vector]
overlap_table[width1<=width2,max_width:=width2]
overlap_table[width1>width2,max_width:=width1]


overlap_table[,length_frac:=100*(intersect_widths/max_width)]

#intersect plot

perc_plot<-ggplot(data=overlap_table,aes(x=length_frac))+
  geom_histogram(fill="lightsalmon1",col="black",alpha=0.9,binwidth = 5) +
  theme_bw() +
  xlab("Omjer presjeka dvaju gena i duljine dužeg gena (%)") +
  ylab("Broj gena")

ggsave(perc_plot,file="between_transcriptomes_gene_intersect_plot.jpg",width = 5,height = 4)

# Finding longer genes from overlaps longer than 20% of the length of the longer gene in the pair and longer + shorter genes from overlaps shorter than 20% of the length of the longer gene in the pair
overlap_table_filtered_more_than_20_length_frac<-overlap_table[length_frac>=20]
overlap_table_filtered_less_than_20_length_frac<-overlap_table[length_frac<20]

overlap_table_filtered_more_than_20_length_frac[,':='(first_longer=ifelse(width1==max_width,T,F),second_longer=ifelse(width2==max_width,T,F))]

reduced_morethan20_longer_rnaspades<-overlap_table_filtered_more_than_20_length_frac[first_longer==T,.SD,.SDcols=c("seqnames1","start1","end1","width1","gene1")]
reduced_morethan20_longer_trinity<-overlap_table_filtered_more_than_20_length_frac[second_longer==T & first_longer!=T,.SD,.SDcols=c("seqnames2","start2","end2","width2","gene2")]
setnames(reduced_morethan20_longer_rnaspades,c("seqnames","start","end","width","gene"))
setnames(reduced_morethan20_longer_trinity,c("seqnames","start","end","width","gene"))
reduced_morethan20_longer_together<-unique(rbind(reduced_morethan20_longer_rnaspades,reduced_morethan20_longer_trinity))

reduced_morethan20_shorter_rnaspades<-overlap_table_filtered_more_than_20_length_frac[first_longer==F,.SD,.SDcols=c("seqnames1","start1","end1","width1","gene1")]
reduced_morethan20_shorter_trinity<-overlap_table_filtered_more_than_20_length_frac[second_longer==F,.SD,.SDcols=c("seqnames2","start2","end2","width2","gene2")]
setnames(reduced_morethan20_shorter_rnaspades,c("seqnames","start","end","width","gene"))
setnames(reduced_morethan20_shorter_trinity,c("seqnames","start","end","width","gene"))
reduced_morethan20_shorter_together<-unique(rbind(reduced_morethan20_shorter_rnaspades,reduced_morethan20_shorter_trinity))


overlap_table_filtered_less_than_20_length_frac[,':='(first_longer=ifelse(width1==max_width,T,F),second_longer=ifelse(width2==max_width,T,F))]

reduced_lessthan20_longer_rnaspades<-unique(overlap_table_filtered_less_than_20_length_frac[first_longer==T,.SD,.SDcols=c("seqnames1","start1","end1","width1","gene1")])
reduced_lessthan20_longer_trinity<-unique(overlap_table_filtered_less_than_20_length_frac[second_longer==T,.SD,.SDcols=c("seqnames2","start2","end2","width2","gene2")])
setnames(reduced_lessthan20_longer_rnaspades,c("seqnames","start","end","width","gene"))
setnames(reduced_lessthan20_longer_trinity,c("seqnames","start","end","width","gene"))
reduced_lessthan20_longer_together<-unique(rbind(reduced_lessthan20_longer_rnaspades,reduced_lessthan20_longer_trinity))

reduced_lessthan20_shorter_rnaspades<-unique(overlap_table_filtered_less_than_20_length_frac[first_longer==F,.SD,.SDcols=c("seqnames1","start1","end1","width1","gene1")])
reduced_lessthan20_shorter_trinity<-unique(overlap_table_filtered_less_than_20_length_frac[second_longer==F,.SD,.SDcols=c("seqnames2","start2","end2","width2","gene2")])
setnames(reduced_lessthan20_shorter_rnaspades,c("seqnames","start","end","width","gene"))
setnames(reduced_lessthan20_shorter_trinity,c("seqnames","start","end","width","gene"))
reduced_lessthan20_shorter_together<-unique(rbind(reduced_lessthan20_shorter_rnaspades,reduced_lessthan20_shorter_trinity))


#Finding unique genes 
diff_rnaspades_less_longer<-anti_join(as.data.table(rnaspades_paf_filtered_granges_reduced)[,.SD,.SDcols=c("seqnames","start","end","width","gene")],reduced_lessthan20_longer_rnaspades)
diff_rnaspades_less_and_more_longer<-anti_join(diff_rnaspades_less_longer,reduced_morethan20_longer_rnaspades)
diff_rnaspades_less_and_more_longer_less_shorter<-anti_join(diff_rnaspades_less_and_more_longer,reduced_lessthan20_shorter_rnaspades)
diff_rnaspades_less_and_more_longer_less_and_more_shorter<-anti_join(diff_rnaspades_less_and_more_longer_less_shorter,reduced_morethan20_shorter_rnaspades)

diff_trinity_less_longer<-anti_join(as.data.table(trinity_paf_filtered_granges_reduced)[,.SD,.SDcols=c("seqnames","start","end","width","gene")],reduced_lessthan20_longer_trinity)
diff_trinity_less_and_more_longer<-anti_join(diff_trinity_less_longer,reduced_morethan20_longer_trinity)
diff_trinity_less_and_more_longer_less_shorter<-anti_join(diff_trinity_less_and_more_longer,reduced_lessthan20_shorter_trinity)
diff_trinity_less_and_more_longer_less_and_more_shorter<-anti_join(diff_trinity_less_and_more_longer_less_shorter,reduced_morethan20_shorter_trinity)


#lncRNAs which I will use

lncrnas_final_reduced_rnaspades<-unique(makeGRangesFromDataFrame(rbind(reduced_morethan20_longer_rnaspades,reduced_lessthan20_longer_rnaspades,reduced_lessthan20_shorter_rnaspades,diff_rnaspades_less_and_more_longer_less_and_more_shorter),keep.extra.columns = T))
lncrnas_final_reduced_rnaspades$gene_id<-1:length(lncrnas_final_reduced_rnaspades)
lncrnas_final_reduced_trinity<-unique(makeGRangesFromDataFrame(rbind(reduced_morethan20_longer_trinity,reduced_lessthan20_longer_trinity,reduced_lessthan20_shorter_trinity,diff_trinity_less_and_more_longer_less_and_more_shorter),keep.extra.columns = T))
#removing the 100% overlaps from one of the transcriptomes (Trinity)
total_overlap<-overlap_table_filtered_more_than_20_length_frac[first_longer==second_longer]
setnames(total_overlap,c("seqnames1","start1","end1","width1","strand1","gene1","seqnames","start","end","width","strand","gene","intersect_widths","max_width","length_frac","first_longer","second_longer"))
lncrnas_final_reduced_trinity_dt<-as.data.table(lncrnas_final_reduced_trinity)
lncrnas_final_reduced_trinity_dt_anti<-anti_join(lncrnas_final_reduced_trinity_dt,total_overlap[,.SD,.SDcols=c("seqnames","start","end","width","strand","gene")])
lncrnas_final_reduced_trinity<-makeGRangesFromDataFrame(lncrnas_final_reduced_trinity_dt_anti,keep.extra.columns = T)
lncrnas_final_reduced_together<-c(lncrnas_final_reduced_rnaspades,lncrnas_final_reduced_trinity)

#some of the 100% overlaps were missed by total_overlap<-overlap_table_filtered_more_than_20_length_frac[first_longer==second_longer] - removing those as well
lncrnas_final_reduced_together_dt<-as.data.table(lncrnas_final_reduced_together)
lncrnas_final_reduced_together_dt<-lncrnas_final_reduced_together_dt[,.SD,.SDcols=c("seqnames","start","end","width","strand")]
lncrnas_final_reduced_together_dt_dup<-lncrnas_final_reduced_together_dt[duplicated(lncrnas_final_reduced_together_dt)]
lncrnas_final_reduced_together_dt_dup_gene<-merge(lncrnas_final_reduced_trinity_dt_anti,lncrnas_final_reduced_together_dt_dup)

lncrnas_final_reduced_trinity_dt_anti<-anti_join(lncrnas_final_reduced_trinity_dt_anti,lncrnas_final_reduced_together_dt_dup_gene)
lncrnas_final_reduced_trinity<-makeGRangesFromDataFrame(lncrnas_final_reduced_trinity_dt_anti,keep.extra.columns = T)
lncrnas_final_reduced_trinity$gene_id<-(length(lncrnas_final_reduced_rnaspades)+1):(length(lncrnas_final_reduced_trinity)+length(lncrnas_final_reduced_rnaspades))
lncrnas_final_reduced_together<-unique(c(lncrnas_final_reduced_rnaspades,lncrnas_final_reduced_trinity))
lncrnas_final_reduced_together_dt<-as.data.table(lncrnas_final_reduced_together)
write.table(lncrnas_final_reduced_together_dt,file="lncrnas_final_reduced_together_dt",row.names=F,qzote=F)

#numbers
all_together_reduced_num<-length(rnaspades_paf_filtered_granges_reduced)+length(trinity_paf_filtered_granges_reduced)

rnaspades_reduced_length<-length(rnaspades_paf_filtered_granges_reduced)

trinity_reduced_length<-length(trinity_paf_filtered_granges_reduced)

same_overlaps_reduced_num<-nrow(total_overlap)

overlap_rnaspades_morethan20_longer_reduced_num<-nrow(reduced_morethan20_longer_rnaspades)

overlap_trinity_morethan20_longer_reduced_num<-nrow(reduced_morethan20_longer_trinity)

overlap_rnaspades_lessthan20_longerreduced_num<-nrow(reduced_lessthan20_longer_rnaspades)

overlap_trinity_lessthan20_longer_reduced_num<-nrow(reduced_lessthan20_longer_trinity)

overlap_together_morethan20_longer_reduced_num<-nrow(reduced_morethan20_longer_together)

overlap_together_lessthan20_longer_reduced_num<-nrow(reduced_lessthan20_longer_together)

overlap_rnaspades_morethan20_shorter_reduced_num<-nrow(reduced_morethan20_shorter_rnaspades)

overlap_trinity_morethan20_shorter_reduced_num<-nrow(reduced_morethan20_shorter_trinity)

overlap_rnaspades_lessthan20_shorter_reduced_num<-nrow(reduced_lessthan20_shorter_rnaspades)

overlap_trinity_lessthan20_shorter_reduced_num<-nrow(reduced_lessthan20_shorter_trinity)

overlap_together_morethan20_shorter_reduced_num<-nrow(reduced_morethan20_shorter_together)

overlap_together_lessthan20_shorter_reduced_num<-nrow(reduced_lessthan20_shorter_together)

rnaspades_unique<-nrow(diff_rnaspades_less_and_more_longer_less_and_more_shorter)

trinity_unique<-nrow(diff_trinity_less_and_more_longer_less_and_more_shorter)

#final number of lncRNA genes
#rnaspades
lncrnas_final_reduced_num_rnaspades<-length(lncrnas_final_reduced_rnaspades)
#trinity
lncrnas_final_reduced_num_trinity<-length(lncrnas_final_reduced_trinity)
#together
lncrnas_final_reduced_num_together<-length(lncrnas_final_reduced_together)

#finding the transcripts

overlaps_isoform_reduced_rnaspades<-findOverlaps(rnaspades_paf_filtered_granges,lncrnas_final_reduced_rnaspades,ignore.strand=T)
overlaps_isoform_reduced_trinity<-findOverlaps(trinity_paf_filtered_granges,lncrnas_final_reduced_trinity,ignore.strand=T)

isoform_overlap_table_rnaspades<-cbind(as.data.table(rnaspades_paf_filtered_granges[queryHits(overlaps_isoform_reduced_rnaspades)]),as.data.table(lncrnas_final_reduced_rnaspades[subjectHits(overlaps_isoform_reduced_rnaspades)]))
setnames(isoform_overlap_table_rnaspades,c("seqnames1","start1","end1","width1","strand1","transcript_name","transcript_length","transcript_start","transcript_end","transcript_strand","contig_length","matches","alignment_block","mapping_quality","perc_matches","gene_tr","seqnames2","start2","end2","width2","strand2","gene","gene_id"))

isoform_overlap_table_trinity<-cbind(as.data.table(trinity_paf_filtered_granges[queryHits(overlaps_isoform_reduced_trinity)]),as.data.table(lncrnas_final_reduced_trinity[subjectHits(overlaps_isoform_reduced_trinity)]))
setnames(isoform_overlap_table_trinity,c("seqnames1","start1","end1","width1","strand1","transcript_name","transcript_length","transcript_start","transcript_end","transcript_strand","contig_length","matches","alignment_block","mapping_quality","perc_matches","gene_tr","seqnames2","start2","end2","width2","strand2","gene","gene_id"))


final_lncrna_list_rnaspades<-unique(isoform_overlap_table_rnaspades$transcript_name)
final_lncrna_list_trinity<-unique(isoform_overlap_table_trinity$transcript_name)
final_lncrna_list_rnaspades_length<-length(final_lncrna_list_rnaspades)
final_lncrna_list_trinity_length<-length(final_lncrna_list_trinity)

final_lncrna_list_together_length<-length(final_lncrna_list_rnaspades) + length(final_lncrna_list_trinity)

#numbers table

lncrna_numbers_table<-data.table(Kategorija=c("svi geni","RNASpades geni","Trinity geni","geni s potpunim preklapanjem", "duži RNASpades geni s više od 20 % preklapaja","duži Trinity geni s više od 20 % preklapaja","duži RNASpades geni s manje od 20 % preklapaja","duži Trinity geni s manje od 20 % preklapaja","svi duži geni s više od 20% preklapanja","svi duži geni s manje od 20% preklapanja","kraæi RNASpades geni s više od 20 % preklapaja","kraæi Trinity geni s više od 20 % preklapaja","kraæi RNASpades geni s manje od 20 % preklapaja","kraæi Trinity geni s manje od 20 % preklapaja","svi kraæi geni s više od 20% preklapanja","svi kraæi geni s manje od 20% preklapanja","RNASpades geni bez preklapanja","Trinity geni bez preklapanja","konaèni RNASpades geni","Konaèni Trinity geni","Konaèni geni","konaèni RNASpades transkripti","konaèni Trinity transkripti","konaèni transkripti"),Broj=c(all_together_reduced_num,rnaspades_reduced_length,trinity_reduced_length,same_overlaps_reduced_num,overlap_rnaspades_morethan20_longer_reduced_num,overlap_trinity_morethan20_longer_reduced_num,overlap_rnaspades_lessthan20_longerreduced_num,overlap_trinity_lessthan20_longer_reduced_num,overlap_together_morethan20_longer_reduced_num,overlap_together_lessthan20_longer_reduced_num,overlap_rnaspades_morethan20_shorter_reduced_num,overlap_trinity_morethan20_shorter_reduced_num,overlap_rnaspades_lessthan20_shorter_reduced_num,overlap_trinity_lessthan20_shorter_reduced_num,overlap_together_morethan20_shorter_reduced_num,overlap_together_lessthan20_shorter_reduced_num,rnaspades_unique,trinity_unique,lncrnas_final_reduced_num_rnaspades,lncrnas_final_reduced_num_trinity,lncrnas_final_reduced_num_together,final_lncrna_list_rnaspades_length,final_lncrna_list_trinity_length,final_lncrna_list_together_length))

lncrna_numbers_table

#writing the lists of consensus transcripts
write.table(final_lncrna_list_rnaspades,file="final_lncrna_list_rnaspades.txt",row.names = F,col.names = F,quote = F,sep = "\t")
write.table(final_lncrna_list_trinity,file="final_lncrna_list_trinity.txt",row.names = F,col.names = F,quote = F,sep = "\t")
#writing the summary table
write.table(lncrna_numbers_table,file="lncrna_numbers_table.txt",row.names = F,quote = F)
#writing the gene-transcript overlap tables for both transcriptomes
write.table(isoform_overlap_table_rnaspades,file="gene_transcript_overlap_table_rnaspades.txt",row.names=F,quote=F)
write.table(isoform_overlap_table_trinity,file="gene_transcript_overlap_table_rnaspades.txt",row.names=F,quote=F)

