args <- commandArgs(trailingOnly = TRUE)  
library(data.table)
library(IRanges)
library(GenomicRanges)
library(ggplot2)
library(Biostrings)
library(BSgenome)
setwd(".")
#analysis of the main characteristics of lncRNA (getting exons, introns, element numbers, element lengths, GC content, splice sites)

#Arguments: 1 - rnaspades filtered mapping (paf), 2 - trinity filtered mapping (paf), 3 - rnaspades filtered transcript list, 4 - rnaspades filtered transcript list,  5 - gene transcript overlap table rnaspades, 6 - gene transcript overlap table trinity, 7 rnaspades filtered exons,  8 - trinity filtered exons, argument 9 - genome assembly (fasta), 10 - protein coding annotation with filtered bacterial scaffolds (gtf), 11 - rnaspades transcriptome (fasta), 12 - trinity transcriptome (fasta)

rnaspades_paf_filtered<-fread(args[1])
rnaspades_paf_filtered<-fread(args[2])
final_lncrna_list_rnaspades<-fread(args[3],header = F)
final_lncrna_list_trinity<-fread(args[4],header = F)
isoform_overlap_table_rnaspades<-fread(args[5])
isoform_overlap_table_trinity<-fread(args[6])
rnaspades_filtered_exons_copy<-fread(args[7])
trinity_filtered_exons_copy<-fread(args[8])
esu_genome<-readDNAStringSet(args[9])
protein_coding_annotation_clean<fread(args[10])
all_rnaspades_sequences<-readDNAStringSet(args[11])
all_trinity_sequences<-readDNAStringSet(args[12])

names(all_trinity_sequences)<-str_extract(names(all_trinity_sequences),"TRINITY_.*_i\\d*")
lncrna_rnaspades_sequences<-all_rnaspades_sequences[names(all_rnaspades_sequences)%in%rnaspades_longest_iso]
lncrna_trinity_sequences<-all_trinity_sequences[names(all_trinity_sequences)%in%trinity_longest_iso]

lncrna_together_sequences<-c(lncrna_rnaspades_sequences,lncrna_trinity_sequences)

rnaspades_paf_filtered_union<-rnaspades_paf_filtered[query_name%in%final_lncrna_list_rnaspades]
trinity_paf_filtered_union<-trinity_paf_filtered[query_name%in%final_lncrna_list_trinity]

#writing all isoforms
rnaspades_paf_filtered_more_iso<-rnaspades_paf_filtered[query_name%in%final_lncrna_list_rnaspades]
trinity_paf_filtered_more_iso<-trinity_paf_filtered[query_name%in%final_lncrna_list_trinity]
write.table(rnaspades_paf_filtered_more_iso,file="rnaspades_paf_filtered_more_iso.txt",row.names=F,quote=F)
write.table(trinity_paf_filtered_more_iso,file="trinity_paf_filtered_more_iso.txt",row.names=F,quote=F)

#isoforms


#wathcing for transcript names given by the assemblers
isoform_overlap_table_rnaspades<-isoform_overlap_table_rnaspades[gene_tr==gene]
isoform_overlap_table_trinity<-isoform_overlap_table_trinity[gene_tr==gene]
isoform_map_table_rnaspades<-isoform_overlap_table_rnaspades[,.SD,.SDcols=c("transcript_name","gene_id")]
isoform_map_table_trinity<-isoform_overlap_table_trinity[,.SD,.SDcols=c("transcript_name","gene_id")]
isoform_map_table_together<-rbind(isoform_map_table_rnaspades,isoform_map_table_trinity)
write.table(isoform_map_table_together,file="isoform_map_table_together.txt",row.names=F,quote=F)
isoforms_number_rnaspades<-isoform_overlap_table_rnaspades[,.(number_of_isoforms=.N),by=c("seqnames2","start2","end2","gene")]
isoforms_number_trinity<-isoform_overlap_table_trinity[,.(number_of_isoforms=.N),by=c("seqnames2","start2","end2","gene")]
#lncRNA isoforms
isoforms_number_together_lncrna<-rbind(isoforms_number_rnaspades,isoforms_number_trinity)
isoforms_number_together_lncrna<-isoforms_number_together_lncrna[,.SD,.SDcols="number_of_isoforms"]
isoforms_number_together_lncrna_num<-isoforms_number_together_lncrna[,.N,by=number_of_isoforms]
isoforms_number_together_lncrna_num[,category:="lncRNA"]
isoforms_number_together_lncrna_num[,proportion:=N/sum(N)]
isoforms_number_together_lncrna_num_copy<-copy(isoforms_number_together_lncrna_num)
isoforms_number_together_lncrna_num<-isoforms_number_together_lncrna_num[number_of_isoforms<10]
isoforms_number_together_lncrna_num<-rbind(isoforms_number_together_lncrna_num,data.table(number_of_isoforms=">9",N=sum(isoforms_number_together_lncrna_num_copy[number_of_isoforms>9]$N),category="lncRNA",proportion=sum(isoforms_number_together_lncrna_num_copy[number_of_isoforms>9]$proportion)))
#protein
prot_transcript_for_isoforms<-protein_coding_annotation[feature=="transcript"]
prot_transcript_for_isoforms[,gene:=str_remove(attribute,"\\.t\\d+")]
prot_transcript_for_isoforms_table<-prot_transcript_for_isoforms[,.(number_of_isoforms=.N),by=gene]
prot_transcript_for_isoforms_table<-prot_transcript_for_isoforms_table[,.SD,.SDcols="number_of_isoforms"]
prot_transcript_for_isoforms_table_num<-prot_transcript_for_isoforms_table[,.N,by=number_of_isoforms]
prot_transcript_for_isoforms_table_num[,category:="protein_coding_genes"]
prot_transcript_for_isoforms_table_num[,proportion:=N/sum(N)]
prot_transcript_for_isoforms_table_num_copy<-copy(prot_transcript_for_isoforms_table_num)
prot_transcript_for_isoforms_table_num<-prot_transcript_for_isoforms_table_num[number_of_isoforms<10]
prot_transcript_for_isoforms_table_num<-rbind(prot_transcript_for_isoforms_table_num,data.table(number_of_isoforms=">9",N=sum(prot_transcript_for_isoforms_table_num_copy[number_of_isoforms>9]$N),category="protein_coding_genes",proportion=sum(prot_transcript_for_isoforms_table_num_copy[number_of_isoforms>9]$proportion)))
isoforms_table_together_num<-rbind(isoforms_number_together_lncrna_num,prot_transcript_for_isoforms_table_num)

isoforms_plot_proportions<-ggplot(data=isoforms_table_together_num,aes(x=number_of_isoforms,y=proportion,fill=category)) +
  geom_bar(col="black",position="dodge",stat="identity",width=0.5,alpha=0.9) +
  theme_bw() +
  xlab("Broj izoformi po genu") +
  ylab("Udio gena") +
  scale_x_discrete(limits=c(1:9,">9")) +
  theme(legend.position = "top") +
  scale_fill_manual(values=c("indianred2","steelblue"), name = "Kategorija", labels=c("lncRNA-kodirajuæi geni","protein-kodirajuæi geni")) +
  theme(legend.title = element_text(size=9,face="bold")) 

ggsave(isoforms_plot_proportions,file="isoforms_plot_proportions.jpg",width = 6,height = 4.2)

#exons and introns
##isoforms with the biggest number of exons per gene - if more isoforms have the maximum number of exons, take the one which is longer



rnaspades_filtered_exons_numb_per_trans<-rnaspades_filtered_exons_copy[,.N,by=transcript]
trinity_filtered_exons_numb_per_trans<-trinity_filtered_exons_copy[,.N,by=transcript]
setnames(rnaspades_filtered_exons_numb_per_trans,c("transcript_name","number_of_exons"))
setnames(trinity_filtered_exons_numb_per_trans,c("transcript_name","number_of_exons"))
isoform_overlap_table_rnaspades_merged<-merge(isoform_overlap_table_rnaspades,rnaspades_filtered_exons_numb_per_trans,by="transcript_name")
isoform_overlap_table_trinity_merged<-merge(isoform_overlap_table_trinity,trinity_filtered_exons_numb_per_trans,by="transcript_name")  
isoform_overlap_table_rnaspades_merged[,is_max_exons_iso:=ifelse(.SD==max(.SD),T,F),.SDcols="number_of_exons",by="gene_id"]
isoform_overlap_table_rnaspades_merged[is_max_exons_iso==T,max_len:=ifelse(.SD==max(.SD),T,F),.SDcols="transcript_length",by="gene_id"]
isoform_overlap_table_trinity_merged[,is_max_exons_iso:=ifelse(.SD==max(.SD),T,F),.SDcols="number_of_exons",by="gene_id"]
isoform_overlap_table_trinity_merged[is_max_exons_iso==T,max_len:=ifelse(.SD==max(.SD),T,F),.SDcols="transcript_length",by="gene_id")]
#in the case that more isoforms have the same number of exons and the same length
isoform_overlap_table_rnaspades_merged<-unique(isoform_overlap_table_rnaspades_merged,by=c("gene_id","number_of_exons","transcript_length"))
isoform_overlap_table_trinity_merged<-unique(isoform_overlap_table_trinity_merged,by=c("gene_id","number_of_exons","transcript_length"))

rnaspades_longest_iso<-unique(isoform_overlap_table_rnaspades_merged[max_len==T]$transcript_name)
trinity_longest_iso<-unique(isoform_overlap_table_trinity_merged[max_len==T]$transcript_name)
#writing isoforms with the biggest number of exons (and length in the case of equal number of exons)
write.table(rnaspades_longest_iso,file="rnaspades_longest_iso",row.names = F,col.names = F,quote = F,sep = "\t")
write.table(trinity_longest_iso,file="trinity_longest_iso",row.names = F,col.names = F,quote = F,sep = "\t")
rnaspades_paf_filtered_union<-rnaspades_paf_filtered_union[rnaspades_paf_filtered_union$query_name%in%rnaspades_longest_iso]
trinity_paf_filtered_union<-trinity_paf_filtered_union[trinity_paf_filtered_union$query_name%in%trinity_longest_iso]
filtered_together_paf_union<-rbind(rnaspades_paf_filtered_union,trinity_paf_filtered_union)
write.table(filtered_together_paf_union,file="filtered_together_paf_union.txt",row.names=T,quote=T)
rnaspades_paf_filtered_union_granges<-makeGRangesFromDataFrame(rnaspades_paf_filtered_union,seqnames.field = "target_name",start.field = "target_start",end.field = "target_end",ignore.strand=T,keep.extra.columns = T)
trinity_paf_filtered_union_granges<-makeGRangesFromDataFrame(trinity_paf_filtered_union,seqnames.field = "target_name",start.field = "target_start",end.field = "target_end",ignore.strand=T,keep.extra.columns = T)
filtered_together_paf_union_granges<-c(rnaspades_paf_filtered_union_granges,trinity_paf_filtered_union_granges)


rnaspades_paf_filtered_union_exons<-rnaspades_filtered_exons[transcript%in%rnaspades_longest_iso]
trinity_paf_filtered_union_exons<-trinity_filtered_exons[transcript%in%trinity_longest_iso]
filtered_together_paf_union_exons<-rbind(rnaspades_paf_filtered_union_exons,trinity_paf_filtered_union_exons)
#only metazoan annotation

protein_coding_annotation[,width:=end-start+1]
protein_coding_annotation[feature!="gene",attribute:=str_extract(attribute,'".*t\\d+"')]
protein_coding_annotation[feature!="gene",attribute:=str_remove_all(attribute,'"')]
protein_coding_annotation[feature!="gene",gene:=str_extract(attribute,".*(?=\\.t\\d*)")]
protein_coding_annotation[feature=="gene",gene:=attribute]
#proteins - unique exons - taking the isoform with the largest number of exons (same as for lncRNA genes)
protein_coding_annotation[feature=="CDS" | feature=="transcript",exon_width_sum:=sum(.SD[feature=="CDS"]$width),by=attribute]
protein_coding_annotation[feature=="CDS" | feature=="transcript",number_of_exons:=nrow(.SD)-1,by=attribute]
protein_coding_annotation[feature=="transcript",most_exons_transcript:=ifelse(.SD==max(.SD),T,F),.SDcols="number_of_exons",by=gene]
protein_coding_annotation[most_exons_transcript==T,max_length:=ifelse(.SD==max(.SD),T,F),.SDcols="exon_width_sum",by=gene]
protein_transcripts<-protein_coding_annotation[most_exons_transcript==T]
#writing the protein transcripts name with the largest number of exons (and length if more have the maximum number of exons)
write.table(protein_transcripts,file="protein_transcripts_longest_iso.txt",row.names=F,quote=F)
protein_exons<-protein_coding_annotation[feature=="CDS" & attribute %in% protein_transcripts$attribute]
#number of exons per transcript
number_of_exons_table_lncrna<-filtered_together_paf_union_exons[,.(number_of_exons=.N),by=transcript]
number_of_exons_table_lncrna_num<-number_of_exons_table_lncrna[,.N,by=number_of_exons]
for(i in 1:9) {
  if(i > max(number_of_exons_table_lncrna_num$number_of_exons)) {
    number_of_exons_table_lncrna_num<-rbind(number_of_exons_table_lncrna_num,data.table(number_of_exons=i,N=0))
  }
}
number_of_exons_table_lncrna_num[,category:="lncRNA"]
number_of_exons_table_lncrna_num[,proportion:=N/sum(N)]
number_of_exons_table_lncrna_num<-rbind(number_of_exons_table_lncrna_num,data.table(number_of_exons=paste(">",9,sep=""),N=sum(number_of_exons_table_lncrna_num[number_of_exons>9]$N),category="lncRNA",proportion=sum(number_of_exons_table_lncrna_num[number_of_exons>9]$proportion)))

protein_exons_table<-protein_exons[,.(number_of_exons=.N),by=attribute][order(-number_of_exons)]
protein_exons_table_num<-protein_exons_table[,.N,by=number_of_exons]
protein_exons_table_num[,category:="protein_coding_gene"]
protein_exons_table_num<-protein_exons_table_num[number_of_exons>1]
protein_exons_table_num[,proportion:=N/sum(N)]
protein_exons_table_num_copy<-copy(protein_exons_table_num)
protein_exons_table_num<-protein_exons_table_num[number_of_exons<=9]
protein_exons_table_num<-rbind(protein_exons_table_num,data.table(number_of_exons=paste(">",9,sep=""),N=sum(protein_exons_table_num_copy[number_of_exons>9]$N),category="protein_coding_gene",proportion=sum(protein_exons_table_num_copy[number_of_exons>9]$proportion)))
exon_number_table_together_num<-rbind(number_of_exons_table_lncrna_num,protein_exons_table_num)
exons_plot_proportion<-ggplot(data=exon_number_table_together_num,aes(x=number_of_exons,y=proportion,fill=category)) +
  geom_bar(col="black",position="dodge",stat="identity",width=0.5,alpha=0.9) +
  theme_bw() +
  xlab("Broj egzona po genu") +
  ylab("Udio gena") +
  theme(legend.position = "top") +
  scale_x_discrete(limits=c(2:9,">9"),drop=F) +
  scale_fill_manual(values=c("indianred2","steelblue"), name = "Kategorija", labels=c("lncRNA-kodirajuæi geni","protein-kodirajuæi geni")) +
  theme(legend.title = element_text(size=9,face="bold")) 



ggsave(exons_plot_proportion,file="exons_plot_proportion.jpg",width=6, height = 4.2)



#lncRNA and protein exon width comparison
protein_exons<-makeGRangesFromDataFrame(protein_exons,keep.extra.columns = T)
filtered_together_paf_union_exons[,width:=end-start+1]
exon_widths_table<-data.table(width=c(filtered_together_paf_union_exons$width,width(protein_exons)),category=c(rep("lncRNA",nrow(filtered_together_paf_union_exons)),rep("protein_coding_genes",length(protein_exons))),feature=rep("egzoni",nrow(filtered_together_paf_union_exons)+length(protein_exons)))
wilcox.test(exon_widths_table[category=="lncRNA"]$width,exon_widths_table[category=="protein_coding_genes"]$width)
#lncRNA introns
rnaspades_filtered_exons_granges<-makeGRangesFromDataFrame(rnaspades_paf_filtered_union_exons,seqnames.field = "chromosome",keep.extra.columns = T)
trinity_filtered_exons_granges<-makeGRangesFromDataFrame(trinity_paf_filtered_union_exons,seqnames.field = "chromosome",keep.extra.columns = T)
together_filtered_exons_granges<-c(rnaspades_filtered_exons_granges,trinity_filtered_exons_granges)
rnaspades_filtered_exons_granges_split<-split(rnaspades_filtered_exons_granges,rnaspades_filtered_exons_granges$transcript)
trinity_filtered_exons_granges_split<-split(trinity_filtered_exons_granges,trinity_filtered_exons_granges$transcript)

rnaspades_paf_filtered_union_granges_split<-split(rnaspades_paf_filtered_union_granges[rnaspades_paf_filtered_union_granges$query_name%in%rnaspades_filtered_exons_granges$transcript],rnaspades_paf_filtered_union_granges[rnaspades_paf_filtered_union_granges$query_name%in%rnaspades_filtered_exons_granges$transcript]$query_name)
trinity_paf_filtered_union_granges_split<-split(trinity_paf_filtered_union_granges[trinity_paf_filtered_union_granges$query_name%in%trinity_filtered_exons_granges$transcript],trinity_paf_filtered_union_granges[trinity_paf_filtered_union_granges$query_name%in%trinity_filtered_exons_granges$transcript]$query_name)
rnaspades_filtered_introns<-unlist(GenomicRanges::setdiff(rnaspades_paf_filtered_union_granges_split,rnaspades_filtered_exons_granges_split))
rnaspades_filtered_introns<-rnaspades_filtered_introns[width(rnaspades_filtered_introns)>1]
trinity_filtered_introns<-unlist(GenomicRanges::setdiff(trinity_paf_filtered_union_granges_split,trinity_filtered_exons_granges_split))
trinity_filtered_introns<-trinity_filtered_introns[width(trinity_filtered_introns)>1]

lncrna_filtered_introns<-c(rnaspades_filtered_introns,trinity_filtered_introns)
lncrna_filtered_introns$transcript<-names(lncrna_filtered_introns)
lncrna_filtered_introns<-lncrna_filtered_introns[width(lncrna_filtered_introns)>3]
introns_table<-as.data.table(lncrna_filtered_introns)
introns_table[,lncrna_transcript:=names(lncrna_filtered_introns)]
write.table(introns_table,file="introns_table_together.txt",row.names=F,quote=F)
#protein introns
protein_exons_split<-split(protein_exons,protein_exons$attribute)
protein_transcripts_granges<-makeGRangesFromDataFrame(protein_transcripts,ignore.strand = T,keep.extra.columns = T)
protein_transcripts_granges_split<-split(protein_transcripts_granges,protein_transcripts_granges$attribute)
protein_exons_split<-protein_exons_split[names(protein_exons_split) %in% names(protein_transcripts_granges_split)]
protein_transcripts_granges_split<-protein_transcripts_granges_split[names(protein_transcripts_granges_split)%in%names(protein_exons_split)]
protein_introns<-unlist(GenomicRanges::setdiff(protein_transcripts_granges_split,protein_exons_split))
protein_introns$transcript<-names(protein_introns)
protein_introns<-protein_introns[width(protein_introns)>3]
#lncRNA and protein intron comparison
intron_width_table<-data.table(width=c(width(lncrna_filtered_introns),width(protein_introns)),category=c(rep("lncRNA",length(lncrna_filtered_introns)),rep("protein_coding_genes",length(protein_introns))),feature=rep("introni",length(lncrna_filtered_introns)+length(protein_introns)))
wilcox.test(intron_width_table[category=="lncRNA"]$width,intron_width_table[category=="protein_coding_genes"]$width)

#transcript lengths

#only metazoan scaffolds
esu_genome<-esu_genome[names(esu_genome)%in%taxonomy_table_nonbact]
lncrna_intron_sequences<-getSeq(esu_genome,lncrna_filtered_introns)
protein_exons_sequences<-getSeq(esu_genome,protein_exons)
names(protein_exons_sequences)<-protein_exons$attribute
exon_table<-as.data.table(protein_exons_sequences)
exon_table[,name:=names(protein_exons_sequences)]
exon_table[,transcripts:=lapply(.SD, paste0, collapse=""), by = name]
exon_table_unique<-unique(exon_table,by="transcripts")
protein_transcript_sequences<-DNAStringSet(exon_table_unique$transcripts)
names(protein_transcript_sequences)<-exon_table_unique$name
protein_transcript_lengths<-sapply(protein_transcript_sequences,length)


lncrna_transcript_lengths<-filtered_together_paf_union_granges$query_length
transcript_length_table<-data.table(width=c(lncrna_transcript_lengths,protein_transcript_lengths),category=c(rep("lncRNA",length(lncrna_transcript_lengths)),rep("protein_coding_genes",length(protein_transcript_lengths))),feature="transkripti")

wilcox.test(transcript_length_table[category=="lncRNA"]$width,transcript_length_table[category=="protein_coding_genes"]$width)

#gene length
protein_coding_genes<-protein_coding_annotation[feature=="gene"]
protein_coding_genes_granges<-makeGRangesFromDataFrame(protein_coding_genes,ignore.strand = T,keep.extra.columns = T)
lncrna_gene_lengths<-c(width(rnaspades_paf_filtered_granges_reduced),width(trinity_paf_filtered_granges_reduced))
protein_coding_gene_lengths<-width(protein_coding_genes_granges)
gene_length_table<-data.table(width=c(lncrna_gene_lengths,protein_coding_gene_lengths),category=c(rep("lncRNA",length(lncrna_gene_lengths)),rep("protein_coding_genes",length(protein_coding_gene_lengths))),feature="geni")
wilcox.test(gene_length_table[category=="lncRNA"]$width,gene_length_table[category=="protein_coding_genes"]$width)  
#Element lengths
element_width_table<-rbind(exon_widths_table,intron_width_table,transcript_length_table,gene_length_table)
element_width_table[,feature:=factor(feature,levels = c("egzoni","introni","transkripti","geni"))]
element_length_boxplot<-ggplot(data=element_width_table,aes(x=category,y=width,fill=category)) +
  geom_boxplot(alpha=0.9) +
  facet_grid(.~feature) +
  scale_y_log10() +
  theme_bw() +
  scale_fill_manual(values=c("indianred2","steelblue"),name="Kategorija",label=c("lncRNA-kodirajuæi geni","protein-kodirajuæi geni")) +
  theme(legend.position = "top") +
  theme(legend.title = element_text(size=9,face="bold")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x =element_blank()) +
  ylab("Duljina (log10 skalirano)") 

ggsave(element_length_boxplot,file="element_length_boxplot.jpg",width=6,height = 4)




#splice sites - problematic because of unstranded ranscriptomes 



lncrna_introns_first_two_bases<-subseq(lncrna_intron_sequences,start=1,end=2)
lncrna_introns_last_two_bases<-reverse(subseq(reverse(lncrna_intron_sequences),start=1,end=2))
splice_sites_table_lncrna<-data.table(first_two=as.character(lncrna_introns_first_two_bases),last_two=as.character(lncrna_introns_last_two_bases))
splice_sites_table_total_lncrna<-splice_sites_table_lncrna[,.N,by=c("first_two","last_two")]
splice_sites_table_total_lncrna<-splice_sites_table_total_lncrna[first_two!="NN" & last_two!="NN"]
all_dinucleotides<-colnames(dinucleotideFrequency(esu_genome))
all_dinucleotide_combinations<-as.data.table(expand.grid(all_dinucleotides,all_dinucleotides))
setnames(all_dinucleotide_combinations,c("first_two","last_two"))
missing_dinucleotides_lncrna<-anti_join(all_dinucleotide_combinations,splice_sites_table_total_lncrna[,.SD,.SDcols=1:2])
zeros_lncrna<-data.table(N=rep(0,nrow(missing_dinucleotides_lncrna)))
missing_dinucleotides_frequencies_lncrna<-cbind(missing_dinucleotides_lncrna,zeros_lncrna)
splice_sites_table_total_lncrna<-rbind(splice_sites_table_total_lncrna,missing_dinucleotides_frequencies_lncrna)
splice_sites_table_total_lncrna[,category:="lncRNA"]
splice_sites_table_total_lncrna[,proportion:=N/sum(N)]
protein_intron_sequences<-getSeq(esu_genome,protein_introns)
protein_introns_first_two_bases<-subseq(protein_intron_sequences,start=1,end=2)
protein_introns_last_two_bases<-reverse(subseq(reverse(protein_intron_sequences),start=1,end=2))
splice_sites_table_protein<-data.table(first_two=as.character(protein_introns_first_two_bases),last_two=as.character(protein_introns_last_two_bases))
splice_sites_table_total_protein<-splice_sites_table_protein[,.N,by=c("first_two","last_two")]
splice_sites_table_total_protein<-splice_sites_table_total_protein[first_two!="NN" & last_two!="NN"]
missing_dinucleotides_protein<-anti_join(all_dinucleotide_combinations,splice_sites_table_total_protein[,.SD,.SDcols=1:2])
zeros_protein<-data.table(N=rep(0,nrow(missing_dinucleotides_protein)))
missing_dinucleotides_frequencies_protein<-cbind(missing_dinucleotides_protein,zeros_protein)
splice_sites_table_total_protein<-rbind(splice_sites_table_total_protein,missing_dinucleotides_frequencies_protein)
splice_sites_table_total_protein[,category:="protein_coding_genes"]
splice_sites_table_total_protein[,proportion:=N/sum(N)]


splice_sites_table_total_together<-rbind(splice_sites_table_total_lncrna,splice_sites_table_total_protein)
splice_sites_table_total_together[,category:=ifelse(category=="lncRNA","lncRNA geni","protein-kodirajuæi geni")]
heatmap_splice_sites<-ggplot(data=splice_sites_table_total_together,aes(x=first_two,y=last_two,fill=proportion)) +
  geom_tile(col="black")  +
  facet_wrap(.~category) +
  scale_fill_gradient(low="white",high="indianred2",trans="sqrt",na.value = "white",breaks=c(0,0.1,0.2,0.3,0.5))  +
  xlab("5' mjesto prekrajanja / 3' reverzni komplement mjesta prekrajanja") +
  ylab("3' mjesto prekrajanja / 5' reverzni komplement mjesta prekrajanja") +
  labs(fill = "Udio introna (sqrt skala)") +
  theme(legend.position="top") +
  theme(legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(0.6,"cm")) +
  theme(legend.title = element_text( size = 8)) +
  theme(axis.title.x = element_text(size=8)) +
  theme(axis.title.y = element_text(size=8)) +
  theme(text = element_text(size=6)) +
  theme(legend.title = element_text(size=9,face="bold")) 



ggsave(heatmap_splice_sites,file="heatmap_splice_sites.jpg")

#GC content

names(all_trinity_sequences)<-str_extract(names(all_trinity_sequences),"TRINITY_.*_i\\d*")
lncrna_rnaspades_sequences<-all_rnaspades_sequences[names(all_rnaspades_sequences)%in%rnaspades_longest_iso]
lncrna_trinity_sequences<-all_trinity_sequences[names(all_trinity_sequences)%in%trinity_longest_iso]

lncrna_together_sequences<-c(lncrna_rnaspades_sequences,lncrna_trinity_sequences)

genome_gc<-alphabetFrequency(unlist(esu_genome),as.prob = T)["C"]+alphabetFrequency(unlist(esu_genome),as.prob = T)["G"]
prot_gc<-alphabetFrequency(protein_transcript_sequences,as.prob = T)[,"C"]+alphabetFrequency(protein_transcript_sequences,as.prob = T)[,"G"]
lncrna_gc<-alphabetFrequency(lncrna_together_sequences,as.prob = T)[,"C"]+alphabetFrequency(lncrna_together_sequences,as.prob = T)[,"G"]
gc_table<-data.table(transcript_gc=c(lncrna_gc,prot_gc),category=c(rep("lncRNA",length(lncrna_gc)),rep("protein_coding_transcripts",length(prot_gc))))



gc_density_plot<-ggplot(data=gc_table,aes(x=transcript_gc,y=..scaled..,fill=category)) +
  geom_density(alpha=0.5) +
  theme_bw() +
  theme(legend.position = "top") +
  geom_vline(xintercept = genome_gc,linetype="dotted",size=1) +
  scale_fill_manual(values=c("indianred2","steelblue"), name = "Kategorija", labels=c("lncRNA","mRNA")) +
  xlab("udio GC tranksripta") +
  theme(axis.title.y =  element_blank()) +
  theme(legend.title = element_text(size=9,face="bold"))


ggsave(gc_density_plot,file="gc_density_plot.jpg",width = 6,height = 4.6)
t.test(lncrna_gc,prot_gc)  
