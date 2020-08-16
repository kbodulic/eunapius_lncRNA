args <- commandArgs(trailingOnly = TRUE)  
library(data.table)
library(IRanges)
library(GenomicRanges)
library(DescTools)
library(ggplot2)
library(gameofthrones)
setwd(".")
#analysis of transposone overlaps - counting the number of "within" overlaps betweem tramsposones and both the inside and the 5 ' region of lncrna-coding genes and prot-coding genes, counting the normalized number of transposone bases which overlap 5' region, exons and introns of those genes (5' region = 1 kb, shorter if the gene is on the begining of a scaffold) + counting the normalized number of bases of transposone classes which overlap 5' region, exons and introns of those genes (5' region = 1 kb, shorter if the gene is on the begining of a scaffold)
#Arguments: 1 - lncrna classification table, 2 - table for mapping gene_id to longest isoforms,  argument 3 - rnaspades filtered exons, argument 4 - trinity filtered exons, argument 5 - paf with every rnaspades isoform, argument 6 - paf with every trinity isoform, argument 7 - repeat annotation, argument 8 - lncrna genes, argument 9 - protein-coding annottaion (filtered for bacterial scaffolds)

lncRNA_classification_table<-fread(args[1])
isoform_map_table_together<-fread(args[2])
exons_clean_lncrna_rnaspades<-fread(args[3])
exons_clean_lncrna_trinity<-fread(args[4])
rnaspades_paf_filtered_more_iso<-fread(args[5])
trinity_paf_filtered_more_iso<-fread(args[6])
repeatmasker<-readRDS(args[7])
lncrnas_final_reduced_together<-fread(args[8])
protein_coding_annotation_clean<-fread[args[9]]

#merging the transcript classification table with the information about transcript gene_ids
lncRNA_classification_table_with_gene_id<-merge(lncRNA_classification_table,isoform_map_table_together)

#calculating every intron and exon of lncRNA genes (not only of the biggest-number-of-exons isoforms) 

exons_clean_lncrna_rnaspades[,gene:=str_extract(transcript,"NODE_\\d*")]
exons_clean_lncrna_trinity[,gene:=str_extract(transcript,".*(?=_i\\d*)")]
exons_clean_lncrna_together<-rbind(exons_clean_lncrna_rnaspades,exons_clean_lncrna_trinity)
rnaspades_paf_filtered_more_iso<-makeGRangesFromDataFrame(rnaspades_paf_filtered_more_iso,keep.extra.columns = T,seqnames.field = "target_name",start.field = "target_start",end.field = "target_end")
trinity_paf_filtered_more_iso<-makeGRangesFromDataFrame(trinity_paf_filtered_more_iso,keep.extra.columns = T,seqnames.field = "target_name",start.field = "target_start",end.field = "target_end")
rnaspades_paf_filtered_more_iso_split<-split(rnaspades_paf_filtered_more_iso[rnaspades_paf_filtered_more_iso$query_name%in%exons_clean_lncrna_rnaspades$transcript],rnaspades_paf_filtered_more_iso[rnaspades_paf_filtered_more_iso$query_name%in%exons_clean_lncrna_rnaspades$transcript]$query_name)
trinity_paf_filtered_more_iso_split<-split(trinity_paf_filtered_more_iso[trinity_paf_filtered_more_iso$query_name%in%exons_clean_lncrna_trinity$transcript],trinity_paf_filtered_more_iso[trinity_paf_filtered_more_iso$query_name%in%exons_clean_lncrna_trinity$transcript]$query_name)
exons_clean_lncrna_rnaspades_gr<-makeGRangesFromDataFrame(exons_clean_lncrna_rnaspades,keep.extra.columns = T)
exons_clean_lncrna_rnaspades_gr<-exons_clean_lncrna_rnaspades_gr[exons_clean_lncrna_rnaspades_gr$transcript%in%rnaspades_paf_filtered_more_iso$query_name]
exons_clean_lncrna_rnaspades_gr_split<-split(exons_clean_lncrna_rnaspades_gr,exons_clean_lncrna_rnaspades_gr$transcript)
rnaspades_filtered_introns_clean<-unlist(GenomicRanges::setdiff(rnaspades_paf_filtered_more_iso_split,exons_clean_lncrna_rnaspades_gr_split))
rnaspades_filtered_introns_clean<-rnaspades_filtered_introns_clean[width(rnaspades_filtered_introns_clean)>1]

exons_clean_lncrna_trinity_gr<-makeGRangesFromDataFrame(exons_clean_lncrna_trinity,keep.extra.columns = T)
exons_clean_lncrna_trinity_gr<-exons_clean_lncrna_trinity_gr[exons_clean_lncrna_trinity_gr$transcript%in%trinity_paf_filtered_more_iso$query_name]
exons_clean_lncrna_trinity_gr_split<-split(exons_clean_lncrna_trinity_gr,exons_clean_lncrna_trinity_gr$transcript)
trinity_filtered_introns_clean<-unlist(GenomicRanges::setdiff(trinity_paf_filtered_more_iso_split,exons_clean_lncrna_trinity_gr_split))
trinity_filtered_introns_clean<-trinity_filtered_introns_clean[width(trinity_filtered_introns_clean)>1]

lncrna_filtered_introns_clean<-c(rnaspades_filtered_introns_clean,trinity_filtered_introns_clean)
lncrna_filtered_introns_clean$transcript<-names(lncrna_filtered_introns_clean)
lncrna_filtered_exons_clean<-c(exons_clean_lncrna_rnaspades_gr,exons_clean_lncrna_trinity_gr)
lncrna_filtered_introns_clean_dt<-as.data.table(lncrna_filtered_introns_clean)
isoform_map_table_together_cp<-copy(isoform_map_table_together)
setnames(isoform_map_table_together_cp,c("transcript","gene_id"))
lncrna_filtered_introns_clean_dt<-merge(lncrna_filtered_introns_clean_dt,isoform_map_table_together_cp)
lncrna_filtered_introns_clean_dt_gr<-makeGRangesFromDataFrame(lncrna_filtered_introns_clean_dt,keep.extra.columns = T)
lncrna_filtered_introns_clean_dt_gr_split<-split(lncrna_filtered_introns_clean_dt_gr,lncrna_filtered_introns_clean_dt_gr$gene_id)
lncrna_filtered_introns_clean_dt_gr_split_red<-reduce(lncrna_filtered_introns_clean_dt_gr_split)
lncrna_filtered_introns_clean_dt_gr_split_red_un<-unlist(lncrna_filtered_introns_clean_dt_gr_split,use.names = T)
lncrna_filtered_introns_clean_dt_gr_split_red_un_dt<-as.data.table(lncrna_filtered_introns_clean_dt_gr_split_red_un)
lncrna_filtered_introns_clean_dt_gr_split_red_un_dt_merge<-merge(lncrna_filtered_introns_clean_dt_gr_split_red_un_dt,lncRNA_classification_table_with_gene_id,by="gene_id")
introns_sum_table<-lncrna_filtered_introns_clean_dt_gr_split_red_un_dt_merge[,.(sum_el=sum(width)),by="class"]
lncrna_filtered_exons_clean_dt<-as.data.table(lncrna_filtered_exons_clean)
lncrna_filtered_exons_clean_dt<-merge(lncrna_filtered_exons_clean_dt,isoform_map_table_together_cp)
lncrna_filtered_exons_clean_dt_gr<-makeGRangesFromDataFrame(lncrna_filtered_exons_clean_dt,keep.extra.columns = T)
lncrna_filtered_exons_clean_dt_gr_split<-split(lncrna_filtered_exons_clean_dt_gr,lncrna_filtered_exons_clean_dt_gr$gene_id)
lncrna_filtered_exons_clean_dt_gr_split_red<-reduce(lncrna_filtered_exons_clean_dt_gr_split)
lncrna_filtered_exons_clean_dt_gr_split_red_un<-unlist(lncrna_filtered_exons_clean_dt_gr_split,use.names = T)
lncrna_filtered_exons_clean_dt_gr_split_red_un_dt<-as.data.table(lncrna_filtered_exons_clean_dt_gr_split_red_un)
lncrna_filtered_exons_clean_dt_gr_split_red_un_dt_merge<-merge(lncrna_filtered_exons_clean_dt_gr_split_red_un_dt,lncRNA_classification_table_with_gene_id,by="gene_id")
exons_sum_table<-lncrna_filtered_exons_clean_dt_gr_split_red_un_dt_merge[,.(sum_el=sum(width)),by="class"]

#calculating every protein intron and exon (not only of the biggest-number-of-exons isoform)
protein_coding_annotation_clean[feature!="gene",attribute:=str_extract(attribute,'".*t\\d+"')]
protein_coding_annotation_clean[feature!="gene",attribute:=str_remove_all(attribute,'"')]
protein_coding_annotation_exons_clean<-protein_coding_annotation_clean[feature=="CDS"]
protein_coding_annotation_transcripts_clean<-protein_coding_annotation[feature=="transcript"]
protein_exons<-makeGRangesFromDataFrame(protein_coding_annotation_exons_clean,keep.extra.columns = T)
protein_exons_split<-split(protein_exons,protein_exons$attribute)
protein_transcripts_clean_granges<-makeGRangesFromDataFrame(protein_coding_annotation_transcripts_clean,ignore.strand = T,keep.extra.columns = T)
protein_transcripts_clean_granges_split<-split(protein_transcripts_clean_granges,protein_transcripts_clean_granges$attribute)
protein_exons_clean_split<-protein_exons_split[names(protein_exonsn_split) %in% names(protein_transcripts_clean_granges_split)]
protein_transcripts_clean_granges_split<-protein_transcripts_clean_granges_split[names(protein_transcripts_clean_granges_split)%in%names(protein_exons_split)]
protein_introns<-unlist(GenomicRanges::setdiff(protein_transcripts_clean_granges_split,protein_exons_split))
protein_introns$transcript<-names(protein_introns)

protein_coding_genes_granges<-makeGRangesFromDataFrame(protein_coding_annotation_clean[feature=="gene"],keep.extra.columns = T)

#repeats
repeatmasker<-readRDS("AllRepeatMaskerFiles_filtered.RDS")
setnames(repeatmasker,c("SW_score","perc.sub","perc.del","perc.ins","query_name","start","end","base_past_end","strand_ignore","repeat_name","repeat_start","repeat_end","repeat_left_bases","ID","V16","repeat_class","repeat_family","species","size","library","group","ltrretID","type","elementtype"))
repeatmasker_esu<-repeatmasker[species=="Eunapius"]
repeatmasker_esu<-repeatmasker_esu[repeat_class!="Low complexity" & repeat_class!="Simple repeat"]
repeatmasker_esu_granges<-makeGRangesFromDataFrame(repeatmasker_esu,seqnames.field = "query_name",keep.extra.columns = T,ignore.strand = T)
#defining 5+ regions (1kb) - if a gene's start coordiante is less than 1000, take the range between 1 and the start coordinate
lncrnas_final_reduced_together_1kb<-lncrnas_final_reduced_together[start(lncrnas_final_reduced_together)>1000]
lncrnas_final_reduced_together_1kb_short<-lncrnas_final_reduced_together[start(lncrnas_final_reduced_together)<=1000]
regions_lncrna_short<-lncrnas_final_reduced_together_1kb_short
end(regions_lncrna_short)<-(start(regions_lncrna_short)-1)
start(lncrnas_final_reduced_together_1kb_short)<-(rep(1,length(lncrnas_final_reduced_together_1kb_short)))
start(regions_lncrna_short)<-(rep(1,length(lncrnas_final_reduced_together_1kb_short)))
new_starts<-start(lncrnas_final_reduced_together_1kb)-1000
start(lncrnas_final_reduced_together_1kb)<-new_starts
lncrnas_final_reduced_together_1kb_final<-c(lncrnas_final_reduced_together_1kb,lncrnas_final_reduced_together_1kb_short)
norm_width_1bb_lncrna<-1000*length(lncrnas_final_reduced_together_1kb)+sum(start(lncrnas_final_reduced_together_1kb_short))

protein_final_reduced_together_1kb<-protein_coding_genes_granges[start(protein_coding_genes_granges)>1000]
protein_final_reduced_together_1kb_short<-protein_coding_genes_granges[start(protein_coding_genes_granges)<=1000]
regions_prot_short<-protein_final_reduced_together_1kb_short
end(regions_prot_short)<-(start(regions_prot_short)-1)
start(protein_final_reduced_together_1kb_short)<-(rep(1,length(protein_final_reduced_together_1kb_short)))
start(regions_prot_short)<-(rep(1,length(protein_final_reduced_together_1kb_short)))
new_starts_prot<-start(protein_final_reduced_together_1kb)-1000
start(protein_final_reduced_together_1kb)<-new_starts_prot
protein_final_reduced_together_1kb_final<-c(protein_final_reduced_together_1kb,protein_final_reduced_together_1kb_short)
norm_width_1bb_prot<-1000*length(protein_final_reduced_together_1kb)+sum(start(protein_final_reduced_together_1kb_short))

#only_1kb_regions within overlaps
regions_lncrna_long<-lncrnas_final_reduced_together_1kb
end(regions_lncrna_long)<-(start(regions_lncrna_long)+999)
regions_lncrna_together<-c(regions_lncrna_long,regions_lncrna_short)
regions_prot_long<-protein_final_reduced_together_1kb
end(regions_prot_long)<-(start(regions_prot_long)+999)
regions_prot_together<-c(regions_prot_long,regions_prot_short)
overlaps_repeats_lncrna<-findOverlaps(repeatmasker_esu_granges,lncrnas_final_reduced_together,ignore.strand=T,type="within")
overlaps_repeats_lncrna_table<-cbind(as.data.table(repeatmasker_esu[queryHits(overlaps_repeats_lncrna)]),as.data.table(lncrnas_final_reduced_together[subjectHits(overlaps_repeats_lncrna)]))
setnames(overlaps_repeats_lncrna_table,c("SW_score","perc.sub","perc.del","perc.ins","query_name","start1","end1","base_past_end","strand_ignore","repeat_name","repeat_start","repeat_end","repeat_left_bases","ID","V16","repeat_class","repeat_family","species","size","library","group","ltrret_ID","type","element_type","seqnames2","start2","end2","width2","strand2","gene","gene_id"))
#within genes overlaps
overlaps_repeats_lncrna_1kb<-findOverlaps(repeatmasker_esu_granges,regions_lncrna_together,ignore.strand=T,type="within")
overlaps_repeats_lncrna_table_1kb<-cbind(as.data.table(repeatmasker_esu[queryHits(overlaps_repeats_lncrna_1kb)]),as.data.table(lncrnas_final_reduced_together_1kb_final[subjectHits(overlaps_repeats_lncrna_1kb)]))
setnames(overlaps_repeats_lncrna_table_1kb,c("SW_score","perc.sub","perc.del","perc.ins","query_name","start1","end1","base_past_end","strand_ignore","repeat_name","repeat_start","repeat_end","repeat_left_bases","ID","V16","repeat_class","repeat_family","species","size","library","group","ltrret_ID","type","element_type","seqnames2","start2","end2","width2","strand2","gene","gene_id"))

id_trans_within<-unique(overlaps_repeats_lncrna_table[,.SD,.SDcols=c("ID","gene_id")])
id_trans_region<-unique(overlaps_repeats_lncrna_table_1kb[,.SD,.SDcols=c("ID","gene_id")])



lncrna_repeat_overlap_category_class<-merge(id_trans_within,lncRNA_classification_table_with_gene_id,by="gene_id")
lncrna_repeat_overlap_category_class_reg<-merge(id_trans_region,lncRNA_classification_table_with_gene_id,by="gene_id")

#counting the number of unique lncRNA genes with 5' or intrinsic transposone overlap (only whole transposones)
num_within_genes<-lncrna_repeat_overlap_category_class[,uniqueN(gene_id),by=class]
num_within_regions<-lncrna_repeat_overlap_category_class_reg[,uniqueN(gene_id),by=class]


#protein transposone overlaps
#5' region overlap

overlaps_repeats_proteins<-findOverlaps(repeatmasker_esu_granges,protein_coding_genes_granges,ignore.strand=T,type="within")
overlaps_repeats_proteins_table<-cbind(as.data.table(repeatmasker_esu[queryHits(overlaps_repeats_proteins)]),as.data.table(protein_coding_genes_granges[subjectHits(overlaps_repeats_proteins)]))
setnames(overlaps_repeats_proteins_table,c("SW_score","perc.sub","perc.del","perc.ins","query_name","start1","end1","base_past_end","strand_ignore","repeat_name","repeat_start","repeat_end","repeat_left_bases","ID","V16","repeat_class","repeat_family","species","size","library","group","ltrret_ID","type","element_type","seqnames2","start2","end2","width2","strand2","source","feature","score","strand_ignore2","frame","attribute"))
#intrinsic overlap
overlaps_repeats_proteins_1kb<-findOverlaps(repeatmasker_esu_granges,regions_prot_together,ignore.strand=T,type="within")
overlaps_repeats_proteins_table_1kb<-cbind(as.data.table(repeatmasker_esu[queryHits(overlaps_repeats_proteins_1kb)]),as.data.table(regions_prot_together[subjectHits(overlaps_repeats_proteins_1kb)]))
setnames(overlaps_repeats_proteins_table_1kb,c("SW_score","perc.sub","perc.del","perc.ins","query_name","start1","end1","base_past_end","strand_ignore","repeat_name","repeat_start","repeat_end","repeat_left_bases","ID","V16","repeat_class","repeat_family","species","size","library","group","ltrret_ID","type","element_type","seqnames2","start2","end2","width2","strand2","source","feature","score","strand_ignore2","frame","attribute"))

id_trans_within_prot<-unique(overlaps_repeats_proteins_table[,.SD,.SDcols=c("ID","attribute")])
id_trans_region_prot<-unique(overlaps_repeats_proteins_table_1kb[,.SD,.SDcols=c("ID","attribute")])
#counting the number of unique protein genes with 5' or intrinsic transposone overlap (only whole transposones)
num_within_genes_prot<-data.table(class="protein-coding gene",V1=nrow(unique(id_trans_within_prot,by="attribute")))
num_within_regions_prot<-data.table(class="protein-coding gene",V1=nrow(unique(id_trans_region_prot,by="attribute")))

#binding the lncRNA and protein tables
genes_numb_within<-rbind(num_within_genes,num_within_genes_prot)
genes_numb_region<-rbind(num_within_regions,num_within_regions_prot)
genes_numb_region<-genes_numb_region[order(class)]
total_num_genes<-data.table(total_number_of_genes=c(lncRNA_classification_table[class=="intergenic_lncRNA",.N],lncRNA_classification_table[class=="intronic_lncRNA",.N],lncRNA_classification_table[class=="overlapping_lncRNA",.N],protein_coding_genes[,.N]))
overlaps_repeat_together_num_table<-cbind(genes_numb_within,genes_numb_region$V1,total_num_genes)
setnames(overlaps_repeat_together_num_table,c("class","genes_with_transposones_inside","genes_with_transposones_region","total_genes"))
overlaps_repeat_together_num_table[,percentage_genes_within:=100*genes_with_transposones_inside/total_genes]
overlaps_repeat_together_num_table[,percentage_genes_region:=100*genes_with_transposones_region/total_genes]
#Summing ovwerlap bases - comparison between lncRNA and proteins by classes of repeats 
#lncRNA
#any overlaps

overlaps_repeats_lncrna_bases<-findOverlaps(repeatmasker_esu_granges,lncrnas_final_reduced_together_1kb_final,ignore.strand=T,type="any")
overlaps_repeats_lncrna_bases_table<-cbind(as.data.table(repeatmasker_esu[queryHits(overlaps_repeats_lncrna_bases)]),as.data.table(lncrnas_final_reduced_together[subjectHits(overlaps_repeats_lncrna_bases)]))
setnames(overlaps_repeats_lncrna_bases_table,c("SW_score","perc.sub","perc.del","perc.ins","query_name","start1","end1","base_past_end","strand_ignore","repeat_name","repeat_start","repeat_end","repeat_left_bases","ID","V16","repeat_class","repeat_family","species","size","library","group","ltrret_ID","type","element_type","seqnames2","start2","end2","width2","strand2","gene","gene_id"))
inter_width_vector<-c()
for (i in 1:nrow(overlaps_repeats_lncrna_bases_table)) {
  inter_width_vector<-c(inter_width_vector,Overlap(c(overlaps_repeats_lncrna_bases_table$start1[i],overlaps_repeats_lncrna_bases_table$end1[i]),c(overlaps_repeats_lncrna_bases_table$start2[i],overlaps_repeats_lncrna_bases_table$end2[i])))
}
overlaps_repeats_lncrna_bases_table[,intersect_width:=inter_width_vector]

#repeats exons overlap

overlaps_repeats_exons_20<-findOverlaps(repeatmasker_esu_granges,lncrna_filtered_exons_clean,ignore.strand=T)
overlaps_repeats_exons_table_20<-cbind(as.data.table(repeatmasker_esu_granges[queryHits(overlaps_repeats_exons_20)]),as.data.table(lncrna_filtered_exons_clean[subjectHits(overlaps_repeats_exons_20)]))
setnames(overlaps_repeats_exons_table_20,c("seqnames1","start1","end1","width1","strand1","SW_score","perc.sub","perc.del","perc.ins","base_past_end","strand_ignore","repeat_name","repeat_start","repeat_end","repeat_left_bases","ID","V16","repeat_class","repeat_family","species","size","library","group","ltrret_ID","type","elementtype","seqnames2","start2","end2","width2","strand2","transcript","gene"))
overlaps_repeats_exons_table_20<-merge(overlaps_repeats_exons_table_20,isoform_map_table_together_cp)
inter_width_vector<-c()
for (i in 1:nrow(overlaps_repeats_exons_table_20)) {
  inter_width_vector<-c(inter_width_vector,Overlap(c(overlaps_repeats_exons_table_20$start1[i],overlaps_repeats_exons_table_20$end1[i]),c(overlaps_repeats_exons_table_20$start2[i],overlaps_repeats_exons_table_20$end2[i])))
}
overlaps_repeats_exons_table_20[,intersect_width:=inter_width_vector]
#repeats introns overlap
overlaps_repeats_introns_20<-findOverlaps(repeatmasker_esu_granges,lncrna_filtered_introns_clean,ignore.strand=T)
overlaps_repeats_introns_table_20<-cbind(as.data.table(repeatmasker_esu_granges[queryHits(overlaps_repeats_introns_20)]),as.data.table(lncrna_filtered_introns_clean[subjectHits(overlaps_repeats_introns_20)]))
setnames(overlaps_repeats_introns_table_20,c("seqnames1","start1","end1","width1","strand1","SW_score","perc.sub","perc.del","perc.ins","base_past_end","strand_ignore","repeat_name","repeat_start","repeat_end","repeat_left_bases","ID","V16","repeat_class","repeat_family","species","size","library","group","ltrret_ID","type","elementtype","seqnames2","start2","end2","width2","strand2","transcript"))
overlaps_repeats_introns_table_20<-merge(overlaps_repeats_introns_table_20,isoform_map_table_together_cp)
inter_width_vector<-c()
for (i in 1:nrow(overlaps_repeats_introns_table_20)) {
  inter_width_vector<-c(inter_width_vector,Overlap(c(overlaps_repeats_introns_table_20$start1[i],overlaps_repeats_introns_table_20$end1[i]),c(overlaps_repeats_introns_table_20$start2[i],overlaps_repeats_introns_table_20$end2[i])))
}
overlaps_repeats_introns_table_20[,intersect_width:=inter_width_vector]
#repeat regions overlap (5' 1 kb)
overlaps_repeats_lncrna_bases_regions<-findOverlaps(repeatmasker_esu_granges,regions_lncrna_together,ignore.strand=T,type="any")
overlaps_repeats_lncrna_bases_table_regions<-cbind(as.data.table(repeatmasker_esu[queryHits(overlaps_repeats_lncrna_bases_regions)]),as.data.table(regions_lncrna_together[subjectHits(overlaps_repeats_lncrna_bases_regions)]))
setnames(overlaps_repeats_lncrna_bases_table_regions,c("SW_score","perc.sub","perc.del","perc.ins","query_name","start1","end1","base_past_end","strand_ignore","repeat_name","repeat_start","repeat_end","repeat_left_bases","ID","V16","repeat_class","repeat_family","species","size","library","group","ltrret_ID","type","element_type","seqnames2","start2","end2","width2","strand2","gene","gene_id"))
inter_width_vector<-c()
for (i in 1:nrow(overlaps_repeats_lncrna_bases_table_regions)) {
  inter_width_vector<-c(inter_width_vector,Overlap(c(overlaps_repeats_lncrna_bases_table_regions$start1[i],overlaps_repeats_lncrna_bases_table_regions$end1[i]),c(overlaps_repeats_lncrna_bases_table_regions$start2[i],overlaps_repeats_lncrna_bases_table_regions$end2[i])))
}
overlaps_repeats_lncrna_bases_table_regions[,intersect_width:=inter_width_vector]


lncrna_overlap_exons_20<-unique(data.table(ID=overlaps_repeats_exons_table_20$ID,gene_id=overlaps_repeats_exons_table_20$gene_id,intersect=overlaps_repeats_exons_table_20$intersect_width))

lncrna_overlap_introns_20<-unique(data.table(ID=overlaps_repeats_introns_table_20$ID,gene_id=overlaps_repeats_introns_table_20$gene_id,intersect=overlaps_repeats_introns_table_20$intersect_width))
lncrna_overlap_regions_20<-unique(data.table(ID=overlaps_repeats_lncrna_bases_table_regions$ID,gene_id=overlaps_repeats_lncrna_bases_table_regions$gene_id,intersect=overlaps_repeats_lncrna_bases_table_regions$intersect_width))


lncrna_overlap_exons_20[,category:="exon"]
lncrna_overlap_introns_20[,category:="intron"]
lncrna_overlap_regions_20[,category:="five_prime"]

lncrna_repeat_overlap_category_20<-rbind(lncrna_overlap_exons_20,lncrna_overlap_introns_20,lncrna_overlap_regions_20)



lncrna_repeat_overlap_category_20_class<-merge(lncrna_repeat_overlap_category_20,lncRNA_classification_table_with_gene_id,by="gene_id")
lncrna_repeat_overlap_category_20_class_merged<-merge(lncrna_repeat_overlap_category_20_class,unique(overlaps_repeats_introns_table_20[,.SD,.SDcols=c("ID","repeat_class")]),by="ID")
lncrna_repeat_overlap_category_20_class_merged[,more:=ifelse(nrow(.SD)>1,T,F),by=c("ID","gene_id")]
lncrna_repeat_overlap_category_20_class_merged[,max_width_element:=ifelse(.SD==max(.SD),T,F),.SDcols="intersect",by=c("gene_id","ID")]
lncrna_repeat_overlap_category_20_class_merged[max_width_element==T,more:=F]
lncrna_repeat_overlap_category_20_class_merged<-lncrna_repeat_overlap_category_20_class_merged[more==F]
write.table(lncrna_repeat_overlap_category_20_class_merged,file="lncrna_repeat_overlap_category_20_class_merged",row.names=F,quote=F)
#Protein

overlaps_repeat_lncrna_num_table_20<-lncrna_repeat_overlap_category_20_class_merged[,.(number_of_bases=sum(intersect)),by=c("repeat_class","category","class")]




overlaps_repeats_proteins_bases<-findOverlaps(repeatmasker_esu_granges,protein_final_reduced_together_1kb_final,ignore.strand=T,type="any")
overlaps_repeats_proteins_bases_table<-cbind(as.data.table(repeatmasker_esu[queryHits(overlaps_repeats_proteins_bases)]),as.data.table(protein_coding_genes_granges[subjectHits(overlaps_repeats_proteins_bases)]))
setnames(overlaps_repeats_proteins_bases_table,c("SW_score","perc.sub","perc.del","perc.ins","query_name","start1","end1","base_past_end","strand_ignore","repeat_name","repeat_start","repeat_end","repeat_left_bases","ID","V16","repeat_class","repeat_family","species","size","library","group","ltrret_ID","type","element_type","seqnames2","start2","end2","width2","strand2","source","feature","score","strand_ignore2","frame","attribute"))
inter_width_vector<-c()
for (i in 1:nrow(overlaps_repeats_proteins_bases_table)) {
  inter_width_vector<-c(inter_width_vector,Overlap(c(overlaps_repeats_proteins_bases_table$start1[i],overlaps_repeats_proteins_bases_table$end1[i]),c(overlaps_repeats_proteins_bases_table$start2[i],overlaps_repeats_proteins_bases_table$end2[i])))
}
overlaps_repeats_proteins_bases_table[,intersect_width:=inter_width_vector]

#repeats exons overlap

overlaps_repeats_protein_exons_20_prot<-findOverlaps(repeatmasker_esu_granges,protein_exons,ignore.strand=T)
overlaps_repeats_protein_exons_table_20_prot<-cbind(as.data.table(repeatmasker_esu_granges[queryHits(overlaps_repeats_protein_exons_20_prot)]),as.data.table(protein_exons[subjectHits(overlaps_repeats_protein_exons_20_prot)]))
setnames(overlaps_repeats_protein_exons_table_20_prot,c("seqnames1","start1","end1","width1","strand1","SW_score","perc.sub","perc.del","perc.ins","base_past_end","strand_ignore","repeat_name","repeat_start","repeat_end","repeat_left_bases","ID","V16","repeat_class","repeat_family","species","size","library","group","ltrret_ID","type","elementtype","seqnames2","start2","end2","width2","strand2","source","feature","score","strand_ignore2","frame","attribute"))
overlaps_repeats_protein_exons_table_20_prot[,attribute:=str_extract(attribute,".*(?=\\.t\\d*)")]
inter_width_vector<-c()
for (i in 1:nrow(overlaps_repeats_protein_exons_table_20_prot)) {
  inter_width_vector<-c(inter_width_vector,Overlap(c(overlaps_repeats_protein_exons_table_20_prot$start1[i],overlaps_repeats_protein_exons_table_20_prot$end1[i]),c(overlaps_repeats_protein_exons_table_20_prot$start2[i],overlaps_repeats_protein_exons_table_20_prot$end2[i])))
}
overlaps_repeats_protein_exons_table_20_prot[,intersect_width:=inter_width_vector]
#repeats introns overlap
overlaps_repeats_protein_introns_20_prot<-findOverlaps(repeatmasker_esu_granges,protein_introns,ignore.strand=T)
overlaps_repeats_protein_introns_table_20_prot<-cbind(as.data.table(repeatmasker_esu_granges[queryHits(overlaps_repeats_protein_introns_20_prot)]),as.data.table(protein_introns[subjectHits(overlaps_repeats_protein_introns_20_prot)]))
setnames(overlaps_repeats_protein_introns_table_20_prot,c("seqnames1","start1","end1","width1","strand1","SW_score","perc.sub","perc.del","perc.ins","base_past_end","strand_ignore","repeat_name","repeat_start","repeat_end","repeat_left_bases","ID","V16","repeat_class","repeat_family","species","size","library","group","ltrret_ID","type","elementtype","seqnames2","start2","end2","width2","strand2","attribute"))
overlaps_repeats_protein_introns_table_20_prot[,attribute:=str_extract(attribute,".*(?=\\.t\\d*)")]
inter_width_vector<-c()
for (i in 1:nrow(overlaps_repeats_protein_introns_table_20_prot)) {
  inter_width_vector<-c(inter_width_vector,Overlap(c(overlaps_repeats_protein_introns_table_20_prot$start1[i],overlaps_repeats_protein_introns_table_20_prot$end1[i]),c(overlaps_repeats_protein_introns_table_20_prot$start2[i],overlaps_repeats_protein_introns_table_20_prot$end2[i])))
}
overlaps_repeats_protein_introns_table_20_prot[,intersect_width:=inter_width_vector]

#repeat regions overlap (5' 1 kb)
overlaps_repeats_protein_regions_20_prot_regions<-findOverlaps(repeatmasker_esu_granges,regions_prot_together,ignore.strand=T,type="any")
overlaps_repeats_protein_regions_20_prot_regions_table<-cbind(as.data.table(repeatmasker_esu[queryHits(overlaps_repeats_protein_regions_20_prot_regions)]),as.data.table(regions_prot_together[subjectHits(overlaps_repeats_protein_regions_20_prot_regions)]))
setnames(overlaps_repeats_protein_regions_20_prot_regions_table,c("SW_score","perc.sub","perc.del","perc.ins","query_name","start1","end1","base_past_end","strand_ignore1","repeat_name","repeat_start","repeat_end","repeat_left_bases","ID","V16","repeat_class","repeat_family","species","size","library","group","ltrret_ID","type","elementtype","seqnames2","start2","end2","width2","strand2","source","feature","score","strand_ignore2","frame","attribute"))
inter_width_vector<-c()
for (i in 1:nrow(overlaps_repeats_protein_regions_20_prot_regions_table)) {
  inter_width_vector<-c(inter_width_vector,Overlap(c(overlaps_repeats_protein_regions_20_prot_regions_table$start1[i],overlaps_repeats_protein_regions_20_prot_regions_table$end1[i]),c(overlaps_repeats_protein_regions_20_prot_regions_table$start2[i],overlaps_repeats_protein_regions_20_prot_regions_table$end2[i])))
}
overlaps_repeats_protein_regions_20_prot_regions_table[,intersect_width:=inter_width_vector]

protein_overlap_exons_20_prot<-unique(data.table(ID=overlaps_repeats_protein_exons_table_20_prot$ID,attribute=overlaps_repeats_protein_exons_table_20_prot$attribute,intersect=overlaps_repeats_protein_exons_table_20_prot$intersect_width))
protein_overlap_introns_20_prot<-unique(data.table(ID=overlaps_repeats_protein_introns_table_20_prot$ID,attribute=overlaps_repeats_protein_introns_table_20_prot$attribute,intersect=overlaps_repeats_protein_introns_table_20_prot$intersect_width))
protein_overlap_regions_20<-unique(data.table(ID=overlaps_repeats_protein_regions_20_prot_regions_table$ID,attribute=overlaps_repeats_protein_regions_20_prot_regions_table$attribute,intersect=overlaps_repeats_protein_regions_20_prot_regions_table$intersect_width))

protein_overlap_exons_20_prot[,category:="exon"]
protein_overlap_introns_20_prot[,category:="intron"]
protein_overlap_regions_20[,category:="five_prime"]
protein_repeat_overlap_category_20<-rbind(protein_overlap_exons_20_prot,protein_overlap_introns_20_prot,protein_overlap_regions_20)

protein_repeat_overlap_category_20[,class:="protein_coding_genes"]


protein_repeat_overlap_category_20_prot_merged<-merge(protein_repeat_overlap_category_20,unique(overlaps_repeats_proteins_bases_table[,.SD,.SDcols=c("ID","repeat_class")]),by="ID")

protein_repeat_overlap_category_20_prot_merged[,more:=ifelse(nrow(.SD)>1,T,F),by=c("ID","attribute")]
protein_repeat_overlap_category_20_prot_merged[,max_width_element:=ifelse(.SD==max(.SD),T,F),.SDcols="intersect",by=c("attribute","ID")]
protein_repeat_overlap_category_20_prot_merged[max_width_element==T,more:=F]
protein_repeat_overlap_category_20_prot_merged<-protein_repeat_overlap_category_20_prot_merged[more==F]
overlaps_repeat_protein_num_table_20<-protein_repeat_overlap_category_20_prot_merged[,.(number_of_bases=sum(intersect)),by=c("class","category","repeat_class")]

overlaps_repeat_together_num_table_20<-rbind(overlaps_repeat_lncrna_num_table_20,overlaps_repeat_protein_num_table_20)

overlaps_repeat_together_num_table_20[,category:=factor(category,levels = c("intron","exon","five_prime"))]
overlaps_repeat_together_num_table_normalized_sum_bases_20<-copy(overlaps_repeat_together_num_table_20)

overlaps_repeat_together_num_table_20[,class:=factor(class,levels = c("intergenic_lncRNA","intronic_lncRNA","overlapping_lncRNA","protein_coding_genes"))]
overlaps_repeat_together_num_table_20[class=="intergenic_lncRNA",class:="intergenske lncRNA"]
overlaps_repeat_together_num_table_20[class=="intronic_lncRNA",class:="intronske lncRNA"]
overlaps_repeat_together_num_table_20[class=="overlapping_lncRNA",class:="preklapajuæe lncRNA"]
overlaps_repeat_together_num_table_20[class=="protein_coding_genes",class:="protein kodirajuæi geni"]

repeats_plot_20<-ggplot(data=overlaps_repeat_together_num_table_20,aes(x=reorder(repeat_class,-number_of_bases),y=number_of_bases,fill=category)) +
  geom_bar(position=position_dodge(preserve = "single"),stat = "identity", width=0.7,col="black") +
  theme_bw() +
  facet_wrap(.~class) +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle = 50, hjust = 1,size=5)) +
  xlab("Obitelj ponavljanja") +
  ylab("Broj baza koje se preklapaju s ponavljanjima") +
  labs(fill = "lncRNA element",size=5) +
  theme(axis.title=element_text(size=8)) +
  theme(legend.title = element_text(size=8,face="bold")) +
  scale_fill_got(discrete = T,option="Martell",name="lncRNA element",labels=c("intron","egzon","5' regija"))
repeats_plot_20
ggsave(repeats_plot_20,file="repeats_plot_20.jpg")

#length sum normalization (
regions_lncrna_together_dt<-as.data.table(regions_lncrna_together)
regions_lncrna_together_dt_class<-merge(regions_lncrna_together_dt,lncRNA_classification_table_with_gene_id,by="gene_id")
regions_sum_table<-regions_lncrna_together_dt_class[,.(sum_el=sum(width)),by=class]

overlapping_exons_sum<-exons_sum_table[class=="overlapping_lncRNA"]$sum_el
overlapping_introns_sum<-introns_sum_table[class=="overlapping_lncRNA"]$sum_el
overlapping_regions_sum<-regions_sum_table[class=="overlapping_lncRNA"]$sum_el
intronic_exons_sum<-exons_sum_table[class=="intronic_lncRNA"]$sum_el
intronic_introns_sum<-introns_sum_table[class=="intronic_lncRNA"]$sum_el
intronic_regions_sum<-regions_sum_table[class=="intronic_lncRNA"]$sum_el
intergenic_exons_sum<-exons_sum_table[class=="intergenic_lncRNA"]$sum_el
intergenic_introns_sum<-introns_sum_table[class=="intergenic_lncRNA"]$sum_el
intergenic_regions_sum<-regions_sum_table[class=="intergenic_lncRNA"]$sum_el
protein_exons_sum<-sum(width(unique(protein_exons)))
protein_introns_sum<-sum(width(unique(protein_introns)))
protein_regions_sum<-norm_width_1bb_prot
location_sum<-data.table(class=rep(c("overlapping_lncRNA","intronic_lncRNA","intergenic_lncRNA","protein_coding_genses"),3),category=rep(c("exon","intron","five_prime"),4),sum_by_location=c(overlapping_exons_sum,overlapping_introns_sum,overlapping_regions_sum,intronic_exons_sum,intronic_introns_sum,intronic_regions_sum,intergenic_exons_sum,intergenic_introns_sum,intergenic_regions_sum,protein_exons_sum,protein_introns_sum,protein_regions_sum))


class_sum<-repeatmasker_esu[,.(total_sum=sum(end-start)),by="repeat_class"]
overlaps_repeat_together_num_table_normalized_sum_bases_20<-merge(overlaps_repeat_together_num_table_normalized_sum_bases_20,class_sum,by="repeat_class")


overlaps_repeat_together_num_table_normalized_sum_bases_20[,number_of_bases:=number_of_bases/total_sum]
overlaps_repeat_together_num_table_normalized_sum_bases_20[,number_of_bases:=as.double(number_of_bases)]
overlaps_repeat_together_num_table_normalized_sum_bases_20<-merge(overlaps_repeat_together_num_table_normalized_sum_bases_20,location_sum,by=c("class","category"))
overlaps_repeat_together_num_table_normalized_sum_bases_20[,number_of_bases:=number_of_bases/location_sum]
overlaps_repeat_together_num_table_normalized_sum_bases_20<-overlaps_repeat_together_num_table_normalized_sum_bases_20[is.na(repeat_class)==F]
overlaps_repeat_together_num_table_normalized_sum_bases_20[,class:=factor(class,levels = c("intergenic_lncRNA","intronic_lncRNA","overlapping_lncRNA","protein_coding_genes"))]
overlaps_repeat_together_num_table_normalized_sum_bases_20[class=="intergenic_lncRNA",class:="Intergenske lncRNA"]
overlaps_repeat_together_num_table_normalized_sum_bases_20[class=="intronic_lncRNA",class:="Intronske lncRNA"]
overlaps_repeat_together_num_table_normalized_sum_bases_20[class=="overlapping_lncRNA",class:="Preklapajuæe lncRNA"]
overlaps_repeat_together_num_table_normalized_sum_bases_20[class=="protein_coding_genes",class:="Protein-kodirajuæi geni"]
overlaps_repeat_together_num_table_normalized_sum_bases_20[,repeat_class:=factor(repeat_class,levels = c("Unknown","LINE","LTR","DNA"))]
overlaps_repeat_together_num_table_normalized_sum_bases_20[,category:=factor(category,levels = c("five_prime","exon","intron"))]
repeats_plot_normalized_sum_bases_20<-ggplot(data=overlaps_repeat_together_num_table_normalized_sum_bases_20,aes(x=repeat_class,y=number_of_bases,fill=category)) +
  geom_bar(position=position_dodge(preserve = "single"),stat = "identity", width=0.7,col="black",alpha=0.8) +
  theme_bw() +
  facet_grid(row=vars(class)) +
  theme(legend.position = "none") + 
  #theme(axis.text.x = element_text(angle = 50, hjust = 1,size=5)) +
  xlab("Razred transpozona") +
  ylab("Normalizirana ukupna duljina") +
  labs(fill = "lncRNA element",size=5) +
  theme(axis.title=element_text(size=10)) +
  theme(axis.text=element_text(size=9))+
  theme(legend.title = element_text(size=10,face="bold")) +
  theme(legend.text = element_text(size=10)) +
  scale_fill_got(discrete = T,option="Martell",name="Element",labels=c("5' regija","egzon","intron")) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(strip.text.y = element_text(size = 7))

ggsave(repeats_plot_normalized_sum_bases_20,file="repeats_plot_normalized_sum_bases_20.jpg",width=7,height = 6)


#Summing ovwerlap bases - comparison between lncRNA and proteins total


overlaps_repeat_lncrna_num_table_20_new<-lncrna_repeat_overlap_category_20_class_merged[,.(number_of_bases=sum(intersect)),by=c("class","category")]
overlaps_repeat_protein_num_table_20_new<-protein_repeat_overlap_category_20_prot_merged[,.(number_of_bases=sum(intersect)),by=c("class","category")]
overlaps_repeat_together_num_table_normalized_sum_bases_20_2<-rbind(overlaps_repeat_lncrna_num_table_20_new,overlaps_repeat_protein_num_table_20_new)
overlaps_repeat_together_num_table_normalized_sum_bases_20_2[,number_of_bases:=as.double(number_of_bases)]
overlaps_repeat_together_num_table_normalized_sum_bases_20_2<-merge(overlaps_repeat_together_num_table_normalized_sum_bases_20_2,location_sum,by=c("class","category"))

overlaps_repeat_together_num_table_normalized_sum_bases_20_2[,class:=factor(class,levels = c("intergenic_lncRNA","intronic_lncRNA","overlapping_lncRNA","protein_coding_genes"))]
overlaps_repeat_together_num_table_normalized_sum_bases_20_2[class=="intergenic_lncRNA",class:="Intergenske lncRNA"]
overlaps_repeat_together_num_table_normalized_sum_bases_20_2[class=="intronic_lncRNA",class:="Intronske lncRNA"]
overlaps_repeat_together_num_table_normalized_sum_bases_20_2[class=="overlapping_lncRNA",class:="Preklapajuæe lncRNA"]
overlaps_repeat_together_num_table_normalized_sum_bases_20_2[class=="protein_coding_genes",class:="Protein-kodirajuæi geni"]

overlaps_repeat_together_num_table_normalized_sum_bases_20_2[,category:=factor(category,levels = c("five_prime","exon","intron"))]
repeats_plot_normalized_sum_bases_20_2<-ggplot(data=overlaps_repeat_together_num_table_normalized_sum_bases_20_2,aes(x=class,y=number_of_bases,fill=category)) +
  geom_bar(position="dodge",stat = "identity",col="black",alpha=0.8) +
  theme_bw() +
  theme(legend.position = "top") + 
  #theme(axis.text.x = element_text(angle = 50, hjust = 1,size=5)) +
  xlab("Klasa gena") +
  ylab("Normalizirana ukupna duljina") +
  labs(fill = "lncRNA element",size=5) +
  theme(axis.title=element_text(size=11)) +
  theme(axis.text=element_text(size=10.3))+
  theme(legend.title = element_text(size=10,face="bold")) +
  theme(legend.text = element_text(size=10)) +
  scale_fill_got(discrete = T,option="Martell",name="Element",labels=c("5' regija","egzon","intron")) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) 

ggsave(repeats_plot_normalized_sum_bases_20_2,file="repeats_plot_normalized_sum_bases_20_2.jpg",width=7.6,height = 3.5)