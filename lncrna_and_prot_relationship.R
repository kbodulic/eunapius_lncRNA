args <- commandArgs(trailingOnly = TRUE)  
library(data.table)
library(IRanges)
library(GenomicRanges)
library(ggplot2)
library(harrypotter)
setwd(".")


#analysis of the relationship between lncRNA and protein-coding genes - distribution of lengths and numbers of overlapping, intronic and intergenic lncrna, intergenic lncrna distances to closest protein-coding genes, writing protein-coding genes with an overlap with a lncrna-coding gene and protein coding genes 1 kb or closer to the closest lncrna-coding gene (for GO analysis)
#Arguments: 1 - protein coding annotation, 2 - lncRNA inron table, 3 - rnaspades paf of all filtered transcripts, trinity paf od all filtered transcripts, filtered exons, isoform-gene overlapping table for rnaspades, isoform-gene overlapping table for trinity, lncRNA gene list, protein longest isoform list



protein_coding_annotation_clean<-fread(args[1])
filtered_together_paf_union<-fread(args[2])
introns_table<-fread(args[3])
rnaspades_paf_filtered_union<-fread(args[4])
trinity_paf_filtered_union<-fread(args[5])
filtered_together_paf_union_exons<-fread(args[6])
lncrnas_final_reduced_together_dt<-fread(args[7])
protein_coding_annotation_longest_isoforms_copy<-fread(args[8])


#taking all protein exons and introns


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



#taking the consensus paf 

filtered_together_paf_union_granges<-makeGRangesFromDataFrame(filtered_together_paf_union,seqnames.field = "target_name",start.field = "target_start",end.field = "target_end",ignore.strand=T,keep.extra.columns = T)



#intronic lncRNA - found inside protein-coding introns
overlaps_lncRNA_prot_within<-findOverlaps(filtered_together_paf_union_granges,protein_coding_genes_granges,ignore.strand=T,type = "within")
overlaps_lncRNA_prot_within_table<-cbind(as.data.table(filtered_together_paf_union_granges[queryHits(overlaps_lncRNA_prot_within)]),as.data.table(protein_coding_genes_granges[subjectHits(overlaps_lncRNA_prot_within)]))
setnames(overlaps_lncRNA_prot_within_table,c("seqnames1","start1","end1","width1","strand1","query_name","query_length","query_start","query_end","strand_ignore2","target_length","matches","alignment_block","mapping_quality","perc_matches","gene","seqnames2","start2","end2","width2","strand2","source","feature","score","strand_ignore2","frame","attribute"))

intronic_lncrna<-unique(overlaps_lncRNA_prot_within_table$query_name)
#overlapping lncRNA - hold protein-coding genes in their introns 




overlaps_prot_lncrna_within<-findOverlaps(protein_coding_genes_granges,filtered_together_paf_union_granges,ignore.strand=T,type="within")

overlaps_prot_lncrna_within_table<-cbind(as.data.table(protein_coding_genes_granges[queryHits(overlaps_prot_lncrna_within)]),as.data.table(filtered_together_paf_union_granges[subjectHits(overlaps_prot_lncrna_within)]))
setnames(overlaps_prot_lncrna_within_table,c("seqnames1","start1","end1","width1","strand1","source","feature","score","strand_ignore1","frame","attribute","seqnames2","start2","end2","width2","strand2","query_name","query_length","query_start","query_end","strand_ignore2","target_length","matches","alignment_block","mapping_quality","perc_matches","gene"))
overlapping_lncrna<-unique(overlaps_prot_lncrna_within_table$query_name)



intronic_lncrna<-intronic_lncrna[intronic_lncrna%in%intersect(overlapping_lncrna,intronic_lncrna)==F]


#intergenic lncRNA - no overlaps with protein-coding genes


lncrna_prot_transcripts_overlaps<-findOverlaps(filtered_together_paf_union_granges,protein_coding_genes_granges,ignore.strand=T)
lncrna_prot_transcripts_overlaps_table<-cbind(as.data.table(filtered_together_paf_union_granges[queryHits(lncrna_prot_transcripts_overlaps)]),as.data.table(protein_coding_genes_granges[subjectHits(lncrna_prot_transcripts_overlaps)]))
setnames(lncrna_prot_transcripts_overlaps_table,c("seqnames1","start1","end1","width1","strand1","query_name","query_length","query_start","query_end","strand_ignore","target_length","matches","aligmnent_block","mapping_quality","perc_matches","gene","seqnames2","start2","end2","width2","strand2","source","feature","score","strand_ignore2","frame","attribute"))

intergenic_lncrna<-unique(filtered_together_paf_union[query_name %in% lncrna_prot_transcripts_overlaps_table$query_name==F]$query_name)





#checking for proteins whose biggest part lies in an overlapping lncRNA, while the whole last exon lies outside of the lncRNA - classifying kncRNA with those proteins as overlapping
#intron exon overlap
rnaspades_paf_filtered_union<-rnaspades_paf_filtered_union[rnaspades_paf_filtered_union$query_name%in%rnaspades_longest_iso]
trinity_paf_filtered_union<-trinity_paf_filtered_union[trinity_paf_filtered_union$query_name%in%trinity_longest_iso]


filtered_together_paf_union<-rbind(rnaspades_paf_filtered_union,trinity_paf_filtered_union)

rnaspades_paf_filtered_union_granges<-makeGRangesFromDataFrame(rnaspades_paf_filtered_union,seqnames.field = "target_name",start.field = "target_start",end.field = "target_end",ignore.strand=T,keep.extra.columns = T)
trinity_paf_filtered_union_granges<-makeGRangesFromDataFrame(trinity_paf_filtered_union,seqnames.field = "target_name",start.field = "target_start",end.field = "target_end",ignore.strand=T,keep.extra.columns = T)
filtered_together_paf_union_granges<-c(rnaspades_paf_filtered_union_granges,trinity_paf_filtered_union_granges)



int_ex_overlap_candidates_intronic<-GenomicRanges::setdiff(filtered_together_paf_union_granges$query_name,intronic_lncrna)
int_ex_overlap_candidates_intronic_overlapping<-GenomicRanges::setdiff(int_ex_overlap_candidates_intronic,overlapping_lncrna)
int_ex_overlap_candidates_intronic_overlapping_intergenic<-GenomicRanges::setdiff(int_ex_overlap_candidates_intronic_overlapping,intergenic_lncrna)

int_ex_overlap_candidates_exons<-makeGRangesFromDataFrame(filtered_together_paf_union_exons[filtered_together_paf_union_exons$transcript %in% int_ex_overlap_candidates_intronic_overlapping_intergenic],seqnames.field = "chromosome",ignore.strand=T,keep.extra.columns = T)

overlaps_candidates_exons_prot_introns<-findOverlaps(int_ex_overlap_candidates_exons,protein_introns,ignore.strand=T)
overlaps_candidates_exons_prot_introns_table<-cbind(as.data.table(int_ex_overlap_candidates_exons[queryHits(overlaps_candidates_exons_prot_introns)]),as.data.table(protein_introns[subjectHits(overlaps_candidates_exons_prot_introns)]))
setnames(overlaps_candidates_exons_prot_introns_table,c("seqnames1","start1","end1","width1","strand1","transcript","seqnames2","start2","end2","width2","strand2","transcript_int"))
overlapping_lncrna<-c(overlapping_lncrna,overlaps_candidates_exons_prot_introns_table$transcript)
intergenic_lncrna<-c(intergenic_lncrna,GenomicRanges::setdiff(int_ex_overlap_candidates_intronic_overlapping_intergenic,unique(overlaps_candidates_exons_prot_introns_table$transcript)))


#transcript number per class

lncRNA_classification_table<-unique(data.table(transcript=c(intronic_lncrna,overlapping_lncrna,intergenic_lncrna),class=c(rep("intronic_lncRNA",length(intronic_lncrna)),rep("overlapping_lncRNA",length(overlapping_lncrna)),rep("intergenic_lncRNA",length(intergenic_lncrna)))))
write.table(lncRNA_classification_table,file="lncRNA_classification_table.txt",row.names=F,quote=F)
number_of_lncrna_classes<-lncRNA_classification_table[,.(number_of_lncRNA=.N),by=class]
number_of_lncrna_classes_prikaz<-copy(number_of_lncrna_classes)
setnames(number_of_lncrna_classes_prikaz,c("Klasa","Broj"))
write.table(number_of_lncrna_classes,file="number_of_lncrna_classes.txt",row.names = F,quote = F)




#transcript lengths by class
lncRNA_classification_table_name<-copy(lncRNA_classification_table)
setnames(lncRNA_classification_table_name,c("query_name","class"))
lncrna_paf_with_classification<-merge(filtered_together_paf_union,lncRNA_classification_table_name,by="query_name")

class_transcripts_table<-data.table(element_length=c(lncrna_paf_with_classification[class=="intergenic_lncRNA"]$query_length,lncrna_paf_with_classification[class=="intronic_lncRNA"]$query_length,lncrna_paf_with_classification[class=="overlapping_lncRNA"]$query_length),class=c(rep("intergenic_lncRNA",length(lncrna_paf_with_classification[class=="intergenic_lncRNA"]$query_length)),rep("intronic_lncRNA",length(lncrna_paf_with_classification[class=="intronic_lncRNA"]$query_length)),rep("overlapping_lncRNA",length(lncrna_paf_with_classification[class=="overlapping_lncRNA"]$query_length))))

class_transcripts_table[,class:=factor(class,levels=c("intergenic_lncRNA","intronic_lncRNA","overlapping_lncRNA"))]
class_transcripts_table[,feature:="transkripti"]


#class exons

filtered_together_paf_union_exons_class<-merge(filtered_together_paf_union_exons,lncRNA_classification_table,by="transcript")

class_exons_table<-filtered_together_paf_union_exons_class[,.SD,.SDcols=c("width","class")]
class_exons_table[,feature:="egzoni"]
setnames(class_exons_table,c("element_length","class","feature"))

#class introns


filtered_together_paf_union_introns_class<-merge(introns_table,lncRNA_classification_table,by="transcript")
class_introns_table<-filtered_together_paf_union_introns_class[,.SD,.SDcols=c("width","class")]
class_introns_table[,feature:="introni"]
setnames(class_introns_table,c("element_length","class","feature"))
#class genes



setnames(lncRNA_classification_table,c("transcript_name","class"))
overlaps_table_rnaspades_with_class<-merge(isoform_overlap_table_rnaspades,lncRNA_classification_table,by="transcript_name")
overlaps_table_trinity_with_class<-merge(isoform_overlap_table_trinity,lncRNA_classification_table,by="transcript_name")

class_genes_rnaspades<-unique(overlaps_table_rnaspades_with_class[,.SD,.SDcols="gene_id"])
class_genes_trinity<-unique(overlaps_table_trinity_with_class[,.SD,.SDcols="gene_id"])

class_genes_table<-rbind(class_genes_rnaspades,class_genes_trinity)
class_genes_table<-class_genes_table[,.SD,.SDcols=c("width1","class")]
setnames(class_genes_table,c("element_length","class"))
class_genes_table[,feature:="geni"]

class_together_table<-rbind(class_transcripts_table,class_exons_table,class_introns_table,class_genes_table)
class_together_table<-class_together_table[,feature:=factor(feature,c("egzoni","introni","transkripti","geni"))]
class_together_table[,class:=factor(class,levels = c("intronic_lncRNA","intergenic_lncRNA","overlapping_lncRNA"))]
class_lengths_plot<-ggplot(data=class_together_table,aes(x=class,y=element_length,fill=class)) +
  geom_boxplot(alpha=0.8) +
  scale_y_log10() +
  theme_bw() +
  facet_grid(.~feature) +
  theme(legend.position = "top") +
  theme(legend.title = element_text(size=9,face="bold")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x =element_blank()) +
  scale_fill_hp(discrete = TRUE, option = "NewtScamander", name = "Kategorija",labels=c("Intronske lncRNA","Intergenske lncRNA","Preklapajuæe lncRNA")) +
  ylab("Duljina (log10 skalirano)") 


ggsave(class_lengths_plot,file="class_lengths_plot.jpg",width=6,height = 4)

#genesum of bases by class
class_genes_table
class_genes_table_sum_width<-class_genes_table[,.(bases_sum=sum(element_length)),by=class]



#intron sum

class_introns_table_sum_width<-class_introns_table[,.(intron_sum=sum(element_length)),by=class]

#exon sum


class_exons_table_sum_widthh<-class_exons_table[,.(exon_sum=sum(element_length)),by=class]

class_table_lengths<-merge(class_genes_table_sum_width,class_exons_table_sum_widthh,by="class")
class_table_lengths<-merge(class_table_lengths,class_introns_table_sum_width,by="class")

class_table_lengths[,intron_percentage:=100*intron_sum/bases_sum]

kruskal.test(x = c(class_exons_table[class=="intergenic_lncRNA"]$element_length,class_exons_table[class=="intronic_lncRNA"]$element_length,class_exons_table[class=="overlapping_lncRNA"]$element_length),g=c(rep(1,length(class_exons_table[class=="intergenic_lncRNA"]$element_length)),rep(2,length(class_exons_table[class=="intronic_lncRNA"]$element_length)),rep(3,length(class_exons_table[class=="overlapping_lncRNA"]$element_length))))

kruskal.test(x = c(class_introns_table[class=="intergenic_lncRNA"]$element_length,class_introns_table[class=="intronic_lncRNA"]$element_length,class_introns_table[class=="overlapping_lncRNA"]$element_length),g=c(rep(1,length(class_introns_table[class=="intergenic_lncRNA"]$element_length)),rep(2,length(class_introns_table[class=="intronic_lncRNA"]$element_length)),rep(3,length(class_introns_table[class=="overlapping_lncRNA"]$element_length))))

kruskal.test(x = c(class_genes_table[class=="intergenic_lncRNA"]$element_length,class_genes_table[class=="intronic_lncRNA"]$element_length,class_genes_table[class=="overlapping_lncRNA"]$element_length),g=c(rep(1,length(class_genes_table[class=="intergenic_lncRNA"]$element_length)),rep(2,length(class_genes_table[class=="intronic_lncRNA"]$element_length)),rep(3,length(class_genes_table[class=="overlapping_lncRNA"]$element_length))))

kruskal.test(x = c(class_transcripts_table[class=="intergenic_lncRNA"]$element_length,class_transcripts_table[class=="intronic_lncRNA"]$element_length,class_transcripts_table[class=="overlapping_lncRNA"]$element_length),g=c(rep(1,length(class_transcripts_table[class=="intergenic_lncRNA"]$element_length)),rep(2,length(class_transcripts_table[class=="intronic_lncRNA"]$element_length)),rep(3,length(class_transcripts_table[class=="overlapping_lncRNA"]$element_length))))


#intergenic lncRNA - distance to genes


#filtering for intergenic lncRNA - taking whole genes
class_gene_together<-unique(rbind(class_genes_rnaspades,class_genes_trinity))
setnames(lncrnas_final_reduced_together_dt,c("seqnames1","start1","end1","width1","strand1","gene","gene_id"))
lncrnas_final_reduced_together_dt_merged<-merge(lncrnas_final_reduced_together_dt,class_gene_together,by=c("seqnames1","start1","end1","width1"),all.x=T)
intergenic_lncrna_dt<-lncrnas_final_reduced_together_dt_merged[class=="intergenic_lncRNA"]
intergenic_lncrna_granges<-makeGRangesFromDataFrame(intergenic_lncrna_dt,seqnames.field = "seqnames1",start.field = "start1",end.field = "end1",keep.extra.columns = T)
intergenic_nearest_gene_nearest<-nearest(intergenic_lncrna_granges,protein_coding_genes_granges)
prot_cod_intergenic_nearest<-protein_coding_genes[intergenic_nearest_gene_nearest]


intergenic_prot_nearest_table<-data.table(seqnames1=as.character(seqnames(intergenic_lncrna_granges)),start1=start(intergenic_lncrna_granges),end1=end(intergenic_lncrna_granges),lncrna=intergenic_lncrna_granges$query_name,seqnames2=prot_cod_intergenic_nearest$seqname,start2=prot_cod_intergenic_nearest$start,end2=prot_cod_intergenic_nearest$end,prot_cod_gene=prot_cod_intergenic_nearest$attribute)

inter_on_contigs_without_genes<-which(is.na(intergenic_prot_nearest_table$seqname2))


intergenic_prot_nearest_table<-intergenic_prot_nearest_table[is.na(seqnames2)==F]
intergenic_prot_nearest_table[,distance:=ifelse(start1>start2,start1-end2-1,end1-start2+1)]

inter_prot_cod_distance_plot<-ggplot(data=intergenic_prot_nearest_table,aes(x=distance)) +
  geom_histogram(color="black",fill="skyblue",binwidth = 250) +
  theme_bw()  +
  xlab("Udaljenost lncRNA-kodirajuæeg gena od najbližeg protein-kodirajuæeg gena") +
  ylab("Broj lncRNA-kodirajuæih gena") +
  coord_cartesian(xlim=c(-10000,10000)) +
  theme(axis.title.x   = element_text(size=8.5)) +
  theme(axis.title.y   = element_text(size=8.5)) 



ggsave(inter_prot_cod_distance_plot,file="inter_prot_cod_distance_plot.jpg",width = 5,height = 3.8)



#overlapping proteins vs non-overlapping proteins - writing protein-coding genes with an overlap with a lncrna-coding gene and protein coding genes 1 kb or closer to the closest lncrna-coding gene (for GO analysis)

#intronic


intronic_prots<-unique(overlaps_lncRNA_prot_within_table[query_name%in%intersect(overlapping_lncrna,intronic_lncrna)==F])$attribute

#overlapping

overlapping_prots<-unique(overlaps_prot_lncrna_within_table$attribute)


#intergenic - 1 kb

intergenic_prots<-unique(intergenic_prot_nearest_table[distance<=1000 & distance >= -1000]$prot_cod_gene)

category_prots<-c(intronic_prots,overlapping_prots,intergenic_prots)
no_category_prots<-GenomicRanges::setdiff(protein_coding_genes$attribute,category_prots)

prot_classification_table<-data.table(gene=c(category_prots,no_category_prots),lncrna_class=c(rep("intronic",length(intronic_prots)),rep("overlapping",length(overlapping_prots)),rep("intergenic_1kb",length(intergenic_prots)),rep("no_lncRNA",length(no_category_prots))))
write.table(prot_classification_table,file="prot_classification_table.txt",row.names=F,quote=F)

protein_coding_annotation_longest_isoforms_copy[,gene:=str_extract(attribute,".*(?=\\.t\\d*)")]
prot_lncrna_overlaps_genes<-protein_coding_annotation_longest_isoforms_copy[gene%in%prot_classification_table[lncrna_class!="no_lncRNA"]$gene]$attribute

prot_lncrna_intronic_overlaps_genes<-protein_coding_annotation_longest_isoforms_copy[gene%in%prot_classification_table[lncrna_class=="intronic"]$gene]$attribute
prot_lncrna_overlapping_overlaps_genes<-protein_coding_annotation_longest_isoforms_copy[gene%in%prot_classification_table[lncrna_class=="overlapping"]$gene]$attribute
prot_lncrna_intergenic_1kb_overlaps_genes<-protein_coding_annotation_longest_isoforms_copy[gene%in%prot_classification_table[lncrna_class=="intergenic_1kb"]$gene]$attribute

write.table(prot_lncrna_overlaps_genes,file="prot_lncrna_overlaps_genes.txt",row.names = F,col.names = F,quote = F,sep = "\t")
write.table(prot_lncrna_intronic_overlaps_genes,file="prot_lncrna_intronic_overlaps_genes.txt",row.names = F,col.names = F,quote = F,sep = "\t")
write.table(prot_lncrna_overlapping_overlaps_genes,file="prot_lncrna_overlapping_overlaps_genes.txt",row.names = F,col.names = F,quote = F,sep = "\t")
write.table(prot_lncrna_intergenic_1kb_overlaps_genes,file="prot_lncrna_intergenic_1kb_overlaps_genes.txt",row.names = F,col.names = F,quote = F,sep = "\t")
