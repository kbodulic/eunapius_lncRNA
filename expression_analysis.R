args <- commandArgs(trailingOnly = TRUE)  
library(data.table)
library(stringr)
library(ggplot2)
library(gameofthrones)
library(harrypotter)
library(tidyr)
setwd(".")

#analyzing the expression - calculating FPKM for every gene, comparing day 1 to day10 expression levels, comparing the expression levels of lncRNAs with a transposone insertions in one of their elements and lncRNAs without transposone insertions, comparing expression levels of protein with and without a relationship to a lncRNA
#Arguments: 1 - counts table for day1 and day10 genes, 2 - lncRNA classification table, 3 - a table which maps longest lncRNA isoform to a gene_id, 4 - a table with transposone and lncRNA overlaps, 5 - paf of lncRNA longest isoforms, 6 - protein classification table based on the relationship with lncRNA genes
read_counts_table<-fread(args[1])
lncRNA_classification_table<-fread(args[2])
isoform_map_table_together<-fread(args[3])
lncrna_repeat_overlap_category_20_class_merged<-fread(args[4])
filtered_together_paf_union<-fread(args[5])
prot_classification_table<-fread(args[6])

#calculating FPKMs
read_counts_table<-read_counts_table[,.SD,.SDcols=c("Geneid","Length","sorted_noPCR_RNAseq_day1_to_Esu_geonme.bam","sorted_noPCR_RNAseq_day10_to_Esu_geonme.bam")]
setnames(read_counts_table,c("gene_id","gene_length","counts_day1","counts_day10"))
total_counts_table<-fread("total_mapped_counts.txt")
setnames(total_counts_table,"counts")
total_counts_table[,experiment:=c("day1","day10")]
read_counts_table<-read_counts_table[counts_day1 !=0]
read_counts_table<-read_counts_table[counts_day10!=0]
read_counts_table[,counts_scaled_milion_day1:=1e+06*counts_day1]
read_counts_table[,counts_scaled_milion_day10:=1e+06*counts_day10]
read_counts_table[,FPKM_day1:=counts_scaled_milion_day1/(total_counts_table[experiment=="day1"]$counts*1e-3*gene_length)]
read_counts_table[,FPKM_day10:=counts_scaled_milion_day10/(total_counts_table[experiment=="day10"]$counts*1e-3*gene_length)]


#sorting day1 FPKMs into quantiles

read_counts_table_lncrna<-read_counts_table[category=="lncRNA"]
read_counts_table_lncrna<-read_counts_table_lncrna[counts_day1 !=0]
read_counts_table_lncrna<-read_counts_table_lncrna[counts_day10!=0]
day1_fpkm<-read_counts_table_lncrna$FPKM_day1

day1_quantiles<-quantile(day1_fpkm,probs=seq(0,1,0.1))
day1_percentiles<-cut(day1_fpkm,breaks = day1_quantiles, labels = paste(seq(0,90,10),seq(10,100,10),sep="-"))

#day 1 - day10 comparison
read_counts_table_lncrna[,percentiles:=day1_percentiles]
read_counts_table_lncrna<-read_counts_table_lncrna[is.na(percentiles)==F]
read_counts_table_lncrna[,percentiles:=paste(percentiles,"%",sep=" ")]
read_counts_table_lncrna[,FPKM_day1_log:=log(FPKM_day1)]
read_counts_table_lncrna[,FPKM_day10_log:=log(FPKM_day10)]

read_counts_table_lncrna[,gene_id:=as.double(gene_id)]



read_counts_table_lncrna_class<-merge(read_counts_table_lncrna,lncRNA_classification_table_with_gene_id,by="gene_id")
lncRNA_classification_table_with_gene_id<-merge(lncRNA_classification_table,isoform_map_table_together)

read_counts_table_lncrna_class[,class:=factor(class,levels = c("intronic_lncRNA","intergenic_lncRNA","overlapping_lncRNA"))]
read_counts_table_lncrna_class[class=="intronic_lncRNA",class:="Intronske lncRNA"]
read_counts_table_lncrna_class[class=="intergenic_lncRNA",class:="Intergenske lncRNA"]
read_counts_table_lncrna_class[class=="overlapping_lncRNA",class:="Preklapajuæe lncRNA"]
fpkm_groups_all_lncrna_plot<-ggplot(data=read_counts_table_lncrna_class,aes(x=FPKM_day1_log,y=FPKM_day10_log,col=percentiles)) +
  geom_point() +
  theme_bw() +
  facet_grid(rows=vars(class)) +
  scale_color_discrete(name="Postotak\nekspresije") +
  theme(legend.title = element_text(size=9,face="bold")) +
  xlab("FPKM knjižnice RNA1 (log10 skalirano)") +
  ylab("FPKM knjižnice RNA10(log10 skalirano)") +
  theme(axis.title = element_text(size=8)) +
  theme(strip.text = element_text(size = 6.5))

 
ggsave(fpkm_groups_all_lncrna_plot,file="fpkm_groups_all_lncrna_plot.jpg",width=6,height = 5)

#class comparison
read_counts_table_longer_lncrna<-read_counts_table_longer[category=="lncRNA"]
read_counts_table_longer_lncrna[,gene_id:=as.double(gene_id)]
read_counts_table_longer_lncrna_class<-merge(read_counts_table_longer_lncrna,lncRNA_classification_table_with_gene_id,by="gene_id",all.x=T)

fpkm_lncrna_class_plot<-ggplot(data=read_counts_table_longer_lncrna_class,aes(x=class,y=FPKM,fill=class)) +

read_counts_table_longer_lncrna_class_day1<-read_counts_table_longer_lncrna_class[experiment=="day1"]
kruskal.test(x = c(read_counts_table_longer_lncrna_class_day1[class=="intergenic_lncRNA"]$FPKM,read_counts_table_longer_lncrna_class_day1[class=="intronic_lncRNA"]$FPKM,read_counts_table_longer_lncrna_class_day1[class=="overlapping_lncRNA"]$FPKM),g=c(rep(1,length(read_counts_table_longer_lncrna_class_day1[class=="intergenic_lncRNA"]$FPKM)),rep(2,length(read_counts_table_longer_lncrna_class_day1[class=="intronic_lncRNA"]$FPKM)),rep(3,length(read_counts_table_longer_lncrna_class_day1[class=="overlapping_lncRNA"]$FPKM))))

read_counts_table_longer_lncrna_class_day10<-read_counts_table_longer_lncrna_class[experiment=="day10"]
kruskal.test(x = c(read_counts_table_longer_lncrna_class_day10[class=="intergenic_lncRNA"]$FPKM,read_counts_table_longer_lncrna_class_day10[class=="intronic_lncRNA"]$FPKM,read_counts_table_longer_lncrna_class_day10[class=="overlapping_lncRNA"]$FPKM),g=c(rep(1,length(read_counts_table_longer_lncrna_class_day10[class=="intergenic_lncRNA"]$FPKM)),rep(2,length(read_counts_table_longer_lncrna_class_day10[class=="intronic_lncRNA"]$FPKM)),rep(3,length(read_counts_table_longer_lncrna_class_day10[class=="overlapping_lncRNA"]$FPKM))))
#expression and transposons

lncrna_repeat_overlap_category_class_merged_copy<-copy(lncrna_repeat_overlap_category_20_class_merged)
lncrna_repeat_overlap_category_class_merged_copy[,sum_intersect:=sum(intersect),by=c("gene_id","ID")]
lncrna_repeat_overlap_category_class_merged_copy[,max_sum:=ifelse(.SD==max(.SD),T,F),.SDcols="sum_intersect",by=gene_id]
lncrna_repeat_overlap_category_class_merged_copy<-lncrna_repeat_overlap_category_class_merged_copy[max_sum==T]


#Transposone insertions comparison 
setnames(lncrna_repeat_overlap_category_class_merged_copy,c("ID","gene_id","intersect","element","transcript","class","repeat_class","more","max_width_element","sum_intersect","max_sum"))
lncrna_repeat_overlap_category_class_merged_copy_merged<-merge(lncrna_repeat_overlap_category_class_merged_copy,read_counts_table_longer_lncrna[,.SD,.SDcols=c("gene_id","experiment","FPKM")],by="gene_id")
isoform_map_table_together_cp2<-copy(isoform_map_table_together)
setnames(isoform_map_table_together_cp2,c("query_name","gene_id"))
filtered_together_paf_union_copy<-copy(filtered_together_paf_union)
filtered_together_paf_union_copy<-merge(filtered_together_paf_union_copy,isoform_map_table_together_cp2)
lncrna_no_repeats<-filtered_together_paf_union_copy[gene_id%in%lncrna_repeat_overlap_category_class_merged_copy_merged$gene_id==F]

lncrna_no_repeats_dt<-data.table(gene_id=lncrna_no_repeats$gene_id)
lncrna_no_repeats_dt_merged<-merge(lncrna_no_repeats_dt,read_counts_table_longer_lncrna[,.SD,.SDcols=c("gene_id","experiment","FPKM")],by="gene_id")
lncrna_no_repeats_dt_merged[,repeat_class:="No_transposons"]
lncrna_repeat_overlap_category_class_merged_copy_merged_with_no_rep<-rbind(lncrna_repeat_overlap_category_class_merged_copy_merged,lncrna_no_repeats_dt_merged,fill=T)
lncrna_repeat_overlap_category_class_merged_copy_merged_with_no_rep[,repeat_class:=factor(repeat_class,levels = c("Unknown","DNA","LINE","LTR","RC","No_transposons"))]
lncrna_repeat_overlap_category_class_merged_copy_merged_with_no_rep[experiment=="day1",experiment:="Prvi dan"]
lncrna_repeat_overlap_category_class_merged_copy_merged_with_no_rep[experiment=="day10",experiment:="Deseti dan"]
lncrna_repeat_overlap_category_class_merged_copy_merged_with_no_rep[,experiment:=factor(experiment,levels = c("Prvi dan","Deseti dan"))]


#lncrna_repeat_overlap_category_class_merged_copy_merged_with_no_rep[,element:=factor(element,levels = levels_for_plot)]

options(scipen=1000)
fpkm_lncrna_repeats_plot<-ggplot(data=lncrna_repeat_overlap_category_class_merged_copy_merged_with_no_rep,aes(x=repeat_class,y=FPKM,fill=element)) +
  geom_boxplot(alpha=0.8) +
  theme_bw() +
  facet_wrap(.~experiment) +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  scale_fill_got(discrete = T,option = "martell",na.value="grey",labels=c("egzon","intron","neodreðeno"),direction = -1) +
  xlab("Razred transpozona") +
  #theme(axis.text.x=element_blank()) +
  theme(legend.position = "none") +
  scale_x_discrete(labels=c("Unknown","DNA","LINE","LTR","Bez transpozona"))



fpkm_lncrna_repeats_plot
ggsave(fpkm_lncrna_repeats_plot,file="fpkm_lncrna_repeats_plot.jpg",width = 6,height = 5)
lncrna_repeat_overlap_category_class_merged_copy_merged_with_no_rep[,log_FPKM:=log10(FPKM)]
lncrna_repeat_overlap_category_class_merged_copy_merged_with_no_rep_2<-lncrna_repeat_overlap_category_class_merged_copy_merged_with_no_rep[experiment=="day1"]
lncrna_repeat_overlap_category_class_merged_copy_merged_with_no_rep[,cate:=paste(element,repeat_class,sep=" ")]
anova_rep<-aov(log_FPKM~cate, lncrna_repeat_overlap_category_class_merged_copy_merged_with_no_rep)
TukeyHSD(anova_rep)



levels_for_plot<-c("five_prime","exon","intron",NA)
lncrna_repeat_overlap_category_class_merged_copy_merged_with_no_rep[,element:=factor(element,levels = levels_for_plot)]

fpkm_lncrna_repeats_barplot<-ggplot(data=lncrna_repeat_overlap_category_class_merged_copy_merged_with_no_rep,aes(x=repeat_class,fill=element)) +
  geom_bar(position="dodge",col="black",alpha=0.8) +
  theme_bw() +
  facet_wrap(.~experiment) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  theme(legend.position = "top") +
  scale_fill_got(discrete = T,option = "martell",na.value="grey",labels=c("5' regija","egzon","intron","Bez transpozona"), name="Element",direction = -1) +
  xlab("Razred transpozona") +
  ylab("Broj lncRNA-kodirajuæih gena") +
  theme(legend.title = element_text(size=9,face="bold")) +
  scale_x_discrete(labels=c("Unknown","DNA","LINE","LTR","Bez transpozona")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x  = element_blank())


fpkm_lncrna_repeats_barplot
ggsave(fpkm_lncrna_repeats_barplot,file="fpkm_lncrna_repeats_barplot.jpg",width=5,height = 4)

lncrna_repeat_overlap_category_class_merged_copy_merged_with_no_rep_cp<-copy(lncrna_repeat_overlap_category_class_merged_copy_merged_with_no_rep)
lncrna_repeat_overlap_category_class_merged_copy_merged_with_no_rep_cp[,cate:=paste(repeat_class,element,sep="_")]
lncrna_repeat_overlap_category_class_merged_copy_merged_with_no_rep_cp[,FPKM:=log10(FPKM)]
lncrna_repeat_overlap_category_class_merged_copy_merged_with_no_rep_cp[is.na(element),element:="no"]
lncrna_repeat_overlap_category_class_merged_copy_merged_with_no_rep_cp_day1<-lncrna_repeat_overlap_category_class_merged_copy_merged_with_no_rep_cp[experiment=="Prvi dan"]
anova_fpkm<-aov(FPKM~cate,data = lncrna_repeat_overlap_category_class_merged_copy_merged_with_no_rep_cp_day1)
TukeyHSD(anova_fpkm)

#overlapping proteisn vs non-overlapping proteins comparison


read_counts_table_longer_prot<-read_counts_table_longer[category=="protein_coding_gene"]
setnames(prot_classification_table,c("gene_id","lncrna_class"))

read_counts_table_longer_prot_merged<-merge(read_counts_table_longer_prot,prot_classification_table,by="gene_id")
read_counts_table_longer_prot_merged[,lncrna_class:=factor(lncrna_class,levels = c("intronic","intergenic_1kb","overlapping","no_lncRNA"))]
read_counts_table_longer_prot_merged[experiment=="day1",experiment:="Prvi dan"]
read_counts_table_longer_prot_merged[experiment=="day10",experiment:="Deseti dan"]
read_counts_table_longer_prot_merged[,experiment:=factor(experiment,levels = c("Prvi dan","Deseti dan"))]
options(scipen=1000)
fpkm_prot_plot<-ggplot(data=read_counts_table_longer_prot_merged,aes(x=lncrna_class,y=FPKM,fill=lncrna_class)) +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(.~experiment) +
  scale_y_log10() +
  theme(legend.title = element_text(size=9,face="bold")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x =element_blank()) +
  scale_fill_hp(discrete = TRUE, option = "NewtScamander", name = "Kategorija",labels=c("Proteini s\nintronskom lncRNA","Proteini s intergenskom\nlncRNA bližom od 1 kb","Proteini u\npreklapajuæoj lncRNA","Ostal\nproteini")) +
  theme(legend.text = element_text(size=9)) 

ggsave(fpkm_prot_plot,file="fpkm_prot_plot.jpg",width = 7,height = 4)

anova_class_prot<-aov(FPKM~lncrna_class,data=read_counts_table_longer_prot_merged)
TukeyHSD(anova_class_prot)