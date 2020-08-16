args <- commandArgs(trailingOnly = TRUE)  
library(data.table)
library(IRanges
library(GenomicRanges)
library(ggplot2
setwd(".")
#mapping filtering - removing transcrips which map on bacterial scaffolds + transcripts with one exon + transcripts with less than 95% identity to the genome + transcripts whose exons overlap protein-coding exons
#Arguments: 1 - mapping results (paf), 2 - taxonomy table (MEGAN), 3 - exonic coordinates (retrived from samtools mpileup and a custom R script), 4 - protein coding annotation (gtf), 5 - transcriptome file name, 6 - first-half-filtering numbers

mapping_paf<-fread(args[1])
taxonomy_table<-fread(args[2])
lncrna_exons<-fread(args[3])
protein_coding_annotation<-fread(args[4],fill=T)[is.na(V4)==F]
#arg 5 only for the name of a table which will be writen
first_half_filtering_info<-fread(args[6],header = F)



setnames(mapping_paf,c("query_name","query_length","query_start","query_end","strand_ignore","target_name","target_length","target_start","target_end","matches","alignment_block","mapping_quality"))

mapping_paf[,numb_same_trans:=nrow(.SD),by=query_name]
mapping_paf[numb_same_trans>1,query_name:=paste(query_name,numb_same_trans,sep="_")]
mapping_paf$numb_same_trans<-NULL


#bacterial scaffolds filtering 



taxonomy_table_nonbact<-unique(taxonomy_table[my_taxonomy=="Metazoa" | my_taxonomy=="Other"]$seqid)

mappings_paf_nobact<-mapping_paf[target_name %in% taxonomy_table_nonbact]


setnames(protein_coding_annotation,c("seqname","source","feature","start","end","score","strand_ignore","frame","attribute"))

protein_coding_annotation<-protein_coding_annotation[seqname%in%taxonomy_table_nonbact]
#writing only non-bacterial scaffolds
write.table(protein_coding_annotation,file="protein_coding_annotation_nonbacteria.txt",row.names=F,quote=F)

#one exon filtering


lncrna_more_than_1_exon<-lncrna_exons[,.N,by=transcript][N!=1]$transcript


mapping_lncrna_paf_more_than_one_exon<-mappings_paf_nobact[query_name%in%lncrna_more_than_1_exon]


#perc matches filtering
mapping_lncrna_paf_more_than_one_exon[,perc_matches:=100*(matches/query_length)]



mapping_lncrna_paf_more_than_one_exon_percmatches_95<-mapping_lncrna_paf_more_than_one_exon[perc_matches>=95]

filtered_exons<-lncrna_exons[transcript%in%lncrna_more_than_1_exon & transcript%in%mapping_lncrna_paf_more_than_one_exon_percmatches_95$query_name]
filtered_exons_granges<-makeGRangesFromDataFrame(filtered_exons,seqnames.field = "chromosome",ignore.strand = T,keep.extra.columns = T)

#protein exons overlap filtering
protein_exons_dt<-protein_coding_annotation[feature=="CDS"]

protein_exons<-makeGRangesFromDataFrame(df=protein_exons_dt,keep.extra.columns = T,ignore.strand = T)





overlaps_lncrnaexon_genomexone<-findOverlaps(filtered_exons_granges,protein_exons)

transcripts_coding_overlap<-unique(filtered_exons_granges[queryHits(overlaps_lncrnaexon_genomexone)]$transcript)



mapping_paf_filtered<-mapping_lncrna_paf_more_than_one_exon_percmatches_95[query_name%in%transcripts_coding_overlap==F]

filtered_transcript_list<-mapping_paf_filtered$query_name

filtered_exons<-filtered_exons[transcript%in%filtered_transcript_list]
write.table(filtered_exons,file=paste(args[3],"filtered",sep="_"),row.names = F,quote = F)

write.table(filtered_transcript_list,file = paste(args[5],"filtered_list",sep="_"),row.names = F,col.names = F,quote = F,sep = "\t")


write.table(mapping_paf_filtered,file = paste(args[5],"filtered_paf",sep="_"),row.names = F,quote = F)


lengths_table<-data.table(filtering=c("Svi transkripti","Uklanjanje rRNA","Filtriranje po duljini","Filtriranje po duljini okvira èitanja","Filtriranje DIAMOND pogodaka","Filtriranje HMMER pogodaka","Filtriranje programom FEELnc","Filtriranje bakterijskih sljedova ","Filtriranje po broju egzona","Filtriranje po mapiraju","Filtriranje po preklapanju s PK egzonima"),number_of_lncRNA_candidates=c(first_half_filtering_info$V1,nrow(mappings_paf_nobact),length(lncrna_more_than_1_exon),nrow(mapping_lncrna_paf_more_than_one_exon_percmatches_95),length(filtered_transcript_list))
)

lengths_plot<-ggplot(data=lengths_table,aes(x=filtering,y=number_of_lncRNA_candidates,fill=filtering)) +
  geom_bar(position="dodge",stat = "identity", width=0.65,col="black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 50, hjust = 1,size=7.5)) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  ylab("Broj transkripata") +
  scale_x_discrete(limits=c("Svi transkripti","Uklanjanje rRNA","Filtriranje po duljini","Filtriranje po duljini okvira èitanja","Filtriranje DIAMOND pogodaka","Filtriranje HMMER pogodaka","Filtriranje programom FEELnc","Filtriranje bakterijskih sljedova ","Filtriranje po broju egzona","Filtriranje po mapiraju","Filtriranje po preklapanju s PK egzonima")) 
  
  
  
  ggsave(lengths_plot,file="filtering_lengths_plot.jpg")