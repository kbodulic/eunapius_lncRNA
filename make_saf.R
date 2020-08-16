#SAF format for featurecounts

args <- commandArgs(trailingOnly = TRUE)  
library(data.table)
setwd(".")
#defining the annotation file for FeatureCounts (SAF)

#Arguments: 1 - protein annotation with filtered bacterial scaffolds (gtf), 2 - lncrna exonic coordinates, 3 - isoform-to-gene mapping table

protein_coding_annotation<-fread[args[1]]
protein_exons<-protein_coding_annotation[feature=="CDS"]

protein_exons_dt_saf<-protein_exons_dt[,.SD,.SDcols=c("attribute","seqnames","start","end","strand_ignore")]
protein_exons_dt_saf[,attribute:=str_extract(attribute,".*(?=\\.t\\d*)")]
setnames(protein_exons_dt_saf,c("GeneID","Chr","Start","End","Strand"))

filtered_together_paf_union_exons<-fread(args[2])
isoform_map_table_together<-fread(args[3])
filtered_together_paf_union_exons_saf<-filtered_together_paf_union_exons[,.SD,.SDcols=c("transcript","chromosome","start","end")]
setnames(filtered_together_paf_union_exons_saf,c("transcript_name","chromosome","start","end"))
filtered_together_paf_union_exons_saf<-merge(filtered_together_paf_union_exons_saf,isoform_map_table_together,by="transcript_name")
filtered_together_paf_union_exons_saf[,strand:="+"]
filtered_together_paf_union_exons_saf$transcript_name<-filtered_together_paf_union_exons_saf$gene_id
filtered_together_paf_union_exons_saf$gene_id<-NULL
setnames(filtered_together_paf_union_exons_saf,c("GeneID","Chr","Start","End","Strand"))


saf_together<-unique(rbind(protein_exons_dt_saf,filtered_together_paf_union_exons_saf))

write.table(saf_together,file="lncrna_and_prot_annotation.saf",row.names = F,quote = F,sep = "\t")