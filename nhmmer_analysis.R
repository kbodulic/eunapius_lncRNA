args <- commandArgs(trailingOnly = TRUE)  
library(data.table)
library(ggplot2)
library(taxize)
setwd(".")
#conservation analysis (outside phylum Porifera) - nhmmer results (rnacentral database without rRNA)
#Arguments: 1 - nhmmer results
rnacentral_nhmmer_results<-fread(args[1],fill = T,header = F)
#nhmmer produces a weird output - tiding up the table
setnames(rnacentral_nhmmer_results,c("target_name","accession","query_name","accession2","qstart","qend","sstart","send","envstart","endnd","slen","starnd","evalue","score","bias","genus","species","RNAtype","RNAtype2","RNAtype3","V24","V25","V26"))
#filtering for significant e-values
rnacentral_nhmmer_results_eval_filt<-rnacentral_nhmmer_results[evalue<=1e-15]
rnacentral_nhmmer_results_eval_filt[grepl("sp",RNAtype),':='(genus=species,species=RNAtype)]
rnacentral_nhmmer_results_eval_filt[grepl("RNA",RNAtype2),RNAtype:=RNAtype2]
rnacentral_nhmmer_results_eval_filt[grepl("RNA",RNAtype3),RNAtype:=RNAtype3]
table(rnacentral_nhmmer_results_eval_filt$RNAtype)
rnacentral_nhmmer_results_eval_filt[,binary_nom:=paste(genus,species,sep=" ")]
#getting higher taxons
uids <- get_uid(unique(rnacentral_nhmmer_results_eval_filt$binary_nom[1:nrow(rnacentral_nhmmer_results_eval_filt)]))
classifications<-classification(uids)



classifications<-classifications[is.na(classifications)==F]
tax_table<-data.table()
for (i in 1:length(classifications)) {
  categories_table<-data.table(species=as.data.table(classifications[[i]])[rank=="species"]$name,phylum=as.data.table(classifications[[i]])[rank=="phylum"]$name,superkingdom=as.data.table(classifications[[i]])[rank=="superkingdom"]$name)
  tax_table<-rbind(tax_table,categories_table)
}

setnames(tax_table,c("binary_nom","phylum","superkingdom"))

#plotting the number of significant hits
rnacentral_nhmmer_results_eval_filt_merged<-merge(rnacentral_nhmmer_results_eval_filt,tax_table,by="binary_nom")
rnacentral_nhmmer_results_eval_filt_merged<-rnacentral_nhmmer_results_eval_filt_merged[is.na(phylum)==F]
rnacentral_nhmmer_results_eval_filt_merged<-rnacentral_nhmmer_results_eval_filt_merged[RNAtype=="lncRNA" | RNAtype=="misc_RNA"]
rnacentral_nhmmer_results_eval_filt_merged[,is_min_eval:=ifelse(.SD==min(.SD),T,F),.SDcols="evalue",by=c("binary_nom","query_name")]
rnacentral_nhmmer_results_eval_filt_merged<-rnacentral_nhmmer_results_eval_filt_merged[is_min_eval==T]

rnacentral_num_plot<-ggplot(data=rnacentral_nhmmer_results_eval_filt_merged,aes(x=phylum,fill=RNAtype)) +
  geom_bar(col="black",width=0.4) +
  theme_bw() +
  scale_fill_manual(values = c("indianred2","orange")) 

ggsave(rnacentral_num_plot,file="rnacentral_num_plot.jpg")

rnacentral_nhmmer_results_eval_filt_merged[,log_eval:=-log10(evalue)]

#plotting the distribution of e-values of significant hits 
rnacentral_eval_plot<-ggplot(data=rnacentral_nhmmer_results_eval_filt_merged,aes(x=phylum,y=log_eval,fill=RNAtype)) +
  geom_boxplot(width=0.4) +
  theme_bw() +
  scale_fill_manual(values = c("indianred2","orange")) 
ggsave(rnacentral_eval_plot,file="rnacentral_eval_plot.jpg")

