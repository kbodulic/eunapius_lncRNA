

args <- commandArgs(trailingOnly = TRUE)  
library(data.table)
library(IRanges)
library(GenomicRanges)
library(ggplot2)
library(harrypotter)
setwd(".")
#Anylsis of blast hits in Porifera - produces a blast hit heatmap and a concatanated fasta files of transcript hit coordinates (part trans - only the parts which are similar, whole trans - whole transcripts)
#Arguments: 1 - blast reults

blast_results_sponge_transcriptomes<-fread(args[1])
setnames(blast_results_sponge_transcriptomes,c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","qlen","sstart","send","slen","evalue","bitscore"))

blast_results_sponge_transcriptomes[,qperclen:=100*(abs(qend-qstart)/qlen)]
blast_results_sponge_transcriptomes[grepl("genome",sseqid),species:=str_extract(sseqid,".*_genome")]
blast_results_sponge_transcriptomes[grepl("transcriptome",sseqid),species:=str_extract(sseqid,".*_transcriptome")]

#E-value filtering
blast_results_sponge_transcriptomes_eval_filt<-blast_results_sponge_transcriptomes[evalue<=1e-15]

blast_results_sponge_transcriptomes_eval_filt[,number_of_species_per_gene:=uniqueN(.SD),.SDcols="species",by=qseqid]
trans_more_than_one_species<-unique(blast_results_sponge_transcriptomes_eval_filt[number_of_species_per_gene>1]$qseqid)

#writing a table of transcrips which significantly hit more than one species
write.table(trans_more_than_one_species,file="trans_more_than_one_species.txt",row.names = F,col.names = F,quote = F,sep = "\t")



#reversing the rev_complement hits and reducing all hits on the same query
blast_results_sponge_transcriptomes_eval_filt[sstart>send,':='(sstart=send,send=sstart)]
blast_gr<-makeGRangesFromDataFrame(blast_results_sponge_transcriptomes_eval_filt,start.field = "qstart",end.field="qend",seqnames.field = "qseqid",ignore.strand=T,keep.extra.columns = T)
blast_gr_q_split<-split(blast_gr,blast_gr$sseqid)
blast_gr_q_split_reduce<-reduce(blast_gr_q_split)
blast_gr_q_split_reduce_unlist<-unlist(blast_gr_q_split_reduce,use.names = T)

blast_gr_q_split_reduce_unlist_dt<-as.data.table(blast_gr_q_split_reduce_unlist)
blast_gr_q_split_reduce_unlist_dt[,sseqid:=names(blast_gr_q_split_reduce_unlist)]
setnames(blast_gr_q_split_reduce_unlist_dt,c("qseqid","qstart","qend","width","strand","sseqid"))
blast_gr_q_split_reduce_unlist_dt_lengths<-merge(blast_gr_q_split_reduce_unlist_dt,unique(blast_results_sponge_transcriptomes_eval_filt[,.SD,.SDcols=c("qseqid","qlen")]),by="qseqid")
blast_gr_q_split_reduce_unlist_dt_lengths[,qperclen:=100*(width/qlen)]
blast_gr_q_split_reduce_unlist_dt_lengths[,sum_width:=sum(.SD),.SDcols="width",by=c("qseqid","sseqid")]
blast_gr_q_split_reduce_unlist_dt_lengths[,sumqperclen:=sum(.SD),.SDcols="qperclen",by=c("qseqid","sseqid")]
blast_gr_q_split_reduce_unlist_dt_lengths[grepl("genome",sseqid),species:=str_extract(sseqid,".*_genome")]
blast_gr_q_split_reduce_unlist_dt_lengths[grepl("transcriptome",sseqid),species:=str_extract(sseqid,".*_transcriptome")]

blast_gr_q_split_reduce_unlist_dt_lengths_one_per_species<-blast_gr_q_split_reduce_unlist_dt_lengths[,.(sumqperclen=max(sumqperclen)),by=c("qseqid","species")]



species_list<-c("Aphrocallistes_vastus_transcriptome","Chondrilla_nucula_transcriptome","Clathrina_sp_transcriptome","Corticium_candelabrum_transcriptome","Crella_elegans_transcriptome","Ephydatia_mulleri_transcriptome","Haliclona_amboinensis_transcriptome","Haliclona_indistincta","Haliclona_oculata_transcriptome","Haliclona_simulans_transcriptome","Haliclona_tubifera_transcriptome","Halisraca_caerulea_transcriptome","Ircinia_fasciculata_transcriptome","Leucosolenia_complicata_transcriptome","Oscarella_carmela_genome","Oscarella_sp_transcriptome","Petrosia_ficiformis_transcriptome","Pseudospongosorites_suberitoides_transcriptome","Stylissa_carteri_genome","Sycon_cilliatum_genome","Xestospongia_testudinaria_genome","Amphimedon_queenslandica_genome")
nonexistant_species<-setdiff(species_list,unique(blast_gr_q_split_reduce_unlist_dt_lengths_one_per_species$species))




blast_gr_q_split_reduce_unlist_dt_lengths_one_per_species<-rbind(blast_gr_q_split_reduce_unlist_dt_lengths_one_per_species,data.table(qseqid=rep("NODE_6217_length_4525_cov_17.657251_g1739_i1",length(nonexistant_species)),species=nonexistant_species,sumqperclen=rep(NA,length(nonexistant_species))))

species_ordered<-blast_gr_q_split_reduce_unlist_dt_lengths_one_per_species[,.N,by=species][order(-N)]$species
label_cols<-ifelse(grepl("genome",species_ordered)==T,"brown2","gray52")

blast_gr_q_split_reduce_unlist_dt_lengths_one_per_species[,species:=str_remove(species,"_genome")]
blast_gr_q_split_reduce_unlist_dt_lengths_one_per_species[,species:=str_remove(species,"_transcriptome")]
species_ordered2<-str_remove(species_ordered,"_genome")
species_ordered2<-unique(str_remove(species_ordered2,"_transcriptome"))
blast_heatmap_plot<-ggplot(data=blast_gr_q_split_reduce_unlist_dt_lengths_one_per_species,aes(x=species,y=qseqid,fill=sumqperclen)) +
  geom_tile() +
  theme_bw() +
  scale_fill_hp(house = "HermioneGranger",name="Omjer duljine BLAST pogotka i\nukupne duljine ulazne lncRNA (%) ",na.value="gray94") +
  ylab("Molekula lncRNA") + 
  theme(legend.position="top") +
  theme(legend.key.size = unit(0.3, "cm"),
        legend.key.width = unit(0.5,"cm")) +
  theme(legend.title = element_text( size = 9, face="bold")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size=9.2)) +
  theme(axis.ticks.y=element_blank()) +
  theme(axis.text.y=element_blank())  +
  scale_x_discrete(limits=species_ordered2) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size=7,colour = label_cols)) 


ggsave(blast_heatmap_plot,file="blast_heatmap_plot.jpg",width = 7,height = 6)


#findig candidatess for multiple alignment

blast_gr_q_split_reduce_unlist_dt_lengths<-blast_gr_q_split_reduce_unlist_dt_lengths[grepl("transcriptome",species)]

blast_gr_q_split_reduce_unlist_dt_lengths[,is_max_sumqperclen:=ifelse(.SD==max(.SD),T,F),.SDcols="sumqperclen",by=c("qseqid","species")]

only_max_sumperceln_width_filtered<-blast_gr_q_split_reduce_unlist_dt_lengths[is_max_sumqperclen==T & sum_width>=150]


#reducing all hits on the subject queries - finding the coordinates of the parts which are similar
sumandwidth_filtered<-merge(blast_results_sponge_transcriptomes_eval_filt,unique(only_max_sumperceln_width_filtered[,.SD,.SDcols=c("qseqid","sseqid","is_max_sumqperclen")]),by=c("qseqid","sseqid"),all.x=T)[is_max_sumqperclen==T]
sumandwidth_filtered<-sumandwidth_filtered[order(evalue)]
sumandwidth_filtered<-unique(sumandwidth_filtered,by=c("qseqid","species"))
sumandwidth_filtered[sstart>30,sstart:=sstart-30]
sumandwidth_filtered[send<slen-30,send:=send+30]
sumandwidth_filtered_gr<-makeGRangesFromDataFrame(sumandwidth_filtered,seqnames.field="sseqid",start.field = "sstart",end.field = "send",ignore.strand = T,keep.extra.columns = T)

sumandwidth_filtered_gr_split<-split(sumandwidth_filtered_gr,sumandwidth_filtered_gr$qseqid)

sumandwidth_filtered_gr_split_range<-range(sumandwidth_filtered_gr_split)                                                                                 

sumandwidth_filtered_gr_split_range_unlist<-unlist(sumandwidth_filtered_gr_split_range)

sumandwidth_filtered_gr_split_range_unlist_table<-as.data.table(sumandwidth_filtered_gr_split_range_unlist)[,.SD,.SDcols=c("seqnames","start","end")]
sumandwidth_filtered_gr_split_range_unlist_table[,transcript:=names(sumandwidth_filtered_gr_split_range_unlist)]
#blast - bed conversion 

sumandwidth_filtered_gr_split_range_unlist_table[,start:=start-1]

sumandwidth_filtered_gr_split_range_unlist_table[,species_num:=nrow(.SD),by=transcript]
sumandwidth_filtered_gr_split_range_unlist_table<-sumandwidth_filtered_gr_split_range_unlist_table[species_num>1]
sumandwidth_filtered_gr_split_range_unlist_table$species_num<-NULL


sumandwidth_filtered_gr_split_range_unlist_table_split<-split(sumandwidth_filtered_gr_split_range_unlist_table,sumandwidth_filtered_gr_split_range_unlist_table$transcript)


for (i in sumandwidth_filtered_gr_split_range_unlist_table_split) {
  write.table(i,file=paste(i$transcript[1],"msa_part_trans.txt",sep="_"),row.names = F, col.names = F, quote = F,sep = "\t")
}

#finding the coordinates of whole transcripts 
sumandwidth_filtered_gr_split_range_unlist_table_for_merge<-copy(sumandwidth_filtered_gr_split_range_unlist_table)

setnames(sumandwidth_filtered_gr_split_range_unlist_table_for_merge,c("sseqid","sstart","send","qseqid"))
sumandwidth_filtered_gr_split_range_unlist_table_merged<-merge(sumandwidth_filtered_gr_split_range_unlist_table_for_merge,unique(blast_results_sponge_transcriptomes[,.SD,.SDcols=c("sseqid","qseqid","slen")],by=c("sseqid","qseqid")))

sumandwidth_filtered_gr_split_range_unlist_table_merged[,':='(sstart=0,send=slen)]
sumandwidth_filtered_gr_split_range_unlist_table_merged$slen<-NULL
sumandwidth_filtered_gr_split_range_unlist_table_merged<-sumandwidth_filtered_gr_split_range_unlist_table_merged[,.SD,.SDcols=c("sseqid","sstart","send","qseqid")]
sumandwidth_filtered_gr_split_range_unlist_table_whole_trans_split<-split(sumandwidth_filtered_gr_split_range_unlist_table_merged,sumandwidth_filtered_gr_split_range_unlist_table_merged$qseqid)

for (i in sumandwidth_filtered_gr_split_range_unlist_table_whole_trans_split) {
  write.table(i,file=paste(i$qseqid[1],"msa_whole_trans.txt",sep="_"),row.names = F, col.names = F, quote = F,sep = "\t")
}

