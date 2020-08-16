args <- commandArgs(trailingOnly = TRUE)  
library(data.table)
library(harrypotter)
setwd(".")
#randfold results analysis - take the lower p-value for each sequence-revcomp pair + make a distribution histogram of p-values + compare distributions of p-values across lncRNA classes
#Arguments: 1 - randfold results, 2 - randfold_revcomplement results, 3 - lncRNA classification table

randfold_res<-fread(args[1],header = F)
randfold_res_revcomp<-fread(args[2],header = F)
lncRNA_classification_table<-fread(args[3])

#taking the lover value out of sequence-revcomp pairs
setnames(randfold_res,c("transcript_name","energy","p_value"))
setnames(randfold_res_revcomp,c("transcript_name","energy","p_value"))
randfold_res[,category:="norm"]
randfold_res_revcomp[,category:="revcomp"]
randfold_together<-rbind(randfold_res,randfold_res_revcomp)
randfold_together_one<-randfold_together[,.(p_value=min(p_value)),by=transcript_name]


#p-value distribution

pval_dist_plot<-ggplot(randfold_together_one,aes(x=p_value)) +
  geom_histogram(fill="indianred2",col="black") +
  theme_bw() +
  xlab("p-vrijednost") +
  ylab("Broj molekula lncRNA")
pval_dist_plot
ggsave(pval_dist_plot,file="pval_dist_plot.jpg")


#p-value distribution comparison across classes of lncRNA genes
randfold_together_one_merged<-merge(randfold_together_one,lncRNA_classification_table,by="transcript_name")

class_pval_plot<-ggplot(randfold_together_one_merged,aes(x=class,y=p_value,fill=class)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_hp(discrete = TRUE, option = "NewtScamander", name = "Klasa",labels=c("Intergenske lncRNA","Intronske lncRNA","Preklapajuæe lncRNA")) +
  theme(legend.position = "top") +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank())  +
  ylab("p-vrijednost") +
  theme(legend.title = element_text(face = "bold",size = 10)) +
  theme(legend.text = element_text(size=9.5))
class_pval_plot

ggsave(class_pval_plot,file="class_pval_plot.jpg",width = 6,height = 5)


randfold_together_one_merged[,log_p_value:=log10(p_value)]
anova_p<-aov(log_p_value~class,data=randfold_together_one_merged)

TukeyHSD(anova_p)
