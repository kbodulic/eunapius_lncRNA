args <- commandArgs(trailingOnly = TRUE)  
library(data.table)
library(stringr)
setwd(".")
get_exons<-function(pile_path) {
  #Loading the resulting table of samtools mpileup -A --output-QNAME

  pile<-fread(pile_path)
  #Removing unneccesary columns to save on RAM
  pile$V3<-NULL
  pile$V4<-NULL
  pile$V6<-NULL
  #The non-removed part of the pileup format contains 4 columns: chromosome (contig) name, base position,
  #status of the base and read name. The base status column contains information on the mapping of every
  #read base to the particular position on the reference. Matches and missmatches are represended by A, T, C,
  #G or N and a, t ,g ,c or n, with capital letters denoting the forward reference strand and lower case
  #letters denoting the reverse reference strand). Futhermore, small reference insertions are represented
  # with a + and a number denoting the insertion lengths, while small deletions are represented with a -
  #followed by a number denoting the deletion length. Introns (longer insertions) are represented with
  # < and >. If there are multiple reads mapping to the same base, the base statuses of those reads are
  #written together (AA for example) and the reads are separated with a comma. Generally, introns are
  #regarded as insertions larger or equal to 40 bp. However, there are a few examples with larger deletions
  #than this which are regarded as deletions (not as introns) - all of these examples correspond to poorly
  #mapped transcripts and shouldn't represent a problem.
  setnames(pile,c("chrom_name","base","base_status","transcript"))
  #For some reason, a small number of rows has NA in the base_status column - the next line will remove
  #those positions. 
  pile<-pile[is.na(base_status)==F]
  #Next five lines remove the unneccesary special characters which rarely come up in the base_status
  #columns. These characters, for example, represent the begining of a read(^), the ending of a read
  #($) and so on. It's important to remove these characters because they will cause a missalignment in
  #the base_status and read_name columns (one special character doesn't necessarily correspond to one
  #read name - important fot the next part of the script.). Also, I can't simply str_remove all of those
  #characters because some of them (+ and - for example) are usually followed by more special characters 
  #- it's better to simply remove those columns and account for that latter. I str_removed ^ and $ because 
  #$ is commonly followed by ^ (end of a read is commonly followed by a begining of a new read (which 
  #causes the first and  last exons of reads to have a shift of 1.
   beg_end<-c("\\^","\\$")
  pile[,base_status:=str_remove_all(base_status,beg_end)]
  chars<-c("\\^","\\$","#","]","!","\\\\","\\?","\\.","\\(","\\[","\\)",";","=","&","@",":","'")
  remove<-c(chars,letters[letters%in%c("a","t","c","g","n")==F],LETTERS[LETTERS%in%c("A","T","G","C","N")==F],0:9)
  pile<-pile[grepl(paste(remove,collapse="|"),base_status)==F]
  
  #Counting the maximum number of characters in the base_status column
  max_letters<-max(nchar(pile$base_status))
  #Counting the maximum number of reads
  max_transcripts<-max(str_count(pile$transcript,","))+1
  #Splitting the base status column into multiple columns containing the base of every read 
  #(each base is in its own column)
  pile[,paste("base_status",1:max_letters,sep="_"):=tstrsplit(base_status,"")]
  #Spliting the read column into multiple columns so that each read is in its own column
  pile[,paste("transcript",1:max_transcripts,sep="_"):=tstrsplit(transcript,",")]
  #Reducing this table into a table with 4 columns (chromosome, base position, base status and read name)
  #with each row corresponding to one transcript
  pile_list<-lapply(1:max_transcripts, function(x) {
    dt_sep<-pile[,.SD,.SDcols=c(1,2,4+x,4+max_letters+x)]
    dt_sep<-dt_sep[as.vector(is.na(dt_sep[,3]))==F & as.vector(is.na(dt_sep[,4]))==F]
    return(dt_sep)
  })
  rm(pile)
  pile_list_table<-rbindlist(pile_list)
  rm(pile_list)
  
  #Sorting the table by chromosome, transcript and base position
  pile_list_table<-pile_list_table[order(chrom_name,transcript_1,base)]
  #For some reason, mpileup occasionally outputs the same transcript mapping on the same base twice 
  #(secondary alignments?), these need to be removed - explanation latter. 
   pile_list_table<-pile_list_table[,red:=1:nrow(pile_list_table)]
  rm_rows<-pile_list_table[transcript_1==data.table::shift(transcript_1,n=1,type="lead") & base==data.table::shift(base,n=1,type="lead")]
  pile_list_table<-fsetdiff(pile_list_table,rm_rows)
  pile_list_table$red<-NULL
  #Removing the intronic bases
  pile_list_table<-pile_list_table[base_status_1!=">" & base_status_1!="<"]
  #Calculating the difference of base positions after removing introns
  pile_list_table[,diffe:=c(0,diff(base))]
  #Now, the difference of adjacent exonic bases should be one. However, since I removed the columns containing
  #special symbols and NA's (in the begining of the script), it's possible that a base will be skipped,
  #which is why I allowed this difference to be up to 15 (in case the special symbols occured in more
  #consecutive columns (according to literature, the minimal intron length is 40-50 bp, so this is fine.)
  pile_list_table[,cond:=ifelse(diffe>=1 & diffe <30,T,F)]
  #Grouping the adjacent exonic bases into exons (represented by the number given by rleid)
  #This is why it was important to remove one of the rows corresponding to the same transcript mapping
  #on the same base multiple times -after sorting, diff of those rows would be 0 and they would come off as
  #exons of length 1
  pile_list_table[,condnumb:=rleid(cond)]
  #Odd rleid numbers will represent the first exon base (exonic boundaries) - the next line removes them
  pile_list_table<-pile_list_table[condnumb%%2==0]
  #Making a table with ranges of the exons - accaunting for the removal of the first exon base by
  #substracting 1 from the minima
  exon_coordinates<-pile_list_table[,.(start=min(base)-1,end=max(base)),by=c("transcript_1","condnumb","chrom_name")][order(-chrom_name,start)]
  exon_coordinates$condnumb<-NULL
  setnames(exon_coordinates,c("transcript","chromosome","start","end"))
  exon_coordinates[,start:=ifelse(start==.SD$start[1] ,start-2,start),by=transcript]
  return(exon_coordinates)
}

write.table(get_exons(args[1]),file=paste(args[1],"exons",sep="_"),quote = F,row.names = F)

#To conclude, this script workts quite well (I tested many transcripts and all of their exons were
#correctly found). There are some extremely small exons, but those correspond to actual mappings, not
#to a bug in the script (I checked manually). However, the script is really time and memory consuming,
#(and quite ugly) which is why I think that the cigar string approach would be much better and easier.
      
      
      
      
      
      
      
      
      
      
      
      
      
