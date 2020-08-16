#!/bin/bash

#Arguments: 1 - transcriptome, 2 - mRNA list (fasta, 3 - genome assembly (fasta), 4 - taxonomy table (MEGAN), 5 - protein annotation (gtf)


#Filtering rRNA (BLASTn) + by length (seqkit)

MAKEBLASTDB=/path/to/makeblastdb
INPUTDB=/path/to/Silva/rRNA/database

$MAKEBLASTDB -in $INPUTDB -dbtype nucl -title sponge_rrna_both_subunits

BLASTN=/path/to/blastn

QUERY=$1
KNOWN_MRNA=$2
GENOME=$3
TAXONOMY_TABLE=$4
PROTEIN_CODING_ANNOTATION=$5

RRNA_BLAST_OUTPUT="${QUERY##*.}"_blast_rRNA.outfmt6

$BLASTN -task blastn -db $INPUTDB -query $QUERY -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore" -num_threads 10 -out $RRNA_BLAST_OUTPUT

#Finding the length of transcripts 
SEQKIT=/path/to/seqkit

$SEQKIT fx2tab -n -l $QUERY > "${QUERY##*.}"_table

#R script for filtering rRNA based on blast results and filtering by length based on seqkit data
#Arguments: 1 - seqkit length table, 2 - BLASTn results (rRNA database)
Rscript filter_rrna_and_length.R "${QUERY##*.}"_table $RRNA_BLAST_OUTPUT

$SEQKIT grep --pattern-file "${QUERY##*.}"_table_without_rrna_more_than_200 $QUERY > "${QUERY##*.}"_without_rrna_more_than_200.fa

#ORF length filtering


LONGORFS=/path/to/longorfs

PREDICT=/path/to/predict

$PREDICT -t "${QUERY##*.}"_without_rrna_more_than_200.fa

CDSFILE="${QUERY##*.}"_without_rrna_more_than_200_transdecoder.cds

#making a table of ORF lengths - stolen from Dunja Glavaš :) 

  grep ">" $CDSFILE | cut -d " " -f 1,3-7 | awk '{ gsub(">", "", $1); print }' | awk '{ gsub("type:", "", $3); print }' | awk '{ gsub("len:", "", $4); print }' | sed 's/(//' | sed 's/),/ /' | sed 's/ /\t/g' > ${CDSFILE}_table.csv
  
#R script for ORF length filtering -  remove transcripts which have no ORFS and transcripts whose longest ORF is shorter than 150 nucleotides (50 amino acids
#Arguments: 1 - rRNA and length filtered transcript table, 2 - Transdecoder.cds table 
Rscript filter_orfs.R "${QUERY##*.}"_table_without_rrna_more_than_200 ${CDSFILE}_table.csv 

$SEQKIT grep --pattern-file "${QUERY##*.}"_table_rrna_length_orfs_filtered "${QUERY##*.}"_without_rrna_more_than_200.fa >  "${QUERY##*.}"_table_rrna_length_orfs_filtered.fa

#Filtering by similarity to known proteins using DIAMOND

DIAMOND=path/to/diamond
DATABASE=/path/to/nr/database

DIAMOND_INPUT="${QUERY##*.}"_table_rrna_length_orfs_filtered.fa

DIAMOND_OUTPUT="${QUERY##*.}"_table_rrna_length_orfs_filtered_diamond_res.outfmt6

$DIAMOND blastx -d $DATABASE -q $DIAMOND_INPUT -o $DIAMOND_OUTPUT -f 6 -p 10 

#R script for filtering by DIAMOND results - remove transcripts which show a significant similarity to known proteeins (NR database)
#Arguments: 1 - List of transcripts fitlered by rRNA, length and ORF length, 2 - DIAMOND results
Rscript filter_diamond.R "${QUERY##*.}"_table_rrna_length_orfs_filtered $DIAMOND_OUTPUT 

$SEQKIT grep --pattern-file "${QUERY##*.}"_table_rrna_length_orfs_filtered_without_diamond "${QUERY##*.}"_table_rrna_length_orfs_filtered.fa > "${QUERY##*.}"_table_rrna_length_orfs_filtered_without_diamond.fa

#hmmer filtering - removing transcripts with a significant similarity to conserved protein domains in Pfam-A

HMMSCAN=/path/to/hmmscan

HMMSCAN_INPUT="${QUERY##*.}"_without_rrna_more_than_200_transdecoder.pep

PFAMDATABASE=/path/to/pfam/database

HMMSCAN_OUTPUT="${QUERY##*.}"_table_rrna_length_orfs_filtered_hmmer_res


$HMMSCAN --cpu 12 --domtblout $HMMSCAN_OUTPUT $PFAMDATABASE $HMMSCAN_INPUT > ${HMMSCAN_OUTPUT}_log.txt

#R script for filtering hmmscan results - remove translated transcripts with a significant similarity to conserved protein domains
#Arguments: 1 - table of transcripts filtered by rRNA, length, ORF length and DIAMOND hits, 2 - hmmscan results 
Rscript filter_hmmer.R "${QUERY##*.}"_table_rrna_length_orfs_filtered_without_diamond $HMMSCAN_OUTPUT 

$SEQKIT grep --pattern-file "${QUERY##*.}"_table_rrna_length_orfs_filtered_without_diamond_no_hmmer "${QUERY##*.}"_table_rrna_length_orfs_filtered_without_diamond.fa > "${QUERY##*.}"_table_rrna_length_orfs_filtered_without_diamond_no_hmmer

#FEELnc - finding lncRNAS based on nucleotide composition

FEELNC_CODPOT=/path/to/feelnc/codpot

$FEELNC_CODPOT -i "${QUERY##*.}"_table_rrna_length_orfs_filtered_without_diamond_no_hmmer -a $KNOWN_MRNA --mode=shuffle

mv feelnc_codpot_out/"${QUERY##*.}"_table_rrna_length_orfs_filtered_without_diamond_no_hmmer.lncrna.fa .

#mapping

MINIMAP2=path/to/minimap2


$MINIMAP2 -x splice -K 10000M -c --secondary=no -t 12 $GENOME "${QUERY##*.}"_table_rrna_length_orfs_filtered_without_diamond_no_hmmer.lncrna.fa -o "${QUERY##*.}"_filtered_before_mapping_to_genome.paf

$MINIMAP2 -ax splice -K 10000M --secondary=no  -t 12 $GENOME "${QUERY##*.}"_table_rrna_length_orfs_filtered_without_diamond_no_hmmer.lncrna.fa -o "${QUERY##*.}"_filtered_before_mapping_to_genome.sam

#converting and sorting

SAMBAMBA=path/to/sambamba
$SAMBAMBA view -t 12 -f bam -S -o "${QUERY##*.}"_filtered_before_mapping_to_genome.bam "${QUERY##*.}"_filtered_before_mapping_to_genome.sam
$SAMBAMBA sort -m 80GB --tmpdir tmp -t 12 "${QUERY##*.}"_filtered_before_mapping_to_genome.bam
$SAMBAMBA flagstat -t 12 "${QUERY##*.}"_filtered_before_mapping_to_genome.sorted.bam > "${QUERY##*.}"_filtered_before_mapping_to_genome_flagstat.txt
rm "${QUERY##*.}"_filtered_before_mapping_to_genome.sam
#get exonic coordinates from mapping file 

SAMTOOLS=/path/to/samtools
$SAMTOOLS mpileup -A --output-QNAME "${QUERY##*.}"_filtered_before_mapping_to_genome.sorted.bam > "${QUERY##*.}"_filtered_before_mapping_to_genome_mpileup

#R script for extracting exonic coordinates from a mpileup file
#Arguments: 1 - Mapping mpileup format (samtools)
Rscript get_exons.R "${QUERY##*.}"_filtered_before_mapping_to_genome_mpileup

#getting the numbers of filtered transcripts
cat "${QUERY##*.}"_table_info_table_rrna_length $(wc -l "${QUERY##*.}"_table_rrna_length_orfs_filtered) $( wc -l "${QUERY##*.}"_table_rrna_length_orfs_filtered_without_diamond) $(wc -l "${QUERY##*.}"_table_rrna_length_orfs_filtered_without_diamond_no_hmmer) $(grep -c ">" "${QUERY##*.}"_table_rrna_length_orfs_filtered_without_diamond_no_hmmer.lncrna.fa) > "${QUERY##*.}"_first_half_filtering_info.txt

#filtering bacterial transcripts, one-exon transcript, transcripts which mapp on the genome with less than 95% identity and transcripts whose exons overtlap with prot exons
#also, producing a plot which shows the filtered number of transcripts by each category in this script
#Arguments: 1 - mapping results (paf), 2 - taxonomy table (MEGAN), 3 - exonic coordinates (retrived from samtools mpileup and a custom R script), 4 - protein coding annotation (gtf), 5 - transcriptome file name, 6 - first-half-filtering numbers
Rscript mapping_filtering.R "${QUERY##*.}"_filtered_before_mapping_to_genome.paf $TAXONOMY_TABLE "${QUERY##*.}"_filtered_before_mapping_to_genome_mpileup_exons $PROTEIN_CODING_ANNOTATION "${QUERY##*.}" "${QUERY##*.}"_first_half_filtering_info.txt










