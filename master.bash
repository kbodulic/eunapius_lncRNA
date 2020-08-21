#!/bin/bash
#master script
#Arguments: 1 - left-paired reads of RNA1, 2 - right-paired reads of RNA1, 3 - left-paired reads of RNA10, 4 - rught-paired reads of RNA10, 5 - known mRNA (fasta),6 -genome assebmly (fasta), 7 - taxonomy table (MEGAN), 8 - protein annotation (gtf), 9 - repeats annotation (repeatmasker), 10 - database of sponge datasets (fasta)

IN11=$1
IN12=$2
IN3=$3
IN4=$4
KNOWN_MRNA=$5
GENOME=$6
TAXONOMY_TABLE=$7
PROTEIN_CODING_ANNOTATION=$8
REPEAT_ANNOTATION=$9
SPONGE_DATABASE=$10

# remove adapters, perform light quality trim (remove bases from ends with quality score under 8) and filter out reads with average quality after trimming below 16:
BBTOOLS=path/to/BBduk
# reference file for adapters
ADAPTERS=/path/to/adapters
LIBNAMEA1=RNA1
LIBNAME=RNA10
IN1_TRIMMED=${LIBNAME1}01_trimfilt.fq.gz
IN2_TRIMMED=${LIBNAME1}02_trimfilt.fq.gz
IN3_TRIMMED=${LIBNAME2}01_trimfilt.fq.gz
IN4_TRIMMED=${LIBNAME2}02_trimfilt.fq.gz
$BBDUK-Xmx60g in1=$IN1 in2=$IN2 out1=$IN1_TRIMMED out2=$IN2_TRIMMED ref=$ADAPTERS threads=4 ktrim=r mink=10 hdist=1 ftm=5 qtrim=rl trimq=8 maq=16
$BBDUK-Xmx60g in1=$IN3 in2=$IN4 out1=$IN3_TRIMMED out2=$IN4_TRIMMED ref=$ADAPTERS threads=4 ktrim=r mink=10 hdist=1 ftm=5 qtrim=rl trimq=8 maq=16

#seqkit statistics
SEQKIT=/path/to/seqkit
PARALLEL=/path/to/parallel

ls | grep TRIMMED | $PARALLEL $SEQKIT fx2tab -q -g -l -n {} > {.}_stats.txt


#transcriptome assembly
RNASPADES=/path/to/rnaSPAdes
$RNASPADES --pe-1 1 $IN1_TRIMMED --pe-2 1 $IN2_TRIMMED  --pe-1 2 $IN3_TRIMMED --pe-2 2 $IN4_TRIMMED -t 20 -m 300 -o . 
TRINITY=/path/to/trinity
IN_LEFT=$(zcat $IN1_TRIMMED $IN3_TRIMMED)
IN_RIGHT=$(zcat $IN2_TRIMMED $IN4_TRIMMED)
$TRINITY --seqType fq  --left $IN_LEFT --right $IN_RIGHT  --max_memory 300G --CPU 24
mv trinity_out_dir/Trinity.fasta . && mv transcripts.fasta transcriptome_rnaspades.fasta && mv Trinity.fasta transcriptome_trinity.fasta
RNASPADES_TRANSCRIPTOME=transcriptome_rnaspades.fasta
TRINITY_TRANSCRIPTOME=transcriptome_trinity.fasta

#transrate metrics
TRANSRATE=/path/to/trasnrate
ls | grep transcriptome | $PARALLEL $TRANSRATE --assembly {} --threads 10

#finding lncaRNA
#Arguments: 1 - transcriptome, 2 - mRNA list (fasta, 3 - genome assembly (fasta), 4 - taxonomy table (MEGAN), 5 - protein annotation (gtf)


ls | grep transcriptome | $PARALLEL ./find_lncrna.bash {} $KNOWN_MRNA $GENOME_ASSEMBLY $TAXONOMY_TABLE $PROTEIN_CODING_ANNOTATION

#naming new variables - filtered paf mappings of rnaspades to genome, filtered paf mapping of trinity to genome, filtered exon coordinates of rnaspades, filtered exonic coordinates of trinity, merged exonic coordinates, protein coding annootation withuout bacterial scaffolds
RNASPADES_PAF_FILTERED=transcriptome_rnaspades_filtered_paf
TRINITY_PAF_FILTERED=transcriptome_trinity_filtered_paf
RNASPADES_FILTERED_EXONS=transcriptomes_rnaspades_filtered_before_mapping_to_genome_mpileup_exons_filtered
TRINITY_FILTERED_EXONS=transcriptomes_trinity_filtered_before_mapping_to_genome_mpileup_exons_filtered
TOGETJER_FILTERED_EXONS=$(cat $RNASPADES_FILTERED_EXONS $TRINITY_FILTERED_EXONS)
PROTEIN_CODING_ANNOTATION=protein_coding_annotation_nonbacteria.txt
# finding the lncRNA consensus - R script which reduces transcripts to genes, finds unique genes and between-transcriptome gene overlaps, takes only unique genes, longer genes from overlaps longer than 20% of the length of the longer gene in the pair and longer + shorter genes from overlaps shorter than 20% of the length of the longer gene in the pair
#Arguments: 1 - filtered rnaspades mapping (paf), 2 - filtered trinity mapping (paf)
Rscript find_consensus.R transcriptome_rnaspades_filtered_paf transcriptome_trinity_filtered_paf
#naming new variables - list of rnaspades transcripts in the consensus, list of rnaspades transcripts in the consensus, gene-transcript overlap table for rnaspades, gene-transcript overlap table for trinity, consensus genedatble 
RNASPADES_LIST_PILTERED=transcriptome_rnaspades_filtered_list
TRINITY_LIST_FILTERED=transcriptome_trinity_filtered_list
GENE_TRANSCRIPT_TABLE_RNASPADES=gene_transcript_overlap_table_rnaspades
GENE_TRANSCRIPT_TABLE_TRINITY=gene_transcript_overlap_table_trinity
GENES_TABLE=lncrnas_final_reduced_together_dt
#analysis of the main characteristics of lncRNA (getting exons, introns, element numbers, element lengths, GC content, splice sites)

#Arguments: 1 - rnaspades filtered mapping (paf), 2 - trinity filtered mapping (paf), 3 - rnaspades filtered transcript list, 4 - rnaspades filtered transcript list,  5 - gene transcript overlap table rnaspades, 6 - gene transcript overlap table trinity, 7 rnaspades filtered exons,  8 - trinity filtered exons, argument 9 - genome assembly (fasta), 10 - protein coding annotation with filtered bacterial scaffolds (gtf), 11 - rnaspades transcriptome (fasta), 12 - trinity transcriptome (fasta)
Rscript main_features.R $RNASPADES_PAF_FILTERED $TRINITY_PAF_FILTERED $RNASPADES_LIST_PILTERED $TRINITY_LIST_FILTERED $ISOFORM_TABLE_RNASPADES $ISOFORM_TABLE_TRINITY   $RNASPADES_FILTERED_EXONS $TRINITY_FILTERED_EXONS  $GENOME $PROTEIN_CODING_ANNOTATION $RNASPADES_TRANSCRIPTOME $TRINITY_TRANSCRIPTOME
#Naming new variables - lncrna paf - only isoforms with the biggest number of exons, rnaspades paf - only isoforms with the biggest number of exons,, trinity paf - only isoforms with the biggest number of exons, list of rnaspades longes isoforms, list of trinity longest isoforms, list of protein longest isoforms, table for mapping gene_id to longest isoforms, lncrna intron table
TOGETHER_LONGEST_ISO_PAF=filtered_together_paf_union.txt
RNASPADES_ALL_ISOFORM_PAF=rnaspades_paf_filtered_more_iso.txt
TRINITY_ALL_ISOFORM_PAF=trinity_paf_filtered_more_iso.tx
RNASPADES_LONGEST_ISO_LIST=rnaspades_longest_iso.txt
TRINITY_LONGEST_ISO_LIST=trinity_longest_iso.txt
PROT_LONGEST_ISO=protein_transcripts_longest_iso.txt
LNCRNA_ISOFORM_MAP_TABLE=isoform_map_table_together.txt
LNCRNA_INTRONS_TABLE=introns_table_together.txt
#analysis of the relationship between lncRNA and protein-coding genes - distribution of lengths and numbers of overlapping, intronic and intergenic lncrna, intergenic lncrna distances to closest protein-coding genes, writing protein-coding genes with an overlap with a lncrna-coding gene and protein coding genes 1 kb or closer to the closest lncrna-coding gene (for GO analysis)
#Arguments: 1 - protein coding annotation (filtered for bacterial scaffolds), 2 - lncRNA inron table, 3 - rnaspades paf of all isoforms , trinity paf od all isoforms, filtered exons, isoform-gene overlapping table for rnaspades, isoform-gene overlapping table for trinity, lncRNA gene list, protein longest isoform list
Rscript lncrna_and_prot_relationship.R $PROTEIN_CODING_ANNOTATION $LNCRNA_INTRONS_TABLE $RNASPADES_ALL_ISOFORM_PAF $TRINITY_ALL_ISOFORM_PAF $TOGETJER_FILTERED_EXONS transcriptome_rnaspades.fasta transcriptome_trinity.fasta $ISOFORM_TABLE_RNASPADES $ISOFORM_TABLE_TRINITY $GENES_TABLE $PROT_LONGEST_ISO
#Naming new variables - a table with classification information for every lncRNA transcript (intergenic, intronic, overlapping) + a table with protein annotation according to the relationship with lncRNA genes
LNCRNA_CLASSIFICATION_TABLE=lncRNA_classification_table.txt
PROTEIN_CLASSIFICATION_TABLE=prot_classification_table.txt
#analysis of transposone overlaps - counting the number of "within" overlaps betweem tramsposones and both the inside and the 5 ' region of lncrna-coding genes and prot-coding genes, counting the normalized number of transposone bases which overlap 5' region, exons and introns of those genes (5' region = 1 kb, shorter if the gene is on the begining of a scaffold) + counting the normalized number of bases of transposone classes which overlap 5' region, exons and introns of those genes (5' region = 1 kb, shorter if the gene is on the begining of a scaffold)
#Arguments: 1 - lncrna classification table, 2 - table for mapping gene_id to longest isoforms,  argument 3 - rnaspades filtered exons, argument 4 - trinity filtered exons, argument 5 - paf with every rnaspades isoform, argument 6 - paf with every trinity isoform, argument 7 - repeat annotation, argument 8 - lncrna genes, argument 9 - protein-coding annottaion (filtered for bacterial scaffolds)
#New varaible names lncRNA genes and transposones overlapping table
Rscript transposone_overlaps_analysis.R $LNCRNA_CLASSIFICATION_TABLE $LNCRNA_ISOFORM_MAP_TABLE $RNASPADES_FILTERED_EXONS $TRINITY_FILTERED_EXONS $RNASPADES_ALL_ISOFORM_PAF $TRINITY_ALL_ISOFORM_PAF $REPEAT_ANNOTATION $GENES_TABLE $PROTEIN_CODING_ANNOTATION 
LNCRNA_TRANSPOSONES_OVERLAP_TABLE=lncrna_repeat_overlap_category_20_class_merged.txt
#Expression analysis

#Mapping RNAseq to genome (RNA1 and RNA10)

BBMAP=/path/to/BBmap
REFINDEX=.


$BBMAP path=$REFINDEX in=$IN1_TRIMMED in2=$IN1_TRIMMED out=$OUT.sam t=10 -Xmx100g  ambiguous=random secondary=f maxindel=100000 unpigz=T out=day1_reads_to_genome.sam

$BBMAP path=$REFINDEX in=$IN3_TRIMMED in2=$IN4_TRIMMED out=$OUT.sam t=10 -Xmx100g  ambiguous=random secondary=f maxindel=100000 unpigz=T out=day10_reads_to_genome.sam

SAMBAMBA=/path/to/sambamba

ls | grep reads_to_genome | $PARALLEL $SAMBAMBA view -t 12 -f bam -S -o {.}.bam {}
ls | grep reads_to_genome.bam | $PARALLEL $SAMBAMBA sort -m 80GB --tmpdir tmp -t 12 {}
ls | grep reads_to_genome.sorted.bam | $PARALLEL $SAMTOOLS flagstat ">" {.}_flagstat.txt

#defining the annotation file for FeatureCounts (SAF) - holds lncRNA and protein genes

#Arguments: 1 - protein annotation with filtered bacterial scaffolds (gtf), 2 - lncrna exonic coordinates, 3 - isoform-to-gene mapping table
Rscript make_saf.R $PROTEIN_CODING_ANNOTATION $TOGETJER_FILTERED_EXONS $LNCRNA_ISOFORM_MAP_TABLE

#Counting the mapped reads

FEATURECOUNTS=/path/to/featurecounts

$FEATURECOUNTS -a lncrna_and_prot_annotation.saf -o day1_day10_rna_to_genome_counts.txt  day1_reads_to_genome.sorted.bam day10_reads_to_genome.sorted.bam -F SAF -p  -T 10

#total number of mapped reads for both libraries

ls | grep reads_to_genome.sorted_flagstat.bam | $PARALLEL cat {} | grep -Po ".*(?= \+ 0 mapped)" ">" {.}_total_reads.txt

ls | grep total_reads.txt | xargs cat > total_counts_table.txt
TOTAL_COUNTS_TABLE=total_counts_table.txt

#analyzing the expression - calculating FPKM for every gene, comparing day 1 to day10 expression levels, comparing the expression levels of lncRNAs with a transposone insertions in one of their elements and lncRNAs without transposone insertions, comparing expression levels of protein with and without a relationship to a lncRNA
#Arguments: 1 - counts table for day1 and day10 genes, 2 - total number of mapped reads table, 3 -  lncRNA classification table, 4 - a table which maps longest lncRNA isoform to a gene_id, 5 - a table with transposone and lncRNA overlaps, 6 - paf of lncRNA longest isoforms, 7 - protein classification table based on the relationship with lncRNA genes

Rscript expression_analysis.R day1_day10_rna_to_genome_counts.txt $TOTAL_COUNTS_TABLE $LNCRNA_CLASSIFICATION_TABLE $LNCRNA_ISOFORM_MAP_TABLE $LNCRNA_TRANSPOSONES_OVERLAP_TABLE $TOGETHER_LONGEST_ISO_PAF $PROTEIN_CLASSIFICATION_TABLE

#conservation analysis inside phylum Porifera

#making fasta files of the-biggest-number-of-exons isoforms
$SEQKIT --pattern-file $RNASPADES_LONGEST_ISO_LIST transcriptome_rnaspades.fasta > rnaspades_longest_iso.fa
$SEQKIT --pattern-file $TRINITY_LONGEST_ISO_LIST transcriptome_trinity.fasta > trinity_longest_iso.fa


cat rnaspades_longest_iso.fa trinity_longest_iso.fa > together_longest_iso.fa

#dustmasker the sponge database + perform BLASTn
./blastn_sponge_script.bash together_longest_iso.fa $SPONGE_DATABASE 64 18

#Anylsis of blast hits - produces a blast hit heatmap and a concatanated fasta files of transcript hit coordinates (part trans - only the parts which are similar, whole trans - whole transcripts)
#Arguments: 1 - blast reults
Rscript blast_analysis.R sponges_blast.outfmt6 

#multiple sequence alignment

mkdir msa_dir && mv TRINITY* msa_dir && mv NODE* msa_dir && cd msa_dir

#converting coordinates into fasta files
./get_msa_seq.bash ../together_longest_iso.fa ../${SPONGE_DATABASE} whole_trans

#reversing the rev_comp hits
mkdir for_R && mv *trans_pre.fa for_R && cd for_R
#reversing the rev_comp hits
#Arguments: 1 - concatanated fasta files of transcripts with similarities 
Rscript reverse_sequences.R

#do a MSA 
CLUSTALW2=/path/to/clustalw2

ls | grep _R.fa | $PARALLEL clustalw2 {}

#analyzing conserved secondary structures in the alignment
RNAZ=/path/to/rnaz

ls | grep aln | $PARALLEL {} -o {.}_rnaz_res.txt

#conservation analysis outside of phylum Porifera

NHMMER=/path/to/nhmmer
NHMMER_DATABASE=/path/to/rnacentral/database
NHMMER_QUERY=together_longest_iso.fa
NHMMER_OUTPUT=together_longest_iso_hmmer_out
$NHMMER -tblout  $NHMMER_OUTPUT --rna --qfasta  --tformat fasta --cpu 20 $NHMMER_QUERY $NHMMER_DATABASE > ${NHMMER_OUTPUT}_log.txt

#conservation analysis (outside phylum Porifera) - nhmmer results (rnacentral database without rRNA)
#Arguments: 1 - nhmmer results
Rscript nhmmer_analysis.R $NHMMER_OUTPUT

#infernal analysis - searching for similarities between the RNA queries and Rfam (conserved RNA domains

INFERNAL_CMSCAN=/path/to/infernal_cmscan
INFERNAL_QUERY=together_longest_iso.fa
INFERNAL_OUTPUT=filtered_lncrnas_lonest_trasncripts_infernal_rfam_tblout.txt
CLANIN=/path/to/clanin
RFAM=path/to/rfam
$INFERNAL_CMSCAN --rfam --cut_ga --nohmmonly --tblout $INFERNAL_OUTPUT  --fmt 2 --clanin $CLANIN  $RFAM $INFERNAL_QUERY > $INFERNAL_OUTPUT 
#no significant hits were found - no R analysis

#secondary structure analysis of whole lncRNAs

cd .. && mkdir secondary_structure && cd secondary_structure && cp together_longest_iso.fa .
#generate secondary structures
RNAfold=/path/to/rnafold

$RNAFOLD -i together_longest_iso.fa -p --jobs 10 

#color the secondary structure
RELPLOT=/path/to/relplot

ls | grep -P -o  ".*(?=dp)" | $PARALLEL $RELPLOT -p  {}ss.ps {}dp.ps ">" {}rss.ps

#secondary structure permutation test
RANDFOLD=/path/to/randfold

$RANDFOLD -s together_longest_iso.fa 100 > randfold_together_longest_iso.txt

#generate secondary structures for reverse complements
mkdir revcomp && cp ../together_longest_iso.fa .
REVSEQ=/path/to/revseq

$REVSEQ together_longest_iso.fa -outseq together_longest_iso_revcomp.fa

#generate secondary structures

$RNAFOLD -i together_longest_iso_revcomp.fa -p --jobs 10 

#color the secondary structure

ls | grep -P -o  ".*(?=dp)" | $PARALLEL $RELPLOT -p  {}ss.ps {}dp.ps ">" {}rss.ps

#secondary structure permutation test

$RANDFOLD -s together_longest_iso_revcomp.fa 100 > randfold_together_longest_iso_revcomp.txt

cp randfold_together_longest_iso_revcomp.txt .. && cd ..

#randfold results analysis - take the lower p-value for each sequence-revcomp pair + make a distribution histogram of p-values + compare distributions of p-values across lncRNA classes
#Arguments: 1 - randfold results, 2 - randfold_revcomplement results, 3 - lncRNA classification table
Rscript randfold_analysis.R randfold_together_longest_iso.txt randfold_together_longest_iso_revcomp.txt $LNCRNA_CLASSIFICATION_TABLE