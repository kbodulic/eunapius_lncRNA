# eunapius_lncRNA
R and bash code for the Master's Thesis "Computational analysis of long non-coding RNA in endemic cave sponge (Eunapius subterraneus)", K. Bodulić, 2020


##File description:

Master script: Master.bash - executes every step of the analysis



##lncRNA extraction pipeline:

filter_rrna_and_length.R - filters rRNA transcripts + transcripts shorter than 200 bp

filter_orfs.R - filters transcripts with the longest ORF longer than 150 bp

filter_diamond.R - filters transcripts with DIAMOND hits to known proteins

filter_hmmer.R - filters transcripts with HMMER hits to known protein domains

mapping_filtering.R - filters transcripts based on mapping to genome (removes transcripts which map to bacterial scaffolds, transcripts with one exon, transcripts with mapping identity which is less than 95% and transcripts whose exons overlap protein-coding exons)

find_consensus.R - finds a miminal-overlapping consensus of lncRNA assembled with rnaSPAdes and Trinity


#lncRNA analysis scripts:

main_features.R - analyses the main features of found lncRNA-coding genes comparing them to protein-coding genes(isoform numbers, exon numbers, exon, intron, transcript and gene lengths, splice sites, GC content)

lncrna_and_prot_relationship.R - analysis the relationship of protein-coding genes and lncRNA coding genes. distribution of lengths and numbers of overlapping, intronic and intergenic lncrna, intergenic lncrna distances to closest protein-coding genes, writing protein-coding genes with an overlap with a lncrna-coding gene and protein coding genes 1 kb or closer to the closest lncrna-coding gene (for the GO analysis)

transposone_overlap.analysis.R - analysis of lncRNA and transposone relationship (analysis of transposone overlaps - counting the number of "within" overlaps betweem 
tramsposones and both the inside and the 5 ' region of lncrna-coding genes and prot-coding genes, counting the normalized number of transposone bases which overlap 5' region, exons and introns of those genes (5' region = 1 kb, shorter if the gene is on the begining of a scaffold) + counting the normalized number of bases of transposone classes which overlap 5' region, exons and introns of those genes (5' region = 1 kb, shorter if the gene is on the begining of a scaffol)

expression_analysis.R - analyses the lncRNA-coding gene expression (calculating FPKM for every gene, comparing day 1 to day10 expression levels, comparing the expression levels of lncRNAs with a transposone insertions in one of their elements and lncRNAs without transposone insertions, comparing expression levels of protein with and without a relationship to a lncRNA)

blastn_sponge_script.bash - constructs a database of all available Porifera datasetsm peforms BLAST on the lncRNA query sequences and the constructed Porifera database

blast_analysis.R - performs the analysis of lncRNA BLAST hits in the Porifera phylum -  produces a blast hit heatmap and a concatanated fasta file of  transcript hit coordinates

nhmmer_analysis.R - performs the analysis of lmncRNA nHMMER hits outside of the Porifera phylum

randfold_analysis.R - analysis the 2D structures of found lncRNAs (takes the lower p-value for each sequence-revcomp pair + makes a distribution histogram of p-values + compares distributions of p-values across lncRNA classes)



##Helper scripts:

exons.R - finds the exonic coordinates of genome-mapped transcripts

make_saf.R - makes the SAF annotation format of lncRNA-coding and protein-coding genes

get_msa_seq.bash - makes a multiple alignmen file of the Porifera lncRNA BLAST hits

reverse_sequences.R - makes a reverse complement of reverse complement BLAST hits - important for multiple alignment

rename.pl - renames the fasta headers according to header list text file
