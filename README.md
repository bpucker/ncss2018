
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2586989.svg)](https://doi.org/10.5281/zenodo.2586989)


# Script collection for the investigation of splice sites

Summary:

These scripts provide functions to analyse splice sites based on a genome sequence (FASTA) and a corresponding annotation (GFF3) provided by the NCBI. Next, RNA-Seq reads can be mapped via STAR and the resulting BAM files can be converted into customized coverage files (COV) which are required for the validation of splice sites. The usage of individual splice sites can be quantified based on the number of mapped and propperly splitted RNA-Seq reads. These scripts generate report files and figures for the documentation of the analysis process.

Detailed instructions:

The following points describe the processing of data for the investigation of splice site combinations based on the NCBI reference sequence and annotation for a given species. Python 2.7.x is required for the execution of these scripts. Some scripts require additional modules including matplotlib, numpy, and scipy. The whole workflow was developed and run on a LINUX system.

•	Genome sequences (FASTA, .fna) and corresponding annotations (GFF3, .gff) were manually downloaded from the NCBI (https://www.ncbi.nlm.nih.gov/genome/browse#!/overview/). Files were named according to species with underlines separating different parts of the name e.g. ‘Homo_sapiens’.

•	Run NCBI_splice_site_check.py on the folder with all genome sequences and corresponding annotations. This script will generate many files necessary for downstream analysis including a text document which lists all splice sites in representative transcripts, a multiple FASTA file containing all representative transcripts, and a FASTA file containing all representative peptide. Representative transcripts/peptides are defined as the combination of exons, which encodes the longest peptide.

•	Run splice_site_divergence_check.py analyses all data prepared in the previous step (--data_dir). Results are written into a single output folder (--output_dir) to allow manual inspection of figures.  Some analyses are designed to reveal a phylogenetic if present. Therefore, species should be supplied based on their order in a tree (--species_file) and not just in alphabetical order. However, this step can be performed without any prior knowledge about the phylogenetic relationships of the involved species.

•	The application of run_BUSCO_across_specs.py assesses the completeness of peptide sets in one folder and can be run independently of most other scripts. This script represents a wrapper for BUSCO v3 in protein mode. It is important to supply an appropriate benchmarking data set to BUSCO which is actual representative for the taxonomic group of interest. It would be easy to modify this script for the assessment of transcript set completeness.

•	RNA-Seq data sets were manually selected at the SRA website (https://www.ncbi.nlm.nih.gov/sra) by selecting the appropriate species and ‘RNA’ as source. If possible, data sets from different tissues/conditions/developmental stages were selected to allow the investigation of many different transcripts. The Run selector was applied for additional filtering to restrict the selection to paired-end data sets.

•	Fastq-dump (sratoolkit v2.9.6-1) was used with the options ‘--skip-technical --read-filter pass --dumpbase --split-3 --clip’ to download and process files of all selected runs.

•	Mapping of the RNA-Seq reads was performed with STAR v2.5.1b wrapped into customized Python scripts to allow automatic processing of all samples on our compute cluster. Data preparation with STAR used the following parameters: “--genomeChrBinNbits 10 --genomeSAsparse 2 --runThreadN 7 --limitGenomeGenerateRAM 70000000000 --genomeSAindexNbases 4”. The actual mapping was performed in the recommended 2-pass mode using the following parameters: “--runThreadN 13 --limitBAMsortRAM 60000000000 --outBAMsortingThreadN 1 --outSAMtype BAM SortedByCoordinate --twopassMode Basic –outFilterMismatchNoverLmax 0.05 --outFilterMatchNminOverLread 0.9”. featureCounts v.1.6.3 was applied for quantification of gene expression based on features annotated in the GFF3 file. Expression was summarized on the gene level while excluding multi mapped reads.

•	Run combine_results_of_all_count_tables.py to merge the expression values from different RNA-Seq samples in a single file with multiple columns. While generating this combined file of raw counts, TPMs and FPKMs are also calculated if a suitable GFF3 file is provided to extract transcript lengths.

•	Different RNA-Seq samples of the same species were mapped in parallel to increase speed, final BAM files were combined via the merge function of samtools. Only the IDs of samples which resulted in valid BAM files were carried on to downstream analyses. This prevented the inclusion of samples which might be mis-labeled at the SRA.

•	The final BAM file was processed with construct_RNA_seq_coverage_file.py was to generate COV files for the validation of annotated splice sites. An installation of bedtools is required for this script. The resulting COV file contains the number of aligned reads per position in the genome sequence.

•	The actual validation of annotated splice sites was performed by check_ncss.py based on the COV file and the previously generated TXT file describing all splice site combinations in the representative transcript. RNA-Seq reads are expected to span introns when aligned to the genome sequence i.e. the intron results in a gap in the alignment. Therefore, the coverage of aligned bases should drop when moving from an exon into an intron. Each splice site is checked for such a drop in the coverage between the two distal exon positions and the two distal intron positions. If at least three reads span this position and the coverage drops by at least 20%, a splice site is considered to be validated.

•	To quantify the total number of aligned RNA-Seq read bases get_RNA_seq_amount.py was applied which calculates this value based on COV files.

•	The initial splice site validation results are summarized and processed by analyze_supported_splite_sites.py. This script also assesses the correlation of available RNA-Seq information with the number of validated splice sites. Therefore, additional information about the number of total bases included in the RNA-Seq data sets is required.

•	In order to reveal specific properties of introns with non-canonical splice site combinations, the lengths of different intron types were compared. Information about introns in representative transcripts is loaded from the initially generated TXT file.




# References
Frey K., Pucker, B. (2019). Animal, fungi, and plant genome sequences harbour different non-canonical splice sites. bioRxiv 616565. https://www.biorxiv.org/content/10.1101/616565v1.full.


Boas Pucker and Samuel Fraser Brockington (2018). Genome-wide analyses supported by RNA-Seq reveal non-canonical splice sites in plant genomes. BMC Genomics. 2018;19(1). https://doi.org/10.1186/s12864-018-5360-z.


Boas Pucker, Daniela Holtgräwe, Bernd Weisshaar (2017). Consideration of non-canonical splice sites improves gene prediction on the Arabidopsis thaliana Niederzenz-1 genome sequence. BMC Research Notes, 10, 667. https://doi.org/10.1186/s13104-017-2985-y. 
