# phigis
PCR primer design

Exome and genome-scale genetic diagnostics has increased the demand for PCR confirmations in the proband and for prenatal and cascade testing in relatives. Genetics Service Exome-based tests, identifying tens to hundreds of targets per week, many of which require the design of cognate PCR amplicons for further testing in other at risk individuals or for prenatal testing. Existing tools are not readily scalable for high throughput primer design: it takes 15 to 30 minutes to use web-driven tools like Primer 3 to design each primer pair, which must then be checked using SNPchecker, an application that is no longer supported. Other high throughput designer tools like SNP Box and Primer Mapper are not optimised for high-throughput, exon-focussed work using human genome reference tools. Simply feed phigis a list of variants in bed file format and primers will be designed in a few minutes using commodity hardware such as a Macbook.

NAME
       phigis

       Primer3 Helper for Indexed Genomes with Ipcress and Samtools

VERSION
       0.1

DESCRIPTION
       This is a wrapper to call shell, samtools and primer3 to make primers
       using a samtools indexed genome, dbSNP and RefFlat co-ordinates and
       gene names.  It reads all variants in a bedfile, finds the
       corresponding exon in RefFlat, designs primers that flank the exon, 
       avoiding common SNPs being withn 8 bases of the 3' end of
       the primers.  If the exon is >450 bases then only 60 base flanking
       sequences of the variant are used.  Primers are checked using the
       ipcress in silico PCR script.

       NB make sure the genome builds are consistent

INPUT

       o   bed file as text file with unix line endings describing variants

       o   genome build with samtools index Repeat sequence marked as lower
           case

       o   Common SNP file in bed format

       o   exon file in bed format

   OUTPUT
   
       o   log file, including ipcress in silico PCR result

       o   primer file contained designed primer specifications

       o   primer3 output file containing additional primers and report on
           primer design

   OTHER REQUIREMENTS
       Bash shell, samtools in path, primer3 in path, exonerate (EBI) in path, hg19
       as samtools indexed fasta, common SNP bed file, exons bed file, Perl
       >=5.12
       primer3 uses themodynamic parameters stored somewhere like
       /usr/local/Cellar/primer3/2.3.7/share/primer3/primer3_config/
       The path may vary depending how you installed primer3, and so may the version number
       
 Additional relative file paths not included
       geneme/genome_fasta.fa : the fasta genome file with lower case repeats
       geneme/genome_fasta.fa.fai : the samtools indexed fasta file
       SNP/SNP_147.bed : the SNP bed file
       RefSeq/RefFlat_coding_exons.bed : the exons file
       test.bed : the bed file describing the variant needing primer design
       
References

Bioinformatics. 2009 25:2078-9. doi: 10.1093/bioinformatics/btp352. Epub 2009 The Sequence Alignment/Map format and SAMtools Li H et al.

Nucleic Acids Res. 2012 e115 Primer3—new capabilities and interfaces Untergasser,A et al.

BMC Bioinformatics 20056:31 Automated generation of heuristics for biological sequence comparison Slater G & Birney E

       
       

