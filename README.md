# phigis
PCR primer design
NAME
       phigis

       Primer3 Helper for Indexed Genomes with Ipcress and Samtools

VERSION
       0.1

DESCRIPTION
       This is a wrapper to call shell, samtools and primer3 to make primers
       using a samtools indexed genome, dbSNP and RefFlat co-ordinates and
       gene names.  It reads all variants in a bedfile, finds the
       corresponding exon in RefFlat, designs primers that flank the exon by
       115 bases, avoiding common SNPs being withn 8 bases of the 3' end of
       the primers.  If the exon is >450 bases then only 150 base flanking
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
       o   log file, including ipcress in silico PCR resulyd

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

