#!/usr/bin/perl
use strict;
use warnings;
use File::Copy;
use File::Path;
use Getopt::Long;
use Cwd ;
use POSIX qw(strftime);

#
=head1 NAME

phigis

Primer3 Helper for Indexed Genomes with Ipcress and Samtools

=head1 VERSION

0.1

=head1 DESCRIPTION


This is a wrapper to call shell, samtools and primer3 to make primers
using a samtools indexed genome, dbSNP and RefFlat co-ordinates and gene names.
It reads all variants in a bedfile, finds the corresponding exon in RefFlat,
designs primers that flank the exon by 115 bases, avoiding common SNPs being
withn 8 bases of the 3' end of the primers.
If the exon is >450 bases then only 150 base flanking sequences of the variant are used.
Primers are checked using the ipcress in silico PCR script.

NB make sure the genome builds are consistent

=head2 INPUT

=over

=item *

bed file as text file with unix line endings describing variants

=item *

genome build with samtools index Repeat sequence marked as lower case

=item *

Common SNP file in bed format

=item *

exon file in bed format

=back

=head2 OUTPUT

=over

=item *

log file, including ipcress in silico PCR resulyd

=item *

primer file contained designed primer specifications

=item *

primer3 output file containing additional primers and report on primer design

=back

=head2 OTHER REQUIREMENTS

Bash shell, samtools in path, primer3 in path, exonerate in path, hg19 as samtools indexed fasta, common SNP bed file, exons bed file, perl >=5.12

=head1 EXMAPLE USAGE

./phigis

=head1 AUTHOR

Graham Taylor, Viapath & King's College London

=head1 Licence

MIT Copyright 2017 GR Taylor

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software


=head1 Program Description

=head2 Global Variables

There are probably too many, not all are used but all are commented. It is easier/lazier to pass them between subroutines as globals.

=head2 SETUP

=over 

=item * 

Setup file locations, hard coded but with options to change to STDN

=item *

Takes hard coded input of date

=item *

Setup output and log files

=item *

Setup path to indexed genome fasta

=item *

Read dbSNP, RefFlat and Variant locus bed file

=item *

Setup padding distance (increments until primers designed) 

=back

=head2 MAIN LOOP

For each padding value, try to make primers using primer3

=over

=item *

find_exon_co_ords: find the exon that includes the SNP of interest

=item *

define_target_range: this is the exon + 15 bases + the padding.  The region of interest.

=item *

write_primer3_file: use samtools with backtick system call to extract the ROI as fasta,remove the fasta header then identify any SNPs in the region populate the boulderIO file used by primer3 with the ROI, excluding the exon + 15 bases and exclude SNPs within 8 bases of the 3' end of the primer

=item *

run_primer3 with backtick system call

The values $target_exon, $primer3_template, $primer_space, $exon_plus and exclusion list are interpoloated.
Primer3 is forbidden to use primers with SNPs 8 bases or less from 3'
Lower case (repeat) is excluded within 10 of the 3'

SEQUENCE_ID\=$target_exon
SEQUENCE_TEMPLATE\=$primer3_template
PRIMER_TASK\=pick_pcr_primers
PRIMER_PICK_LEFT_PRIMER\=1
PRIMER_PICK_INTERNAL_OLIGO\=0
PRIMER_PICK_RIGHT_PRIMER\=1
PRIMER_OPT_SIZE\=22
PRIMER_MIN_SIZE\=18
PRIMER_MAX_SIZE\=28
PRIMER_LOWERCASE_MASKING=10
PRIMER_MAX_NS_ACCEPTED\=1
PRIMER_MIN_THREE_PRIME_DISTANCE=3
PRIMER_PRODUCT_SIZE_RANGE\=100\-600
P3_FILE_FLAG\=1
SEQUENCE_TARGET\=$primer_space,$exon_plus 
SEQUENCE_EXCLUDED_REGION=$exclusion_list

=item *

check_SNPs: check location of SNPs witihn ROI

=item *

IF primer3 succeeds product size will be > 0
extract_primers_and_amplicon will write out primers and other details to tsv file

=item *

ELSE add 30 to the padding and try again up to padding = 185.

=item *

Primers are checked using the in silico PCR program ipcress and an ipcress file
is added to the log  The amplicon site should be unique and the product size should
be identical to the predicted product size.

=back

=cut
########################################################################
#global variables
########################################################################
our $date ; #current year, month day hour min
our $variant ; # bed file of variant
our $input_file ; # path to bed file list for variants
our $input ; #hard coded as test.bed, easy to ammend
our @file_list ; #list of variants in $input_file
our $output_dir ; #path to output
our $line ; # generic parsing variable
our @line ; # generic parsing variable
our $genome ; #genome build as fasta indexed using samttols faidx
our $dbSNP ; #name of dbSNP bedfile
our @dbSNP ; #common SNPs slurped into array
our $dbSNP_list ; #intended to be a subset of dbSNP, but in fact just a copy
our $exons ; #hard coded as RefFlat coding exons
our $exons_list ;  #intended to be a subset of RefFlat, but in fact just a copy
our @exons ; #RefFlat slurped into array
our $exon_start; #start of the exon
our $exon_end ;#end of the exon
our $padding ;#distance beyond exon boundaries plus 15 base that is used to select primers
our $primer_space ; #padding minus 30 bases to leave 15 clear bases between primer and ROI
our $exon_size ; #exon size is bases
our $target_exon ; #the exon that includes the variant
our $variant_chr ; #chromosome of the variant (bed)
our $variant_start; #start base of the variant (bed)
our $variant_end ; # end base of the variant (bed)
our $variant_ID ; # 4th column of bed file contains variant ID
our @variant_parse ; #to split variant bed into chr, start, end
our $upstream ; #the distance from the start of the ROI to the upstream SNP;
our $downstream ; #the distance from the start of the ROI to the downstream SNP ;
our $upstream_start ; #the start of the region of interest around the SNP called by samtools for primer3
our $downstream_end ; #the end of the region of interest around the SNP
our $roi_fasta ; #the ROI as fasta, extracted by samtools
our $len_roi_fasta ; #length of the ROI
our $left_flank ;
our $right_flank ;
our $amplicon ;
our @subset_SNP ;  #common SNPs slurped into array not actually subsetted as planned
our @subset_SNP_parse ; # used to split out values from @subset_SNP
our $SNP_start ; # start of SNP locus in bed file
our $SNP_end ; #end of SNP locus in bed file
our @SNP ; #SNPS within the regon of interest
our @subset_exons ; # same as $exons
our @subset_exons_parse ;
our $boilerplate ; #the values passed to samtools to FASTA extraction
our $primer3_template ; #this is the sequence used by primer3 SEQUENCE_TEMPLATE\=
our $primer3_file ; #this is the hard coded primer3 input file
our $fh_primer3_file ; #filehandle for $primer3_file
our $filename ; #used to name files
our $primer3_output ; #used to store the output from primer3
our $exclusion_list ; #loci of SNPs within the ROI.  
our $primer_log ; #this is report of the design process for audit purposes
our $fh_primer_log ; #the is the file handle for the report
our $primer_db ; #this is the primer output as tsv ready for the primer db
our $fh_primer_db ; #file handle for primer output as tsv ready for the primer db
our $fh_fasta ; #file handle for fasta output of amplicon design
our $fasta ; #this file can be pasted into BLAT for a quick view of the primers
our $primer_file_name ; #this is the name given to the primer file
our $primer_name_parse ;# used to extract fields from primer file
our @primer_name_parse ; # used to extract fields from primer file
our $gene ; #extracted field from primer file: gene name
our $cds ;  #extracted field from primer file: coding sequence number
our $left_primer0 ;  #extracted field from primer file: left primer
our $right_primer0 ;  #extracted field from primer file: right primer
our $product_size ; #PCR product size
our $tm ; #annealing/melting temp.
our $ipcress_check ; #ipcress output
########################################################################
#setup time, input and output locations
# to be replaced with getopt eventually
########################################################################
#get date
$date = strftime "%F %R", gmtime ;
$date =~ tr/ \-:/___/ ; #date is year, month, day, hours, minutes
print "date = $date\n"; 
#print "\tEnter path to input file\n";
#$input_file = <STDIN>;
#setup log and report files
#print "\tEnter output directory\n";
#$output_dir = <STDIN> ;
#chomp $output_dir ;
$output_dir = "primer_design_output/" ;
mkpath($output_dir) ;
$primer_log ="$output_dir"."primer_log_$date.txt" ;
$primer_db = "primers_designed\_$date.tsv" ;
open $fh_primer_log, '>', "$primer_log" ;
open ($fh_primer_db, '>', "$primer_db");
print $fh_primer_db "amplicon\tVariant_ID\tgene\tcds\ttarget exon\tTM\tF\tR\tproduct size\tpadding\tdate\n" ;
print "\tEnter path to input file\n\t";
$input_file = <STDIN>;
chomp $input_file;
if (-e $input_file )
	{
	print "\tinput file \= $input_file\n";
	print $fh_primer_log "\tinput file \= $input_file\n";
	}
else
	{	
	print "\t$input_file not found\n";
	die;
	}
print "\toutput directory \= $output_dir\n" if -e $output_dir;
print $fh_primer_log "\toutput directory \= $output_dir\n" if -e $output_dir;
########################################################################
#setup indexed genome (fasta+index), SNP (bed) and exon (bed) sources
# to be replaced with getopt eventually
########################################################################
#print "\tenter path to reference genome fasta\n"; #edit to soft code path to genome
#$genome = <STDIN>;
$genome = "genome/ucsc.hg19.karyotypic.fa"; #relative path in this distro
chomp $genome ;
if (-e $genome )
	{
	print "\tindex reference genome \= $genome \n" ; 
	print $fh_primer_log "\tindex reference genome \= $genome \n" ;
	}
else
	{	
	print "\t$genome not found\n";
	die;
	}
#print "\tenter path to dbSNP bed file\n"; #edit to soft code path to common SNPs
#$dbSNP = <STDIN> ;
$dbSNP = "SNP/SNP_147.bed" ;  #relative path in this distro
chomp $dbSNP ;
if (-e $dbSNP )
	{
	print "\tcommon SNPs \t$dbSNP\n" ;
	print $fh_primer_log "\tcommon SNPs \t$dbSNP";
	}
else
	{	
	print "\t$dbSNP not found\n";
	die;
	}
#print "\tenter path to reflat exons bed file\n"; #edit to soft code path to common exons
#$exons = <STDIN> ;
#relative path in this distro
$exons = "RefSeq/RefFlat_coding_exons.bed";
chomp $exons ;
if (-e $exons )
	{
	print "\texon reference file \= $exons\n" ;
	print $fh_primer_log "\texon reference file \= $exons\n" ;
	}
else
	{	
	print "\t$exons not found\n";
	print $fh_primer_log "\t$exons not found\n";
	die;
	}
########################################################################
#make the files into arrays
########################################################################
open ($input, "<",$input_file) ;
chomp(@file_list = <$input>) ;
open ($exons_list, "<", $exons);
chomp (@exons = <$exons_list>);
open ($dbSNP_list, "<", $dbSNP);
chomp (@dbSNP = <$dbSNP_list>) ;
########################################################################
# main program: loop through the variants to design primers
########################################################################
foreach $variant(@file_list)
	{
	print "\n\n\n\n\t\t\t***** start loop ******\n\t\t\t***********************\n";
	print "\n\n###########################################\n";
	print "$variant\n";
	print "###########################################\n\n";
	if ($variant !~ m/chr\w*\t\d*\t\d*/)
		{
		print "\t $variant not in bed format\n";
		die ;
		}
	print $fh_primer_log "\t$variant\n\n";
	$padding = 35 ; #from 95
	$exon_size = 10 ;
	while ($padding <186 ) #&& $exon_size <= 450)
		{
		undef @SNP ;
		undef @subset_SNP;
		undef @subset_exons;
		$padding = $padding + 30 ;
		$primer_space = $padding -30 ; #from 60
		find_exon_co_ords();
		if ($exon_size >=450 )
			{
			find_genome_co_ords();
			define_target_range();
			write_primer3_file();
			run_primer3();
			check_SNPs() ;
			if ($product_size > 0)
				{
				extract_primers_and_amplicon();
				last ;
				}
			}
		else
			{
			define_target_range();
			write_primer3_file();
			run_primer3();
			check_SNPs() ;
			}
		if ($product_size > 0)
			{
			extract_primers_and_amplicon();
			last ;
			}
		if ($product_size == 0 && $padding > 186)
			{
			retry_as_genomic();
			last ;
			}
		}
	}
exit ;


########################################################################
#subroutines
########################################################################
sub retry_as_genomic
	{
	$padding = 35;
	while ($padding <186)
		{
		$padding = $padding + 30 ;
		$primer_space = $padding -30 ;
		find_genome_co_ords();
		define_target_range();
		write_primer3_file();
		run_primer3();
		check_SNPs() ;
		if ($product_size > 0)
			{
			extract_primers_and_amplicon();
			last ;
			}
		}
	if ($padding >186 && $product_size == 0 )
		{
		print "\n\tfailed for $variant\n" ;
		print $fh_primer_log "\t$variant\n\n";
		print $fh_primer_db "$target_exon\t$variant_ID\tfail\t \t \t \t \t0\t$padding\t$date\n" ;
		}
	}
########################################################################
#find exon co-ordinates
########################################################################
sub find_exon_co_ords
	{
	print "**  sub find_exon_co_ords finding exon co-ordinates for $variant **\n";
	#parse $variant
	@variant_parse = split("\t", $variant) ;
	$variant_chr = $variant_parse[0]; 
	$variant_start = $variant_parse[1] ; 
	$variant_end = $variant_parse[2];
	$variant_ID = $variant_parse[3];
	#check variant is in bed format
	if ($variant_chr !~/chr[X|Y|\d+]/)
		{
			print "\t$variant_chr not bed compliant\n";
			die;
		}
	if ($variant_end < $variant_start )
		{
			print "\t $variant_end not greater than or equal to $variant_start\n";
			die ;
		}
	#subset @dbSNP and @exons by chromosome
	@subset_SNP = @dbSNP ; #not used: was going to select chr-specific lists
	@subset_exons = @exons ; #not used: was going to select chr-specific lists
	#search exons
	foreach $line (@subset_exons)
		{
		@subset_exons_parse = split("\t", $line) ;
		next unless ($subset_exons_parse[0] eq $variant_chr) ;
		$exon_start = $subset_exons_parse[1] ;
		$exon_end = $subset_exons_parse[2];
		$exon_size = ($exon_end - $exon_start) ;
		if ($variant_start > ($exon_start - 15) && $variant_end < ($exon_end + 15))
			{
			$target_exon = "$variant_chr\_$exon_start\-$exon_end" ;
			print "\ttarget exon = $target_exon \n" ;
			print $fh_primer_log "\ttarget exon = $target_exon\n" ;
			#extract_file_name
			$primer_file_name = $line ; #this is used for the primer file name
			#if exon is too large, take 150 base flanks of the SNP
#			if (($exon_size) > 450 || $padding == 185) #or if $exon_start > $variant_end + 15 use_genomic_not_ccds
			last; # may need t move this
			}
		}
	#print "\t$target_exon\n";
	}
########################################################################
#find_genome_co_ords
########################################################################
sub find_genome_co_ords
	{
	print "**  sub find_genome_co_ords finding genome co-ordinates for $variant **\n";
	#parse $variant
	@variant_parse = split("\t", $variant) ;
	$variant_chr = $variant_parse[0]; 
	$variant_start = $variant_parse[1] ; 
	$variant_end = $variant_parse[2];
	$variant_ID = $variant_parse[3];
	#check variant is in bed format
	if ($variant_chr !~/chr[X|Y|\d+]/)
		{
			print "\t$variant_chr not bed compliant\n";
			die;
		}
	if ($variant_end < $variant_start )
		{
			print "\t $variant_end not greater than or equal to $variant_start\n";
			die ;
		}
	#subset @dbSNP and @exons by chromosome
	@subset_SNP = @dbSNP ; #not used: was going to select chr-specific lists
	#ROI no longer chasing splice sites or exons, so can go close to variant
	$exon_start = ($variant_start - 10) ; 	#this is no longer an exon
	$exon_end = ($variant_end + 10) ;		#may rename variable
	$target_exon = "$variant_chr\_$exon_start\-$exon_end" ;
	$exon_size = $exon_end - $exon_start ;
	$target_exon = "$variant_chr\_$exon_start\-$exon_end" ;
	print "\ttarget exon = $target_exon \n" ;
	print $fh_primer_log "\ttarget exon = $target_exon\n" ;
	}
########################################################################
#define target range
########################################################################
sub define_target_range
	{
	my $start ;
	print "** sub define_target_range defining target range **\n";
	#size of exon
	$exon_size = 31 + $exon_end-$exon_start ;
	print "\ttarget exon \+ 15 base flanks size = $exon_size\n";
	print $fh_primer_log "\ttarget exon \+ 15 base flanks size = $exon_size\n";
	print "\tpadding = $padding\n";
	print $fh_primer_log "\tpadding = $padding\n";
	$upstream_start= $exon_start - $padding ; #padding increments to 185
	$downstream_end = $exon_end + $padding ; #padding increments to 185
	$len_roi_fasta = 1 + $downstream_end - $upstream_start ;
	undef $exclusion_list ;
	$exclusion_list ="";
	print "\tprimer3 sequence template length \= $len_roi_fasta\n" ;
	print $fh_primer_log "\tprimer3 sequence template length \= $len_roi_fasta\n" ;
	#get fasta file of the range
	$boilerplate = "$genome $variant_chr:$upstream_start\-$downstream_end" ;
	print "\tsamtools faidx $boilerplate\n";
	print $fh_primer_log "\tsamtools faidx $boilerplate\n";
	#call samtools for fasta
	$roi_fasta = `samtools faidx $boilerplate`;
	print "$roi_fasta\n";
	print $fh_primer_log "$roi_fasta\n" ;
	#search for upstream SNPs
	foreach $line (@subset_SNP)
		{
		@subset_SNP_parse = split("\t", $line) ;
		next unless ($subset_SNP_parse[0] eq $variant_chr) ; 
		$SNP_start = $subset_SNP_parse[1] ;
		#$upstream = $upstream_start - $SNP_start ;
		#$downstream = $SNP_start- $downstream_end ;
		if ($SNP_start > $upstream_start &&  $SNP_start <($exon_start - 15) )
			{
			$upstream = $SNP_start - $upstream_start ;
			push @SNP, "$line\t$exon_start\t$upstream\tbases from start\n" ;
			$upstream = $upstream -8 ;
			$exclusion_list = "$upstream".',8 '."$exclusion_list" ; 
			print "\t$line\t$exon_start\t$upstream\tbases from start\n";
			print $fh_primer_log "\t$line\t$exon_start\t$upstream\tbases from start\n";
			}
		if ($downstream_end > ($SNP_start) && ($exon_end + 15) < $SNP_start )
			{
			$downstream = $SNP_start - $upstream_start ;
			push @SNP, "$line\t$exon_end\t$downstream\tbases downstream\n" ;
			$exclusion_list = "$downstream".',8 '."$exclusion_list" ;
			print "\t$line\t$exon_end\t$downstream\tbases downstream\n";
			print $fh_primer_log "\t$line\t$exon_end\t$downstream\tbases downstream\n";
			}
		}
		#reomve line breaks
		$roi_fasta =~ s/\n//g ;
		#reomve fasta header with substr
		$primer3_template = substr $roi_fasta, -"$len_roi_fasta" ; 
		#determine count and location of lower case for addition to exclusion list
		#make into loop and add to upstream and downstream SNP exclude
		while (my $count = ($primer3_template =~ m/([a|c|g|t|n]*)/g))
			{
			if (length $1 > 0)
				{
				print "count = $count\t$1\n" ;
				}
			my $length = $+[$count] - $-[$count];
			if ($length > 0)
				{
			print "start = $-[$count] finish = $+[$count] length = $length\n";
				if ($length > 10)
					{
					$length = $length - 10 ;
					}
				if (($-[$count] + $upstream_start) < $exon_start)
					{
					$start = ($-[$count] + $upstream_start) ;
					print "$start is less than $exon_start ";
					$exclusion_list = "$-[$count]".','."$length "."$exclusion_list";
					}
				elsif (($-[$count] + $upstream_start) > $exon_end)
					{
					$start = ($-[$count] + $exon_end ) ;
					print "$start is greater than $exon_end ";
					$start = ($-[$count] + 10) ;
					$exclusion_list = "$start".','."$length "."$exclusion_list"; 
					}
				}
			}
		#will this catch multiple?
		##################################################################################
		#to add: loop for multiple lower case instances
		##################################################################################
		print "exclusion list = $exclusion_list\n";
		print $fh_primer_log "exclusion list = $exclusion_list\n";
	}
########################################################################
#write primer3 file
sub write_primer3_file
	{
	print "** sub write_primer3_file writing primer3 file **\n";
	#convert $roi_fasta to primer3 sequencing template
	#reomve line breaks
	#reomve fasta header with substr
	$primer3_template = substr $roi_fasta, -"$len_roi_fasta" ; 
	#Hard coded Primer3 input parameters with $primer3_template interpolated
	my $exon_plus = $exon_size + 30 ;
	$primer3_file = "SEQUENCE_ID\=$target_exon
SEQUENCE_TEMPLATE\=$primer3_template
PRIMER_TASK\=pick_pcr_primers
PRIMER_PICK_LEFT_PRIMER\=1
PRIMER_PICK_INTERNAL_OLIGO\=0
PRIMER_PICK_RIGHT_PRIMER\=1
PRIMER_OPT_SIZE\=22
PRIMER_MIN_SIZE\=18
PRIMER_MAX_SIZE\=28
PRIMER_LOWERCASE_MASKING=1
PRIMER_MAX_NS_ACCEPTED\=1
PRIMER_MIN_THREE_PRIME_DISTANCE=3
PRIMER_PRODUCT_SIZE_RANGE\=100\-600
P3_FILE_FLAG\=1
SEQUENCE_TARGET\=$primer_space,$exon_plus 
SEQUENCE_EXCLUDED_REGION=$exclusion_list
PRIMER_EXPLAIN_FLAG\=1
PRIMER_THERMODYNAMIC_PARAMETERS_PATH\=/usr/local/Cellar/primer3/2\.3\.7/share/primer3/primer3_config/
\=\n";
	#write out primer3 input file
	$filename = "$output_dir"."$target_exon\.txt";
	print "\tfilename = $filename\n";
	open ($fh_primer3_file, '>' , $filename) or die "Could not open file '$filename' $!" ;
	print $fh_primer3_file "$primer3_file" ;
	close ($fh_primer3_file) ;
	}
########################################################################
#run primer3
########################################################################
sub run_primer3
	{
	print "** sub run_primer3 running primer3 **\n";
	#output 
	$primer3_output = `primer3_core <  $filename`;
	$filename = "$output_dir"."$target_exon\_primer3_output\.txt" ;
	print "\toutput for primer3 \= $filename\n";
	open ($fh_primer3_file, '>', $filename) or die "Could not open file '$filename' $!" ;
	if ($filename !~ /PRIMER_LEFT_NUM_RETURNED/)
		{
		$product_size = 0 ;
		}
	print $fh_primer3_file "$primer3_output";
	print $fh_primer_log "$primer3_output";
	close ($fh_primer3_file) ;
	}
########################################################################
#check SNPs
########################################################################
sub check_SNPs
	{
	print "\n\t** checking SNPs **\n";
	#use $upstream and $downstream to locate the SNPs
	#upstream downstream  are in @SNP ;
	my @upstream_SNP = grep (/start/, @SNP) ;
	my @downstream_SNP = grep (/down/, @SNP) ;
	##upstream check
	#take the primer locations from $primer3_output
	print "\tupstream primers\n";
	my @report ;
	my $count = 0 ;
	my $left_primer ;
	my $left_position ;
	my @left_position_split ;
	my @right_position_split ;
	my $right_primer ;
	my $right_position ;
	my $left_primer_3_prime ;
	my $right_primer_3_prime ;
	my @left_primer_3_prime ;
	my @right_primer_3_prime ;
	my $SNP_rs ;
	my $SNP_position ;
	my $SNP_distance ;
	open (my $in, '<', $filename) or die "Could not open file '$filename' $!" ;
	while ($line = <$in>)
		{
		push (@report, $line);
		}
	foreach $line (@report)
		{
		if ($line =~ /PRIMER_LEFT_._SEQUENCE=(\w*)/ )
			{
			if ($count == 0)
				{
				$left_primer0 = $1 ;
				}
			print "left primer $count $1\n"; 
			}
		if ($line =~ /PRIMER_LEFT_.=(\d*,\d*)/ )
			{
			print "left primer position $count $1\n";
			@left_position_split = split /,/, $1 ;
			if ($count == 0)
				{
				$left_flank = $left_position_split[0] + $upstream_start;
				}
			$left_primer_3_prime = ($left_position_split[0] + $left_position_split[1]) ;
			print "left primer  $count 3' location = $left_primer_3_prime\n" ;
			push @left_primer_3_prime , "$count\t$left_primer_3_prime";
			if (scalar @upstream_SNP > 0)
				{
				foreach $line(@upstream_SNP)
					{
					@line = split ("\t", $line) ;
					$SNP_rs = $line[3] ;
					$upstream = $line[7] ;
					$SNP_distance = $left_primer_3_prime - $upstream ;
					print "$count\t$SNP_rs\t$upstream\tdistance from 3'=$SNP_distance\n";
					}
				}
			}
		if ($line =~ /PRIMER_RIGHT_._SEQUENCE=(\w*)/ )
			{
			if ($count == 0)
				{
				$right_primer0 = $1 ;
				}
			print "right primer $count $1\n";
			}
		if ($line =~ /PRIMER_RIGHT_.=(\d*,\d*)/ )
			{
			print "right primer position $count $1\n";
			@right_position_split = split /,/, $1 ;
			if ($count == 0)
				{
				$right_flank = $right_position_split[0] + $upstream_start;
				}
			$right_primer_3_prime = ($right_position_split[0] - $right_position_split[1]) ;
			print "right primer 3' location = $right_primer_3_prime\n" ;
			push @right_primer_3_prime , "$count\t$right_primer_3_prime";
			if (scalar @downstream_SNP > 0)
				{
				foreach $line(@downstream_SNP)
					{
					@line = split ("\t", $line) ;
					$SNP_rs = $line[3] ;
					$downstream = $line[7] ;
					$SNP_distance = $downstream - $right_primer_3_prime ;
					print "$count\t$SNP_rs\t$downstream\tdistance from 3'=$SNP_distance\n";
					}
				}
			}
		if ($line =~ /PRIMER_LEFT_0_TM=(\S*)/)
			{
			$tm = $1 ;
			}
		if ($line =~ /PRIMER_RIGHT_0_TM=(\S*)/)
			{
			if ($1 > $tm)
				{
				$tm = $1 ;
				}
			}
		if ($line =~ /PRIMER_PAIR_._PRODUCT_SIZE=(\d*)/)
			{
			print "product size $count $1\n";
			if ($count == 0)
				{
				$product_size = $1;
				}
			$count ++ ;
			}
		}
	}
########################################################################
#extract primers and amplicon
########################################################################
sub extract_primers_and_amplicon
	{
	#not sure why this lists are written to working directory
	#so moved to output directory using File::Copy
	move "$target_exon\.for" , "$output_dir"."$target_exon\.for" ;
	move "$target_exon\.rev" , "$output_dir"."$target_exon\.rev" ;
	print "\n\t** extracting primers and amplicon **\n";
	@primer_name_parse = split "\t", $primer_file_name ;
	$amplicon = $target_exon ;
	$amplicon =~ s/(.*)(_.*-.*)/$1_$left_flank-$right_flank/ ;
	$primer_file_name = $primer_name_parse[3] ;
	@primer_name_parse = split "_", $primer_file_name ;
	$gene = $primer_name_parse[0] ;
	$cds = $primer_name_parse[2] ;
	print "amplicon:\t$amplicon\n";
	print "Variant_ID:\t$variant_ID\n";
	print "gene:\t$gene\n";
	print "cds:\t$cds\n";
	print "exon co-ords\t$target_exon\n";
	print "melting temp:\t$tm\n" ; 
	print "Forward:\t$left_primer0\n";
	print "Reverse:\t$right_primer0\n";
	print "product size:\t$product_size\n";
	print "padding:\t$padding\n";
	print "date:\t$date\n";
	print $fh_primer_db "$amplicon\t$variant_ID\t$gene\t$cds\t$target_exon\t$tm\t$left_primer0\t$right_primer0\t$product_size\t$padding\t$date\n" ;
	open (my $out, '>', 'ipcress.txt');
	print $out "$variant_ID $left_primer0 $right_primer0 50 700";
	$ipcress_check = `ipcress -i ipcress.txt -m 1 -P 1 -p 0 -s genome/ucsc.hg19.karyotypic.fa` ;
	close $out ;
	print "$ipcress_check\n" ; 
	print $fh_primer_log "$ipcress_check\n";
	#reset variables
	$target_exon = "not defined"; $gene = "X"; $cds = ""; $tm = 0; $left_primer0 = "N";
	$right_primer0 = "N"; $product_size = 0; $upstream_start = 0; $left_flank = 0; $right_flank = 0;
	#tidy files
	}
########################################################################
#END
########################################################################

