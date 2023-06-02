#!/usr/bin/env perl
use strict;
use warnings;

use Data::Dumper;
use Getopt::Long;
use List::MoreUtils qw/all mesh/;

use FindBin qw( $RealBin );
use lib  "$RealBin/";

my $inputFASTQ;
my $inputFASTQ_R1;
my $inputFASTQ_R2;
my $outputFASTQ;
my $outputFASTQ_R1;
my $outputFASTQ_R2;
my $samtools_bin;
my $minimap2_bin;
my $GRCh38_plus_SARSCoV2_fasta;
my $threads = 4;
GetOptions (
	'inputFASTQ:s' => \$inputFASTQ,
	'inputFASTQ_R1:s' => \$inputFASTQ_R1,
	'inputFASTQ_R2:s' => \$inputFASTQ_R2,
	'outputFASTQ:s' => \$outputFASTQ,
	'outputFASTQ_R1:s' => \$outputFASTQ_R1,
	'outputFASTQ_R2:s' => \$outputFASTQ_R2,
	'samtools_bin:s' => \$samtools_bin, 	
	'minimap2_bin:s' => \$minimap2_bin, 
	'GRCh38_plus_SARSCoV2_fasta:s' => \$GRCh38_plus_SARSCoV2_fasta, 	
	'threads:s' => \$threads, 	
);

die "Please provide either --inputFASTQ (for Nanopore data) or --inputFASTQ_R1 and --inputFASTQ_R2 (for Illumina data)" unless($inputFASTQ xor ($inputFASTQ_R1 or $inputFASTQ_R2));

if($inputFASTQ)
{
	die "Please provide existing --inputFASTQ" unless($inputFASTQ and -e $inputFASTQ);
	die "Please do not provide a gzipped input FASTQ (--inputFASTQ)" if($inputFASTQ =~ /\.gz$/);
	die "Please provide --outputFASTQ" unless($outputFASTQ);
}

if($inputFASTQ_R1 or $inputFASTQ_R2)
{
	die "Please provide existing --inputFASTQ_R1" unless($inputFASTQ_R1 and -e $inputFASTQ_R1);
	die "Please provide existing --inputFASTQ_R2" unless($inputFASTQ_R2 and -e $inputFASTQ_R2);
	die "Please do not provide a gzipped input FASTQ (--inputFASTQ_R1)" if($inputFASTQ_R1 =~ /\.gz$/);
	die "Please do not provide a gzipped input FASTQ (--inputFASTQ_R2)" if($inputFASTQ_R2 =~ /\.gz$/);
	die "Please provide --outputFASTQ_R1" unless($outputFASTQ_R1);
	die "Please provide --outputFASTQ_R2" unless($outputFASTQ_R2);
}

$samtools_bin = find_path('samtools_bin', $samtools_bin, 'samtools');
$minimap2_bin = find_path('minimap2_bin', $minimap2_bin, 'minimap2');

$GRCh38_plus_SARSCoV2_fasta = $RealBin . '/GRCh38_full_analysis_set_plus_decoy_hla_and_SARS-CoV-2.fa' unless(defined $GRCh38_plus_SARSCoV2_fasta);

unless(-e $GRCh38_plus_SARSCoV2_fasta)
{
	die "Path to file GRCh38_full_analysis_set_plus_decoy_hla_and_SARS.fa could not be found - please put either in the same directory as the file filterHumanGenome.pl, or explicitly specify a path via parameter --GRCh38_plus_SARSCoV2_fasta. See GitHub for details."
}

unless((-e $minimap2_bin) and (-x $minimap2_bin))
{
	die "Could not determine the path to minimap2 binary - candidate: '$minimap2_bin', but does not seem to work. This should not happen, please contact the authors of the program.";
}
unless((-e $samtools_bin) and (-x $samtools_bin))
{
	die "Could not determine the path to samtools binary - candidate: '$samtools_bin', but does not seem to work. This should not happen, please contact the authors of the program.";
}

print STDERR "filterHumanGenome.pl:\n";
print STDERR "\t", "samtools path: $samtools_bin\n";
print STDERR "\t", "minimap2 path: $minimap2_bin\n";
print STDERR "\n";


my $minimap2_cmd;
if($inputFASTQ)
{
	print STDERR "Now processing file '$inputFASTQ' (Nanopore mode)\n";
	$minimap2_cmd = qq($minimap2_bin -t $threads -a -x map-ont $GRCh38_plus_SARSCoV2_fasta $inputFASTQ | $samtools_bin view -F 0x4 -);	
}
else
{
	print STDERR "Now processing file '$inputFASTQ_R1' and '$inputFASTQ_R2' (Illumina paired-end mode)\n";
	$minimap2_cmd = qq($minimap2_bin -t $threads -a -x sr $GRCh38_plus_SARSCoV2_fasta $inputFASTQ_R1 $inputFASTQ_R2 | $samtools_bin view -F 0x4 -);	
}
	
my %reads_viral;
my %reads_potentially_human;
open(MINIMAP2, "$minimap2_cmd |") or die "Cannot open pipe to $minimap2_cmd";
while(<MINIMAP2>)
{
	my $line = $_;
	chomp($line);
	my @fields = split(/\t/, $line);
	my $readID = $fields[0];
	my $reference = $fields[2];
	if($reference eq 'MN908947.3')
	{
		$reads_viral{$readID}++;
	}
	else
	{
		$reads_potentially_human{$readID}++;
	}
}
close(MINIMAP2);

my $reads_removed_potentiallyHuman = 0;
my $reads_removed_notAlignedToViralGenome = 0;
my $reads_notRemoved = 0;
my $bases_removed_potentiallyHuman = 0;
my $bases_removed_notAlignedToViralGenome = 0;
my $bases_notRemoved = 0;
my @filterReadTuples;
if($inputFASTQ)
{
	push(@filterReadTuples, [$inputFASTQ, $outputFASTQ]);
}
else
{
	push(@filterReadTuples, [$inputFASTQ_R1, $outputFASTQ_R1]);
	push(@filterReadTuples, [$inputFASTQ_R2, $outputFASTQ_R2]);
}
foreach my $filterReadTuple (@filterReadTuples)
{
	open(FASTQIN, '<', $filterReadTuple->[0]) or die "Cannot open $filterReadTuple->[0]";
	open(FASTQOUT, '>', $filterReadTuple->[1]) or die "Cannot open $filterReadTuple->[1]";
	while(<FASTQIN>)
	{
		my $readID_line = $_;
		chomp($readID_line);
		if($readID_line)
		{
			die "Weird - read ID should begin with '\@' - line $. of $filterReadTuple->[0] - '$readID_line'" unless(substr($readID_line, 0, 1) eq '@');
			my $readID = substr($readID_line, 1);
			$readID =~ s/\s.+//;
			my $read_sequence = <FASTQIN>;
			my $plus = <FASTQIN>;
			my $read_qualities = <FASTQIN>;
			die "FASTQ corruption around line $. of $filterReadTuple->[0] - expect a plus character" unless(substr($plus, 0, 1) eq '+');
			die "FASTQ corruption around line $. of $filterReadTuple->[0] - expect length(read) == length(qualities)" unless(length($read_sequence) == length($read_qualities));
			if(exists $reads_potentially_human{$readID})
			{
				$reads_removed_potentiallyHuman++;
				$bases_removed_potentiallyHuman += length($read_sequence);
			}
			else
			{
				if($reads_viral{$readID})
				{
					$reads_notRemoved++;
					$bases_notRemoved += length($read_sequence);
					print FASTQOUT $readID_line, "\n", $read_sequence, $plus, $read_qualities;
				}
				else
				{
					$reads_removed_notAlignedToViralGenome++;
					$bases_removed_notAlignedToViralGenome += length($read_sequence);
				}
			}
		}
	}
	close(FASTQIN);
	close(FASTQOUT);
}
print STDERR "Done. $reads_notRemoved remaining reads ($reads_removed_potentiallyHuman removed as potentially human, $reads_removed_notAlignedToViralGenome removed as not aligning to viral genome)\n";

print join("\t", $reads_removed_potentiallyHuman, $reads_removed_notAlignedToViralGenome, $reads_notRemoved, $bases_removed_potentiallyHuman, $bases_removed_notAlignedToViralGenome, $bases_notRemoved), "\n";

sub find_path
{
	my $id = shift;
	my $supplied_value = shift;
	my $forWhich = shift;
	
	my $this_bin_dir = $FindBin::RealBin;
	my %paths_ini;
	my $paths_ini = $this_bin_dir . '/paths.ini';
	if(-e $paths_ini)
	{
		open(INI, '<', $paths_ini) or die "Cannot open $paths_ini";
		while(<INI>)
		{
			chomp;
			next unless($_);
			$_ =~ s/[\n\r]//g;
			next if($_ =~ /^\s+$/);
			die "Invalid line format in $paths_ini -- expect all lines to be either empty or key=value pairs" unless($_ =~ /^(.+)=(.*)$/);
			my $id = $1;
			my @alts = split(/,/, $2);
			$paths_ini{$id} = \@alts;
		}
		close(INI);
	}
		
	if(defined $supplied_value)
	{
		die "Command-line supplied value/file for parameter $id not existing" unless(-e $supplied_value);
		return $supplied_value;
	}
	
	if(exists $paths_ini{$id})
	{
		foreach my $alternative (@{$paths_ini{$id}})
		{
			if(-e $alternative)
			{
				return $alternative;
			}
		}
	}	

	if($forWhich)
	{
		my $which_output = `which $forWhich`;
		$which_output =~ s/[\n\r]//g;
		if($which_output and (-e $which_output))
		{
			return $which_output;
		}
	}
	
	die "I couldn't figure out a path for ${id}. Order of precedence: check for parameter --${id}; check paths.ini in $this_bin_dir for a key named $id; parse, if command string defined, the output of the command 'which ${forWhich}' ('which ' means that the command string is not defined).";
}

