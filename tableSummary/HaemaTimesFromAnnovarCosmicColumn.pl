#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $nc = 12;
my $help;
GetOptions(
	'column=i' => \$nc,
	'help!'    => \$help,
);
&usage if $help;
my $file = shift or &usage; 

open IN, $file or die($!);
while (<IN>){
	chomp;
	my @a = split;
	my $cosmic = $a[ $nc - 1 ];
	if($cosmic eq "."){
		$a[ $nc - 1 ] = 0;
	}elsif($cosmic =~ /\D(\d+)\(haematopoietic_and_lymphoid_tissue\)/){
		$a[ $nc - 1 ] = $1;
	}elsif($cosmic =~ /cosmic|COSMIC|Cosmic/){
		$a[ $nc - 1 ] = "CosmicHaema";
	}else{
		$a[ $nc - 1 ] = 0.5;
	}
	print join "\t", @a;
	print "\n";
}
close IN;

sub usage {
	die(
		qq/
USAGE: HaemaTimesFromAnnovarCosmicColumn.pl [-c ColumnNumber (default:11)] <var.tab.txt>

Note:
	1. This program extracts the frequencies of class haematopoietic_and_lymphoid_tissue from Cosmic database;
	2. The original Cosmic field of the input files should be consistent with ANNOVAR VCF output;
	3. "." means no Cosmic record, and will be converted into 0;
	4. If the mutations were found in diseases other than haematopoietic_and_lymphoid_tissue, the converted value would be 0.5
\n/
	);
}

