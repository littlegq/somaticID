#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
# Author: Qiang Gong <qgong@coh.org>

my $chrColumn = 1;
my $posFormat = "VCF";
my ( $help, $norVar );
$norVar = "./normal95.var.txt.gz";

GetOptions(
    'posFormat=s' => \$posFormat,
    'chrColumn=i' => \$chrColumn,
    'norVar=s'    => \$norVar,
    'help!'       => \$help,
);
&usage if $help;
&usage unless $norVar;
my $file = shift or &usage;

my %nfr;
open IN, '-|', '/usr/bin/gunzip', '-c', $norVar or die($!);
while (<IN>) {
    chomp;
    my @a = split;
    my $k;
    if ( $posFormat eq "VCF" ) {
        $k = $a[1];
    }
    elsif ( $posFormat eq "ANNOVAR" ) {
        $k = $a[0];
    }
    else {
        die(
            "ERROR: Unrecognized posFormat: $posFormat.
			Must be either VCF or ANNOVAR\n"
        );
    }
    $nfr{$k} = $a[2];
}
close IN;

open IN, $file or die($!);
while (<IN>) {
    chomp;
    my @a = split;
    my $k;
    if ( $posFormat eq "VCF" ) {
        $k = join ":", @a[ $chrColumn - 1 .. $chrColumn + 2 ];
    }
    else {
        $k = join ":", @a[ $chrColumn - 1 .. $chrColumn + 3 ];
    }
    my $v = 0;
    $v = $nfr{$k} if exists $nfr{$k};
    print join "\t", @a, "$v\n";
}
close IN;

sub usage {
    die(
        qq/
USAGE: $0 [OPTIONS] -n <normalVariantFile> <inputVariantFile>

OPTIONS:
	-p   STR  Format of genomic position, should be either in VCF or ANNOVAR style, 1-based (Default: VCF)
			VCF-style:
				chr	pos	ref	alt
			ANNOVAR-style:
				SNV:
				chr	pos	pos	ref	alt
				INSERSION:
				chr	pos	pos	-	substr\(alt,1\)
				DELETION:
				chr	pos+1 pos+length\(ref\)-1	substr\(ref,1\)	-
	-c   INT  Number of the chromosome column (Default: 1)
	-nor STR  normal variantion file whith the following column order:
			1:  Position in VCF style \(chr:pos:ref:alt\)
			2:  Position in ANNOVAR style \(chr:pos1:pos2:ref:alt\)
			3:  Frequency observed in normal samples
			4+: Other information such as dbSNP.
	-h      Show this page

CONTACT: Qiang Gong <qgong\@coh.org>
\n/
    );
}

