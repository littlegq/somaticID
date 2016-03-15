#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Std;

our ($opt_f,$opt_h,$opt_d);
$opt_d = 1 ;
getopt("f:hd:");
&usage if ($opt_h);
&usage unless ($ARGV[0]);

my ($last,$lastL)=qw(0 0);
my $fie = $opt_f || 2 ;
$fie -- ;

open IN, $ARGV[0];
my ($l,@a) ;
$l = <IN>;
@a = split/\s+/,$l;
&upd();

while(1) {
    if ( $l !~ /^#/ && abs($a[$fie]-$last) <= $opt_d ) {
        while (abs($a[$fie] - $last) <= $opt_d ) {
            &upd();
        }
            &upd();
    }
    else {
        print "$lastL\n";
        &upd();
    }
}

sub upd (){
    if (eof(IN)) {
		unless ($l=~/^#/){
			print "$l" if ($a[$fie] - $last > $opt_d );
		}
        last;
    }
	if($l=~/^#/){
		print "##FILTER=<ID=adjcent$opt_d,Description=\"Mutations that locate within $opt_d bp\">\n" 
			if $l =~ /^#CHROM/;
	}else{
	    $last = $a[$fie];
	}
	$lastL = $l;
	chomp $lastL;
	$l = <IN>;
	@a = split/\s+/,$l unless $l=~/^#/;
}

sub usage (){
    die qq(
#===============================================================================
#
#        USAGE:  ./rmAdjError.pl  <STDIN>  
#                -f : filed of position, start from 1, default :2 
#                -d : addjacent distance, default : 1
#                -h : display this hlep
#
#  DESCRIPTION:  Remove all adjacent error from sorted VCF file
#
#       AUTHOR:  Wang yu , wangyu.big at gmail.com
#    INSTITUTE:  BIG.CAS
#      CREATED:  04/15/2010 
#      Modified by Qiang Gong on 2/10/2016
#===============================================================================

)
}
