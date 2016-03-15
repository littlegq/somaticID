#!/usr/bin/perl
use strict;
use warnings;

## Author: Qiang Gong <gongqiang.big@gmail.com>
## Transfer ANNOVAR-annoated VCF files into tab-deliminated files


my $file = shift or die("Usage: $0 <*.nonsilent.vcf>\n");

open IN, $file or die($!);
while (<IN>){
	chomp;
	my @a = split;
	
