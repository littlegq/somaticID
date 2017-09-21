#!/usr/bin/perl
use strict;
use warnings;


my $pfx = shift or die("Usage: $0 <prefix>
	Note: Only works on *.vcf from somaticID.pl
	      Only report non-silent mutations
		  Exclude common SNPs (>1%)\n");

my @header = qw(Chr Pos Ref Alt Gene AAC MutType dbSNP PopFreq CADD_phred Comic ClinVar MutGroup Somatic_Prob);
foreach my $f ("$pfx.somatic.snp.vcf", "$pfx.somatic.indel.vcf"){
	&vartable($f, "Somatic");
}
foreach my $f ("$pfx.non-somatic.snp.vcf", "$pfx.non-somatic.indel.vcf"){
	&vartable($f, "Non-Somatic");
}
foreach my $f ("$pfx.LowQualityMutations.snp.vcf"){
	&vartable($f, "LowQuality");
}

sub vartable {
	my $mut_group = $_[1];
	open IN, $_[0] or die($!);
	while (<IN>){
		chomp;
		next if /^#/;
		my @a = split /\t/;
		my $n_sample = $#a - 8;
		if( @header == 14 ){
			for (my $i=1; $i<= $n_sample; $i++){
				push @header, ("Sample$i.Ref","Sample$i.Var","Sample$i.VAF");
			}
			print join "\t", @header;
			print "\n";
		}
		print STDERR "WARNING: Unexpected column number:\n$_\n" 
			if $#header - 13 < $n_sample * 3;
		my @s = split /;/, $a[7];
		foreach my $site (@s){
			if($site=~/=(\S+)$/){
				$site = $1;
			}
		}
		$s[9] = $s[6] if $s[9] eq ".";
		my $gene = $s[7];
		my $aac = ".";
		my $som_prob = ".";
		$som_prob = $s[55] if @s == 56;
		unless($s[10] eq "."){
			my @pcs = split /,/, $s[10];
			my $aap = 0;
			foreach my $pc (@pcs){
				my @p = split /:/, $pc;
				if($pc =~ /:p\.\D*(\d+)\D/){
					my $cp = $1;
					next unless $cp > $aap;
					$aap = $cp;
					$aac = $p[4];
					$gene = $p[0];
				}else{
					print STDERR "WARNING: Cannot get amino acid changing information from $pc\n";
					$aac = $pc;
				}
			}
		}
		next if $gene =~ /^LOC|^ENS|^FAM|^GOL|^PRA|^NBP|^POT|^DEF|^OR2|^MUC|^KRT|^WAS|^ANK|^TRI|^OR1|^FRG/;
		my @rdout;
		foreach my $rd (@a[9..$#a]){
			my @ds = split /:/, $rd;
			if(@ds == 14){
				push @rdout, join "\t", @ds[4..6];
			}else{
				push @rdout, ".\t.\t.";
			}
		}
		print join "\t", @a[0,1,3,4], $gene, $aac, @s[9,16,12,35,15,14], $mut_group, $som_prob, @rdout;
		print "\n";
	}
	close IN;
}






		
