#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

## This is a pipeline to identify somatic mutations when there is no matched normal data
## Author: Qiang Gong <gongqiang.big@gmail.com>

## Initialize and obtain arguments
my $pfx;
my $nbam;
my $samtools      = "samtools";
my $varscan       = "~/bin/VarScan.v2.3.6.jar";
my $somaticID_dir = $0;
$somaticID_dir =~ s/[^\/]+$//;
$somaticID_dir = "." unless $somaticID_dir;
my $vcf_filter        = "$somaticID_dir/varscan_vcf_filter.pl";
my $non_silent_filter = "$somaticID_dir/nonsilent_filter_for_annovarVCF.pl";
my $rmAdj             = "$somaticID_dir/rm_adj_error_VCF.pl";
my $rfc               = "$somaticID_dir/RF_classifier.R";
my $mod_dir           = "$somaticID_dir/RDA";
die( "Could not find directory $somaticID_dir/RDA\n" ) unless -d $mod_dir;

## Check necessary files
die( "Could not find RF_classifier.R in the directory $somaticID_dir\n" )
  unless -f $rfc;
die(
"Could not find executable nonsilent_filter_for_annovarVCF.pl in the directory $somaticID_dir\n"
) unless -x $non_silent_filter;
die(
"Could not find executable rm_adj_error_VCF.pl in the directory $somaticID_dir\n"
) unless -x $rmAdj;

my $java       = "java";
my $rscript    = "Rscript";
my $tmpdir     = "/tmp";
my $memory     = "16G";
my $mkfifo     = "/usr/bin/mkfifo";
my $annovardir = "~/bin/annovar";
my $annovar_protocol =
"refGene,cytoBand,popfreq_max_20150413,nci60,clinvar_20150629,cosmic70,snp138,dbnsfp30a,dbnsfp31a_interpro,dbscsnv11";
my $annovar_operation = "g,r,f,f,f,f,f,f,f,f";
my $buildver          = "hg19";
my $ref;
my ( $novar, $noanno, $nopred );
my $adj_d = 0;
my $selfmodel;
my $help;
GetOptions(
    "normal=s"     => \$nbam,
    "output=s"     => \$pfx,
    "reference=s"  => \$ref,
    "samtools=s"   => \$samtools,
    "varscan=s"    => \$varscan,
    "annovardir=s" => \$annovardir,
    "buildver=s"   => \$buildver,
    "java=s"       => \$java,
    "rscript=s"    => \$rscript,
    "tmpjava=s"    => \$tmpdir,
    "memory=s"     => \$memory,
    "mkfifo=s"     => \$mkfifo,
    "novarcall!"   => \$novar,
    "noanno!"      => \$noanno,
    "nopred!"      => \$nopred,
    "selfmodel!"   => \$selfmodel,
    "rmAdjErr=i"   => \$adj_d,
    "help!"        => \$help
) or &Usage;

&Usage if $help;
my $n_tum = $#ARGV + 1;
&Usage unless $n_tum > 0;
my $tbams = join " ", @ARGV;
my $ltime = localtime();
print "[$ltime] There are $n_tum tumors.\n" if $nbam;
print "[$ltime] We have normal data!\n"     if $nbam;

if($selfmodel){
	$mod_dir = "$somaticID_dir/self_trained_RDA";
	die( "Could not find directory $mod_dir\n" ) unless -d $mod_dir;
}

# create a directory for all the output files 
my $od        = "$pfx.somaticID";
system("mkdir $od") unless -d $od;

unless ($noanno) {
    $ltime = localtime();
    print "[$ltime] Checking/Downloading ANNOVAR database files\n";
    &annovar_database_file_check;
}

unless ( $novar or $noanno ) {
	$ltime = localtime();
    print "[$ltime] Calling Variants using SAMtools and VarScan programs\n";
    &make_mpileup_fifo;
    system(
        "$java -Djava.io.tmpdir=$tmpdir -Xmx$memory -jar \\
		$varscan mpileup2snp $pfx.mpileup.fifo --min-var-freq 0.01 \\
		--p-value 0.2 --output-vcf 1 > $od/$pfx.raw.snp.vcf"
    );
    system("rm $pfx.mpileup.fifo");
    &make_mpileup_fifo;
    system(
        "$java -Djava.io.tmpdir=$tmpdir -Xmx$memory -jar \\
		$varscan mpileup2indel $pfx.mpileup.fifo --min-var-freq 0.01 \\
		--p-value 0.2 --output-vcf 1 > $od/$pfx.raw.indel.vcf"
    );
    system("rm $pfx.mpileup.fifo");

    $ltime = localtime();
    print "[$ltime] Filtering variants based on the parameters:\n";
    print " " x 28, "--min-coverage 10\n";
    print " " x 28, "--min-reads2   4\n";
    print " " x 28,
      "--strand-filter 1 (Ignore variants with >90% support on one strand)\n";
    print " " x 28, "--min-var-freq 0.05\n";
    print " " x 28, "--p-value 0.01\n";
    system(
"$vcf_filter -i $od/$pfx.raw.snp.vcf --min-var-freq 0.05 -o $od/$pfx.flt.snp.vcf"
    );
    system(
"$vcf_filter -i $od/$pfx.raw.indel.vcf --min-var-freq 0.05 -o $od/$pfx.flt.indel.vcf"
    );

    if ( $adj_d > 0 ) {
        $ltime = localtime();
        print "[$ltime] Removing adjacent errors locating within $adj_d bp:\n";
        system(
            "mv $od/$pfx.flt.snp.vcf $od/$pfx.adjUnRM.snp.vcf;
			$rmAdj -d $adj_d $od/$pfx.adjUnRM.snp.vcf > $od/$pfx.flt.snp.vcf;
			mv $od/$pfx.flt.indel.vcf $od/$pfx.adjUnRM.indel.vcf;
			$rmAdj -d $adj_d $od/$pfx.adjUnRM.indel.vcf > $od/$pfx.flt.indel.vcf"
        );
    }
}

my ( %mean_baf, %sumbaf, %nbaf );
unless ($noanno) {

    die("Missing file: $od/$pfx.flt.snp.vcf\n")   unless -f "$od/$pfx.flt.snp.vcf";
    die("Missing file: $od/$pfx.flt.indel.vcf\n") unless -f "$od/$pfx.flt.indel.vcf";
    $ltime = localtime();
    print "[$ltime] Annotating the variants using ANNOVAR\n";
    print " " x 28, "ANNOVAR version: 2015-12-14\n";
    print " " x 28, "ANNOVAR protocol: $annovar_protocol\n";
    system(
        "$annovardir/table_annovar.pl $od/$pfx.flt.snp.vcf $annovardir/humandb/ \\
		-buildver $buildver -out $od/$pfx.flt.snp -remove \\
		-protocol $annovar_protocol \\
		-operation $annovar_operation -nastring . -vcfinput"
    );
	system("rm $od/$pfx.flt.snp.avinput");
    system(
        "$annovardir/table_annovar.pl $od/$pfx.flt.indel.vcf $annovardir/humandb/ \\
		-buildver $buildver -out $od/$pfx.flt.indel -remove \\
		-protocol $annovar_protocol \\
		-operation $annovar_operation -nastring . -vcfinput"
    );
	system("rm $od/$pfx.flt.indel.avinput");

    $ltime = localtime();
    print
"[$ltime] Calculating mean B-allele frequency (BAF) of the surrounding 40Mb regions\n";
    open IN, "$od/$pfx.flt.snp.$buildver\_multianno.vcf" or die($!);
    
	# whole chromosome BAF if insufficient local data
    my ( %sumbaf_whole_chr, %nbaf_whole_chr );
    while (<IN>) {
        chomp;
        next if /^#/;
        my @a = split;
        next unless $a[7] =~ /snp138=rs/;
        my $pfm;
        if ( $a[7] =~ /;PopFreqMax=([^;]+);/ ) {
            $pfm = $1;
        }
        next if $pfm eq ".";
        next if $pfm < 0.01;
        for ( my $j = 9 ; $j < $n_tum + 9 ; $j++ ) {
            next if $a[$j] =~ /^\.\/\./;
            my $k = join "\t", $a[0], $j;
            $sumbaf_whole_chr{$k} += &baf( $a[$j] );
            $nbaf_whole_chr{$k}++;
        }
        for (
            my $mb10 = int( $a[1] / 10**7 ) - 1 ;
            $mb10 <= int( $a[1] / 10**7 ) + 2 ;
            $mb10++
          )
        {
            $mb10 = $mb10 * 10**7;
            for ( my $j = 9 ; $j < $n_tum + 9 ; $j++ ) {
                next if $a[$j] =~ /^\.\/\./;

                # Throw low depth sites
                next if &total_depth( $a[$j] ) < 20;
                my $k = join "\t", $a[0], $mb10, $j;
                $sumbaf{$k} += &baf( $a[$j] );
                $nbaf{$k}++;
            }
            $mb10 += 10**7;
        }
    }
    close IN;
    foreach my $k ( keys %nbaf ) {
        my @kf = split /\t/, $k;
        my $k2 = join "\t", @kf[ 0, 2 ];
        $mean_baf{$k2} =
          sprintf( "%.3f", $sumbaf_whole_chr{$k2} / $nbaf_whole_chr{$k2} );

        # do not use local mean baf if there is limited probe SNP data
        $mean_baf{$k} = sprintf( "%.3f", $sumbaf{$k} / $nbaf{$k} )
          if $nbaf{$k} >= 8;
    }
    undef %sumbaf;
    undef %nbaf;
    undef %sumbaf_whole_chr;
    undef %nbaf_whole_chr;

    $ltime = localtime();
    print
      "[$ltime] Filtering out silent variants that do not change proteins\n";
    system(
"$non_silent_filter $od/$pfx.flt.snp.$buildver\_multianno.vcf > $od/$pfx.flt.snp.nonsilent.vcf"
    );
    system(
"$non_silent_filter $od/$pfx.flt.indel.$buildver\_multianno.vcf > $od/$pfx.flt.indel.nonsilent.vcf"
    );

    $ltime = localtime();
    print "[$ltime] Tabulating the VCF files\n";
    &tabulate( "$od/$pfx.flt.snp.nonsilent.vcf",   "$od/$pfx.snp.tab" );
    &tabulate( "$od/$pfx.flt.indel.nonsilent.vcf", "$od/$pfx.indel.tab" );
}
unless ($nopred) {
	$ltime = localtime();
	print "[$ltime] Classify somatic and non-somatic mutations\n";
	system("$rscript $rfc $mod_dir $od/$pfx.snp.tab $n_tum $od/$pfx");
	&som_ind;   # classify indels by simple cutoffs
	
	# Resume VCF files from the predicted results
	my @sfxs = qw(
	      somatic non-somatic CommonSNPs
	      LowQualityMutations UnknownMutationType
	    );
	foreach my $sf (@sfxs) {
	    my $infile = "$od/$pfx.$sf.txt";
		my $ori_vcf = "$od/$pfx.flt.snp.nonsilent.vcf";
		my $outfile = "$od/$pfx.$sf.snp.vcf";
		&re_VCF($sf, $infile, $ori_vcf, $outfile);
		system("rm $infile");
	}
	foreach my $sf (@sfxs[0,1]) {
	    my $infile = "$od/$pfx.$sf.indel.txt";
		my $ori_vcf = "$od/$pfx.flt.indel.nonsilent.vcf";
		my $outfile = "$od/$pfx.$sf.indel.vcf";
		&re_VCF($sf, $infile, $ori_vcf, $outfile);
		system("rm $infile");
	}
}


sub make_mpileup_fifo {
    system("$mkfifo $pfx.mpileup.fifo");
    if ($nbam) {
        system(
			"$samtools mpileup -B -q 37 -f $ref $tbams $nbam > $pfx.mpileup.fifo &"
        );
    }
    else {
        system(
            "$samtools mpileup -B -q 37 -f $ref $tbams > $pfx.mpileup.fifo &"
		);
    }
}

# check the existance of database files for ANNOVAR
sub annovar_database_file_check {
    my @anno_ptc = split /,/, $annovar_protocol;
    foreach my $ap (@anno_ptc) {
        next if -f "$annovardir/humandb/$buildver\_$ap.txt";
        
		# Try three times for downloading before quit the program
        for ( my $i = 0 ; $i < 3 ; $i++ ) {
            system(
				"$annovardir/annotate_variation.pl -buildver $buildver \\
					-downdb -webfrom annovar $ap $annovardir/humandb/"
            ) unless -f "$annovardir/humandb/$buildver\_$ap.txt";
        }
        die( "
Cannot download the ANNOVAR annoation file $annovardir/humandb/$buildver\_$ap.txt
Please try to download it manually and re-run the program.
Downloading guide: http://annovar.openbioinformatics.org/en/latest/user-guide/download/\n"
        ) unless -f "$annovardir/humandb/$buildver\_$ap.txt";
    }
}

# Transfer ANNOVAR-annoated VCF files into tab-deliminated files, 
# removing dbSNP sites with PopFreqMax >= 0.01
sub tabulate {
    my $out_tab_file = $_[1];
    open( my $fh, '>', $out_tab_file )
      or die("Could not open file $out_tab_file $!");
    my $header = join "\t",
      "SampleID", "VarID", "MutType", "PopFreqMax", "nci60", "clinvar", 
	  "cosmic70", "snp138", "SIFT", "Polyphen2_HDIV", "LRT", "MutationTaster", 
	  "MutationAssessor", "FATHMM", "PROVEAN", "VEST3", "CADD_phred", "DANN", 
	  "fathmm_MKL_coding", "MetaSVM", "MetaLR", "integrated_fitCons", "GERP",
      "phyloP7way_vertebrate", "phyloP20way_mammalian", "phastCons7way_vertebrate", 
	  "phastCons20way_mammalian", "SiPhy_29way_logOdds", "Interpro_domain", 
	  "dbscSNV_ADA", "dbscSNV_RF", "RefDP",    "VarDP",    "VAF",    "BAF_s40M",
      "T2_RefDP", "T2_VarDP", "T2_VAF", "T2_BAF_s40M", "N_RefDP",  "N_VarDP",  "N_VAF";
    print $fh "$header\n";
    open IN, $_[0] or die($!);
    while (<IN>) {
        chomp;
        next if /^#/;
        my @a = split;

        # ignore variants with more than one mutant alleles
        next if $a[4] =~ /,/;
        my $var_id = join "|", @a[ 0, 1, 3, 4 ];
        
		# For some features, replace "." with "0" to avoid being replaced as "NA".
        $a[7] =~ s/;snp138=\.;/;snp138=0;/;
        $a[7] =~ s/;Interpro_domain=\.;/;Interpro_domain=0;/;
        $a[7] =~ s/;Interpro_domain=\S[^;]+;/;Interpro_domain=1;/;
        $a[7] =~ s/;PopFreqMax=\.;/;PopFreqMax=0;/;
        $a[7] =~ s/;nci60=\.;/;nci60=0;/;
        $a[7] =~ s/;cosmic70=\.;/;cosmic70=0;/;
        $a[7] =~ s/;clinvar_20150629=\.;/;clinvar_20150629=0;/;

        # Add splicing to the missing exonic function field
        if ( $a[7] =~ /;Func.refGene=splicing;\S+;ExonicFunc.refGene=\.;AA/ ) {
            $a[7] =~
              s/;ExonicFunc.refGene=\.;AA/;ExonicFunc.refGene=splicing;AA/;
        }

        # Reverse MutationTaster score for "N" ("polymorphism")
        # and "P" ("polymorphism_automatic") (1)
        if ( $a[7] =~
            /;MutationTaster_score=(\S[^;]+);MutationTaster_pred=[NP];/ ) {
            $a[7] =~ s/;MutationTaster_score=/;MutationTaster_score=one_minus/;
        }
        my @infos = split /;/, $a[7];
        my @feas;
        push @feas, $var_id;

        # The index of @infos is very sensitive to the defination of ANNOVAR protocol
        foreach my $info (
            @infos[ 9,  12 .. 17, 19, 23, 25, 27, 29, 31,
            33, 35, 36, 37, 39, 41, 43, 45 .. 53 ]
          ) {
            if ( $info =~ /=(\S+)$/ ) {
                $info = $1;
            }
            $info = "NA" if $info eq ".";

            # only get clinical association types for ClinVar
            if ( $info =~ /^CLINSIG\\x3d([^\\]+)\\x3bCLNDBN/ ) {

                # 1: unknown, untested, non-pathogenic
                # 2: probable-non-pathogenic
                # 3: probable-pathogenic
                # 4: pathogenic, drug-response, histocompatibility, other
                my $match = $1;
                $info = 1 if $match =~ /unknown|untested|non-pathogenic/;
                $info = 2 if $match =~ /probable-non-pathogenic/;
                $info = 3 if $match =~ /probable-pathogenic/;
                $info = 4
                  if $match =~ /drug-response|histocompatibility|other/
                      || $match =~ /^pathogenic/
                      || $match =~ /[,\|]pathogenic/;
            }

            # Reverse MutationTaster score for "N" ("polymorphism")
            # and "P" ("polymorphism_automatic") (1)
            if ( $info =~ /^one_minus(\S+)$/ ) {
                $info = 1 - $1;
            }

            # Calculate total cosmic occurance as a single feature
            if ( $info =~ s/^ID\S+OCCURENCE\\x3d// ) {
                my $total_occ = 0;
                my @cosmic_occs = split /,/, $info;
                foreach my $cosmic_occ (@cosmic_occs) {
                    if ( $cosmic_occ =~ /^(\d+)\(\S+\)$/ ) {
                        $total_occ += $1;
                    }
                    else {
                        die("Incorrect cosmic70 OCCURENCE format: $info\n");
                    }
                }
                $info = $total_occ;
            }
            $info = 1 if $info =~ /^rs\d+$/;
            push @feas, $info;
        }

        # Get read depth info for each sample:
        my $dp_n = "NA\tNA\tNA";
        $dp_n = &depth_info_3( $a[$#a] ) if $nbam && $a[$#a] !~ /^\.\/\./;
        for ( my $i = 9 ; $i < $n_tum + 9 ; $i++ ) {
            next if $a[$i] =~ /^0\/0:/ || $a[$i] =~ /^\.\/\./;
            my $tum_sam_ID = $i - 8;
            $tum_sam_ID = "$pfx.$tum_sam_ID";
            my $dp_t1       = &depth_info_3( $a[$i] );
            my $baf_s40m    = &baf_sur40m( $a[0], $a[1], $i );
            my $dp_t2       = "NA\tNA\tNA";
            my $baf_s40m_t2 = "NA";
            my $n_d20 = 0;   # Number of other tumor samples with >= total depth
            my $max_t2_depth = 0;    # Highest depth of the other tumors

            # 1) If there is only one other tumor sample, 
			#    the read depth will be used here;
            # 2) If there are >1 other tumor sample, but <2 with DP>=20,
            #    then the one with highest Depth will be used;
            for ( my $j = 9 ; $j < $n_tum + 9 ; $j++ ) {
                next if $i == $j;
                next if $a[$j] =~ /^\.\/\./;
                my $tot_dp = &total_depth( $a[$j] );
                $n_d20++ if $tot_dp >= 20;
                last if $n_d20 > 1;
                if ( $tot_dp > $max_t2_depth ) {
                    $dp_t2 = &depth_info_3( $a[$j] );
                    $baf_s40m_t2 = &baf_sur40m( $a[0], $a[1], $j );
                }
                $max_t2_depth = $tot_dp;
            }

            # 3) If there are >=2 other tumor samples with DP>=20,
            #    then the one with largest difference in variant allele frequencies
            #    would be selected as the conterpart tumor sample
            if ( $n_d20 > 1 ) {
                my $max_vaf_diff = 0;
                for ( my $j = 9 ; $j < $n_tum + 9 ; $j++ ) {
                    next if $i == $j;
                    my $vaf_d = &vaf_diff( $a[$i], $a[$j] );
                    if ( $vaf_d > $max_vaf_diff ) {
                        $dp_t2 = &depth_info_3( $a[$j] );
                        $baf_s40m_t2 = &baf_sur40m( $a[0], $a[1], $j );
                    }
                    $max_vaf_diff = $vaf_d;
                }
            }
            print $fh join "\t", $tum_sam_ID, @feas, $dp_t1, $baf_s40m, $dp_t2,
              $baf_s40m_t2, "$dp_n\n";
        }
    }
    close $fh;
}

# Get absolute difference of variant allele frequencies of two samples from VCF depth column
sub vaf_diff {
    my @rds = split /:/, $_[0];
    my $vaf1 = $rds[5] / ( $rds[4] + $rds[5] );
    @rds = split /:/, $_[0];
    my $vaf2 = $rds[5] / ( $rds[4] + $rds[5] );
    my $r = $vaf1 - $vaf2;
    $r = -$r if $r < 0;
    return $r;
}

# calculate BAF of from VCF depth string
sub baf {
    my @rds = split /:/, $_[0];
    my $r = $rds[5] / ( $rds[4] + $rds[5] );
    $r = 1 - $r if $r > 0.5;
    return $r;
}

# obtain mean BAF of the surrounding 40Mb regions
sub baf_sur40m {
    my ( $bs_chr, $bs_pos, $bs_col ) = @_;
    $bs_pos = int( ( $bs_pos + 5 * 10**6 ) / 10**7 ) * 10**7;
    my $k = join "\t", $bs_chr, $bs_pos, $bs_col;
    my $r = "NA";
    if ( exists $mean_baf{$k} ) {
        $r = $mean_baf{$k};
    }
    else {
        my @kf = split /\t/, $k;
        my $k2 = join "\t", @kf[ 0, 2 ];
        if ( exists $mean_baf{$k2} ) {
            $r = $mean_baf{$k2};
        }
        else {
            print STDERR "WARNING: Undefined baf region: $k\n";
        }
    }
    return $r;
}

# Get total depth info from VCF depth column
sub total_depth {
    my @rds = split /:/, $_[0];
    my $r = $rds[4] + $rds[5];
}

# Get ref-depth, variant-depth and VAF info from VCF depth column
sub depth_info_3 {
    my @rds = split /:/, $_[0];
    my $rvaf = sprintf( "%.3f", $rds[5] / ( $rds[4] + $rds[5] ) );
    my $r = join "\t", @rds[ 4, 5 ], $rvaf;
}

# Resume VCF files from the predicted results
sub re_VCF {
	my ($sfx, $in_f, $in_vcf, $out_f) = @_;
    my %match;
	my $is_snp = 1;
    open( OUT, '>', $out_f )
          or die "cannot open file $out_f: $!";
    open IN, $in_f or die($!);
    while (<IN>) {
        chomp;
        next if $_ =~ /SampleID/;
        my @a = split;
		$is_snp = 0 if $#a == 41;
        $a[1] =~ s/"//g;
		if($is_snp){
			$match{ $a[1] } = $a[48];  # somatic probability
		}else{
			$match{ $a[1] } = 1;
		}
    }
    close IN;
    open IN, $in_vcf or die($!);
    while (<IN>) {
        chomp;
        if (/^##/) {
            print OUT "$_\n";
            next;
        }
        if (/^#/) {
            print OUT
"##INFO=<ID=somatic_score,Number=.,Type=String,Description=\"Confidence score that the mutation is a somatic\">\n"
	if $is_snp && $sfx =~ /somatic/;
            print OUT
"##FILTER=<ID=MultiVarAllele,Description=\"Mutation that have multiple variant alleles\">\n"
	if $is_snp;
            print OUT
"##FILTER=<ID=$sfx,Description=\"Only keep $sfx mutations in this VCF after mutation classification\">\n";
            print OUT "$_\n";
            next;
        }

        # ignore variants with more than one mutant alleles
        my @a = split;
        next if $a[4] =~ /,/;
        my $var_id = join "|", @a[ 0, 1, 3, 4 ];
        if ( exists $match{$var_id} ) {
			$a[7] .= ";somatic_score=$match{$var_id}" 
				if $is_snp && $sfx =~ /somatic/;
            print OUT join "\t", @a;
			print OUT "\n";
        }
    }
    close IN;
    close OUT;
}

# indel classifier by snp and pop freq
sub som_ind {
    open( SOM, '>', "$od/$pfx.somatic.indel.txt" );
    open( NSO, '>', "$od/$pfx.non-somatic.indel.txt" );
	open IN, "$od/$pfx.indel.tab" or die($!);
	while (<IN>){
		chomp;
		if (/^SampleID/){
			print SOM "$_\n";
			print NSO "$_\n";
			next;
		}
		my @a = split;

		# PopFreq==0 and non-SNP
		if($a[3] == 0 && $a[7] == 0) {
			print SOM "$_\n";
		}else{
			print NSO "$_\n";
		}
	}
	close IN;
	close SOM;
	close NSO;
}

sub Usage {
    die(
        qq/
USAGE: somaticID.pl [OPTIONS] -o prefix -ref genome.fa tumor1.bam [tumor2.bam ...]

NOTE: 1. Each BAM file should be sorted and represent alignment data from one sample. 
      2. All the samples should be from the same individual.

OPTIONS:
	--output|-o STRING         Prefix of the outputed files
	--reference|-ref *.fa      Positon of reference genome as a single FASTA file
	--normal *.bam             If the matched normal sample presents, regular somatic mutations calling will be performed using VarScan 2
	--rmAdjErr INT             Remove adjcent errors that locate within INT bp range [0]
	--novarcall                Let the program start from annotation, requires VCF files from VarScan 2
	--noanno                   Let the program start from somatic mutation identification, requires VCF files annotated by specific protocol
	--nopred                   Do not predict somatic mutations at this time (only generate TAB files for within-study model training)
	--selfmodel                Use self-trained models instead of pre-trained models
	--samtools|-s STRING       Path to SAMtools program [samtools]
	--varscan STRING           Path to VarScan.xxx.jar [~\/bin\/VarScan.v2.3.6.jar]
	--annovardir STRING        Path to the ANNOVAR program (tested version: 2015-12-14) directory [~\/bin\/annovar]
	--mkfifo|-mk STRING        Path to mkfifo [\/usr\/bin\/mkfifo]
	--java|-j STRING           Path to java [java]
	--rscript STRING           Path to Rscript [Rscript]  
	--memory|-me xxG\/M         Max memory for java [16G]
	--tmpjava|-t STRING        Position of a new tmp dir for java [\/tmp]
	--buildver|-b hgXX         Build version of human genome reference sequences [hg19]
	--help|-h                  Print this page
\n/
    );
}

