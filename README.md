# somaticID

```somaticID``` is a computational pipeline to identify somatic mutations when there is no matched normal data. It calls variants mutations from BAM files with [VarScan 2](http://varscan.sourceforge.net/), annotates the variants with [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/) and identify non-synonymous somatic variants based on random forest models.

```somaticID``` works with two conditions. To identify somatic mutations from samples from a single case, or multiple cases without matched normal data, a pre-trained model can be used; for datasets from multiple cases with similar conditions (tumor content, sequencing method, etc), and some of them have available matched normal data, new models can be trained with the cases with normal data, and then be used to identify somatic mutations from other cases. 

For manual, please use
```Bash
	./somaticID.pl -h
	./somaticID.train.pl -h
```
### Condition 1: Identify somatic mutations based on pre-trained models in single cases

##### Input:
```Bash
  Case1.Tumor1.bam
  Case1.Tumor2.bam
  ...
```

```Case1``` could be replaced with any other names representing a single patient, and ```Tumor1```, ```Tumor2``` ... could be any names for tumor samples from the same case. Please note the PCR duplication should be already marked or removed from the input BAM files. For details, please refer to [Picard’s MarkDuplicates](http://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates) or [SAMtools’ rmdup](http://www.htslib.org/doc/samtools.html).

##### Command:
```Bash
  ./somaticID.pl [OPTIONS] -o Case1 -ref genome.fa Case1.Tumor1.bam Case1.Tumor2.bam [...]
```
All the output files will be put in a subdirectory ```Case1.somaticID```, including: 
* The non-silent somatic variations:
```Bash
  Case1.somaticID/Case1.somatic.indel.vcf
  Case1.somaticID/Case1.somatic.snp.vcf
  ...
```
* The non-silent, non-somatic variations: 
```Bash
  Case1.somaticID/Case1.non-somatic.indel.vcf
  Case1.somaticID/Case1.non-somatic.snp.vcf
  Case1.somaticID/Case1.CommonSNPs.snp.vcf
  ...
```
* SNPs with low-quality or unknown mutation types : 
```Bash
  Case1.somaticID/Case1.LowQualityMutations.snp.vcf
  Case1.somaticID/Case1.UnknownMutationType.snp.vcf
  ...
```
* Other intermediate files that may be useful: 
```Bash
  Case1.somaticID/Case1.adjUnRM.indel.vcf
  Case1.somaticID/Case1.adjUnRM.snp.vcf
  Case1.somaticID/Case1.flt.indel.hg19_multianno.txt
  Case1.somaticID/Case1.flt.indel.hg19_multianno.vcf
  Case1.somaticID/Case1.flt.indel.nonsilent.vcf
  Case1.somaticID/Case1.flt.indel.vcf
  Case1.somaticID/Case1.flt.snp.hg19_multianno.txt
  Case1.somaticID/Case1.flt.snp.hg19_multianno.vcf
  Case1.somaticID/Case1.flt.snp.nonsilent.vcf
  Case1.somaticID/Case1.flt.snp.vcf
  Case1.somaticID/Case1.indel.tab
  Case1.somaticID/Case1.raw.indel.vcf
  Case1.somaticID/Case1.raw.snp.vcf
  Case1.somaticID/Case1.snp.tab
  ...
```
### Condition 2: Train models using new data and identify somatic mutations
##### Input:
```Bash
  Case1.Tumor1.bam
  Case1.Tumor2.bam
  Case1.Normal.bam
  Case2.Tumor1.bam
  Case2.Tumor2.bam
  Case2.Normal.bam
  Case3.Tumor1.bam
  Case3.Tumor2.bam
  ...
```
##### A three-step precedure:
1) For each case, identify mutations from BAM files without somatic mutation prediction. Option ```--nopred``` is used to prevent calling somatic mutations with pre-trained models. Note that for samples with matched normal samples, option ```--normal``` was used.

```Bash
  for case_id in "Case1" "Case2"
  do
    ./somaticID.pl [OPTIONS] -o $case_id -ref genome.fa \
      --nopred --normal $case_id.Normal.bam \
      $case_id.Tumor1.bam $case_id.Tumor2.bam
  done
  for case_id in "Case3"
  do
    ./somaticID.pl [OPTIONS] -o $case_id -ref genome.fa \
      --nopred $case_id.Tumor1.bam $case_id.Tumor2.bam
  done
```

2) Train new models based on cases with match normal data

```Bash
	./somaticID.train.pl --tabfiles Case1.tab,Case2.tab
```

3) Identify the somatic mutations based on self-trained models with option ```--selfmodel```

```Bash
	./somaticID.pl [OPTIONS] -o Case3 -ref genome.fa --noanno \
		--selfmodel Case3.Tumor1.bam Case3.Tumor2.bam
```
The output files for each case will be put in a seperate subdirectory, with similar name suffixes as noted above.

### Evaluation Mode
If ```somaticID``` is run with matched normal data, an evaluation process will performed and a CSV file will be generated to show the performance of somatic identification in given sample.
##### Input:
```Bash
  Case1.Tumor1.bam
  Case1.Tumor2.bam
  Case1.Normal.bam
  ...
```
##### Command:
```Bash
  ./somaticID.pl [OPTIONS] -o Case1 -ref genome.fa --normal Case1.Normal.bam \
  	Case1.Tumor1.bam Case1.Tumor2.bam [...]
```
The resulted CSV files will be ```Case1.somaticID/Case1.evaluation.txt```.
