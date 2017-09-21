#!/bin/sh

output="pNOS.Javeed.Cus.SomVar.txt"

echo "CaseID" > $output
header="Chr	Pos	Ref	Alt	Gene	AAC	MutType	dbSNP	PopFreq	CADD_phred	Cosmic	ClinVar	MutGroup	SomaticProb	T.RD	T.VD	T.VAF"
echo $header >> $output
for name in `cat Javeed.pNOS.SampleID.r170518.txt SampleID.r170824.txt`
do
	echo $name >> $output 
	CusVarTable.somatic.pl $name.somaticID/$name | tail -n +2 >> $output
done

AddSampleIDColumn.pl $output > $output.1
msort -k b2,n3  $output.1 | HaemaTimesFromAnnovarCosmicColumn.pl - > $output
addNormalVarStat.pl -c 2 $output > $output.1
awk '$19<=2 && ($10=="." || $10<0.01)' $output.1 > $output
