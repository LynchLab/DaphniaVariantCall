#! /bin/usr/perl -w

#ref_genome index path and file
$work_dir="/N/u/xw63/Carbonate/daphnia/genome_index";
$ref_genome="$work_dir/PA42.4.1";
$pbs_file="./Sequence_Coverage.pbs";

$SampleID="PA2013"; 
$DATA_DIR="/N/dc2/scratch/xw63/$SampleID/Bwa";
$tmp_DIR=$DATA_DIR."/tmp";
$MaxNumberofSamples=125;
$emailaddress='ouqd@hotmail.com';

$BWA="~/bwa-0.7.17/bwa";
$JAVA="java -XX:ParallelGCThreads=4 -Xmx32g -Xms32g -Djavaio.tmpdir=$tmp_DIR -jar";
$PICARD="$JAVA /N/soft/rhel6/picard/2.8.1/picard.jar";
$samtools="/N/soft/rhel6/samtools/1.3.1/bin/samtools";
$Trimmomatic="$JAVA ~/Trimmomatic-0.36/trimmomatic-0.36.jar";
$GATK="$JAVA /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar";
$bam="/N/soft/rhel6/bamUtil/1.0.13/bam";

# The paths to the software used in this pipeline
# You must first make sure you have all these software installed and they are all functional

#Now we find the mpileup files and produce a batch file for them

open OUT1, ">$pbs_file" or die "cannot open file: $!";
print OUT1 
"#!/bin/bash 
#PBS -N Sequence_Coverage
#PBS -k o
#PBS -l nodes=1:ppn=16,walltime=12:00:00
#PBS -l vmem=100gb
#PBS -M $emailaddress
#PBS -m abe
#PBS -j oe

# Updated on 05/28/2018

set +x 
module rm gcc
module load gcc
module load gsl  

module load samtools
set -x
cd $work_dir
set +x
echo ===============================================================
echo 0. Computing the average coverage of NGS sequencing
echo ===============================================================
set -x
date";

$n=0;
$n1=0;
while ($n<=$MaxNumberofSamples+1) {
	$n=$n+1;
	$nstr001= sprintf ("%03d", $n-1);
	$OUTPUT="$DATA_DIR/$SampleID-$nstr001";
	$INPUT="$DATA_DIR/$SampleID-$nstr001-RG_Sorted_dedup_realigned_Clipped.bam";
	
	if(-e "$OUTPUT.mpileup"){ 
		$n1=$n1+1;	
		#print ", Okay, a mpileup file is found! lets make a mapgd pro file:$OUTPUT.proview\n"; 
		print "$n1: $OUTPUT.mpileup\n";
		print OUT1 "\n
time $GATK -T DepthOfCoverage -R $ref_genome.fasta -o $OUTPUT -I $INPUT &\n";			
	}
}
print OUT1 
"\nwait\n
date

echo ===============================================================
echo =============Task completed.===================
echo ===============================================================
";			
if ($n1>0)
{
	print "\n$n1 mpileup files are found in $DATA_DIR.\n\n";
	print "
	============================================================
	Type the following command to produce the mapgd proview files: 

	  qsub $pbs_file
	============================================================\n\n\n";
}
else
{
	print "\n\nNo mpileup file is found in $DATA_DIR.\n\n\n";
}

