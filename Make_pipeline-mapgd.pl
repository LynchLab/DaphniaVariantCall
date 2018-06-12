#! /bin/usr/perl -w

#ref_genome index path and file
$workdir="/N/u/xw63/Carbonate/daphnia/DaphniaVariantCall";
$ref_genome="PA42.4.1";

# The adapter file: an example (Bioo_Adapters.fa) can be found in the same directory ($workdir)
#$Adapters="/PATH/TO/Adapters.fa";
$Adapters="/N/u/xw63/Carbonate/daphnia/Bioo_Adapters.fa";

$SampleID="PA2013"; 
$DATA_DIR="/N/dc2/scratch/xw63/$SampleID/Hisat";
$HeaderFile="$DATA_DIR/PA42.header";
$MaxNumberofSamples=125;
$emailaddress='ouqd@hotmail.com';

# The paths to the software used in this pipeline
# You must first make sure you have all these software installed and they are all functional

#Now we find the mpileup files and produce a batch file for them

open OUT1, ">./mapgd_proview.pbs" or die "cannot open file: $!";
print OUT1 
"#!/bin/bash 
#PBS -N mapgd-proview
#PBS -k o
#PBS -l nodes=1:ppn=16,walltime=6:00:00
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
cd $DATA_DIR
set +x
echo ===============================================================
echo 0. Make a header file
echo ===============================================================
set -x
samtools view -H $DATA_DIR/PA2013-001-RG_Sorted_dedup_realigned_Clipped.bam > $HeaderFile
set +x
echo ===============================================================
echo 1. Make a pro file of nucleotide-read quartets -counts of A, C, G, and T, from the mpileup files of the clones.
echo ===============================================================
set -x";

$n=0;
$n1=0;
while ($n<=$MaxNumberofSamples+1) {
	$n=$n+1;
	$nstr001= sprintf ("%03d", $n-1);
	$OUTPUT="$DATA_DIR/$SampleID-$nstr001";
	print "$nstr001:$OUTPUT.mpileup";
if(-e "$OUTPUT.mpileup"){ 
	print ", Okay, a mpileup file is found! lets make a mapgd pro file:$OUTPUT.pro.txt\n"; 
	$n1=$n1+1;	
	print OUT1 "\nmapgd proview -i $OUTPUT.mpileup -H $HeaderFile > $OUTPUT.pro.txt &\n";			
}
else
 { 
	print ", Ops, this file is not found! \n"; 
 } 
}
print OUT1 "\nwait\n";			

print "\n
============================================================
Type the following command to make the mapgd proview: 

  qsub ./mapgd_proview.pbs
============================================================\n\n";
