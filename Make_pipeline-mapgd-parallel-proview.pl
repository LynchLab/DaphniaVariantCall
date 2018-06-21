#! /bin/usr/perl -w

#ref_genome index path and file
#$workdir="/N/u/xw63/Carbonate/daphnia/DaphniaVariantCall";
#$ref_genome="PA42.4.1";

# The adapter file: an example (Bioo_Adapters.fa) can be found in the same directory ($workdir)
#$Adapters="/PATH/TO/Adapters.fa";
#$Adapters="/N/u/xw63/Carbonate/daphnia/Bioo_Adapters.fa";

$SampleID="PA2013"; 
$DATA_DIR="/N/dc2/scratch/xw63/$SampleID/Bwa";
$MaxNumberofSamples=125;
$emailaddress='ouqd@hotmail.com';
$HeaderFile="$DATA_DIR/PA42.header";

# The paths to the software used in this pipeline
# You must first make sure you have all these software installed and they are all functional

#Now we find the mpileup files and produce a batch file for them

open OUT1, ">./mapgd-parallel-proview.pbs" or die "cannot open file: $!";
print OUT1 
"#!/bin/bash 
#PBS -N mapgd-parallel-proview
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
cd $DATA_DIR
set +x
echo ===============================================================
echo 0. Make a header file
echo ===============================================================
time samtools view -H PA2013-001-RG_Sorted_dedup_realigned_Clipped.bam > $HeaderFile
echo ===============================================================
echo 1. Make a pro file of nucleotide-read quartets -counts of A, C, G, and T, from each of the mpileup files of the clones.
echo ===============================================================
set -x
date";

$n=0;
$n1=0;
$n2=0;
while ($n<=$MaxNumberofSamples+1) {
	$n=$n+1;
	$nstr001= sprintf ("%03d", $n-1);
	$OUTPUT="$SampleID-$nstr001";
	
	if(-e "$DATA_DIR/$OUTPUT.mpileup"){ 
		$n1=$n1+1;	
		if(-e "$DATA_DIR/$OUTPUT.proview"){
			print "$DATA_DIR/$OUTPUT.proview already exist. if you want to reproduce it, please delete it first.\n"; 
		}
		else
		{		
			$n2=$n2+1;	
			#print ", Okay, a mpileup file is found! lets make a mapgd pro file:$OUTPUT.proview\n"; 
			print "$n1: $OUTPUT.mpileup\n";
			print OUT1 "\n			
time /N/dc2/projects/daphpops/Software/MAPGD-0.4.26/bin/mapgd proview -i $OUTPUT.mpileup -H $HeaderFile > $OUTPUT.proview &\n";	
		}			
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
print "\n$n1 mpileup file(s) are found in $DATA_DIR.\n\n";
}
else
{
	print "No mpileup file is found in $DATA_DIR.\n\n\n";
}

if ($n2>0)
{
	print "\n$n2 mpileup file(s) have no mapgd proview file(s).\n\n";
	print "
	============================================================
	Type the following command to produce the mapgd proview files: 

	  qsub ./mapgd-parallel-proview.pbs
	============================================================\n\n\n";
}
else
{
	print "All mpileup file(s) already have proview file(s) exist.\n\n\n";
	print "If you want to reproduce it, please delete it first.\n"; 
}

