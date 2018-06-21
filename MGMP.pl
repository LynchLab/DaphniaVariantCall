#! /bin/usr/perl -w
#Save all your raw reads in the DATA_DIR in a sub dir named as SampleID/fastq/
#Name you files like: 
#SampleID-001-R1.fastq
#SampleID-001-R2.fastq
#SampleID-002-R1.fastq
#SampleID-002-R2.fastq
#......
#SampleID-100-R1.fastq
#SampleID-100-R2.fastq

#alignment tool: bwa, hisat, or novoalign

$aln=lc($ARGV[0]);

if ($aln eq "")
{
	print "\n\nYou have not selected a alignment tool. 
	The first input args is the alignment tool, 
	it must be:  bwa, hisat, or novoalign.\n\n"; 
	exit
}
print "\n\nYou have selected the alignment tool: $aln\n\n";

if ($aln ne "bwa"&&$aln ne "hisat"&&$aln ne "novoalign")
{
	print "\nThe 1st input (args) is invalid. 
	The 1st args is the alignment tool, 
	it must be: bwa, hisat, or novoalign.\n\n"; 
	exit
}

$DATA_DIR=$ARGV[1];

if(!(-e (glob($DATA_DIR))[0]))
{
	print "\nThe 2nd input (args) is the data directory. The data directory is not found:
				$DATA_DIR
	
	"; 
	exit
}

print "\n The data directory is:
				$DATA_DIR
	"; 

$SampleID=$ARGV[2]; 

if ($SampleID eq "")
{
	print "\n\nPlease input a population/Sample ID (The 3rd args).\n\n"; 
	exit
}

$tmp_DIR=$DATA_DIR."/tmp";
$MaxNumberofSamples=125;
$emailaddress='ouqd@hotmail.com';

#ref_genome index path and file
$work_dir="/N/u/xw63/Carbonate/daphnia/genome_index";
$ref_genome="PA42.4.1";

# The adapter file: an example (Bioo_Adapters.fa) can be found in the  directory:
#$Adapters="/PATH/TO/Adapters.fa";
$Adapters="/N/u/xw63/Carbonate/daphnia/Adapters/Bioo_Adapters.fa";



# The paths to the software used in this pipeline
# You must first make sure you have all these software installed and they are all functional
 
$BWA="~/bwa-0.7.17/bwa";
$JAVA="java -XX:ParallelGCThreads=4 -Xmx32g -Xms32g -Djavaio.tmpdir=$tmp_DIR -jar";
$PICARD="$JAVA /N/soft/rhel6/picard/2.8.1/picard.jar";
$samtools="/N/soft/rhel6/samtools/1.3.1/bin/samtools";
$Trimmomatic="$JAVA ~/Trimmomatic-0.36/trimmomatic-0.36.jar";
$GATK="$JAVA /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar";
$bam="/N/soft/rhel6/bamUtil/1.0.13/bam";
$hisat="/N/soft/rhel7/hisat2/2.1.0/hisat2";
$novoalign="/N/soft/rhel6/novoalign/novocraft/novoalign";
$novoindex="/N/soft/rhel6/novoalign/novocraft/novoindex";
$novoindex_ref_genome= "$ref_genome.ndx";

#Now we find the raw reads and produce pbs files for them

open OUT1, ">./qsub_all_pbs-$aln.sh" or die "cannot open file: $!";

$n=0;
$n1=0;
$ppn=16;

while ($n<=$MaxNumberofSamples+1) {
	$n=$n+1;
	$nstr001= sprintf ("%03d", $n-1);
	$walltime=(($n==1)?"01":(($aln eq "novoalign")?"48":"24"));
	$Sample=$DATA_DIR."/".$SampleID."-".$nstr001;
	$Sample_R1=$DATA_DIR."/fastq/".$SampleID."-".$nstr001."-R1";
	$Sample_R2=$DATA_DIR."/fastq/".$SampleID."-".$nstr001."-R2";
	$OUTPUT_DIR=$DATA_DIR."/".ucfirst($aln);
	$OUTPUT=$OUTPUT_DIR."/".$SampleID."-".$nstr001;
if(-e $Sample_R1.".fastq" && -e $Sample_R2.".fastq"){ 
	#print ", Okay, this pair-end reads fastq file is found! lets make a pbs file:"; 
	$n1=$n1+1;	
	$pbsfile=$DATA_DIR."/pbs/$aln-".$SampleID."-".$nstr001.".pbs";
	print $n1.": ";
	print $SampleID."-".$nstr001."-R1/R2.fastq  -->  ";
	print $pbsfile."\n";
	print OUT1 "\nqsub ".$pbsfile;	

	if ($aln eq "bwa")
	{	
#print "Alignment tool is bwa";

$Aln_comand="	
time $BWA mem -t $ppn -M -k 30 $ref_genome.fasta $Sample_R1-paired.fq $Sample_R2-paired.fq > $OUTPUT-paired.sam &
	
time $BWA mem -t $ppn -M -k 30 $ref_genome.fasta $Sample_R1-unpaired.fq > $OUTPUT-R1-unpaired.sam & 
	
time $BWA mem -t $ppn -M -k 30 $ref_genome.fasta $Sample_R2-unpaired.fq > $OUTPUT-R2-unpaired.sam &";
$Combine_comand="time $PICARD MergeSamFiles I=$OUTPUT-paired.sam I=$OUTPUT-R1-unpaired.sam I=$OUTPUT-R2-unpaired.sam O=$OUTPUT.sam"

	}	
	if ($aln eq "hisat")
	{
#	print "Alignment tool is hisat";
$Aln_comand="
time $hisat --no-spliced-alignment -p $ppn -q -x $ref_genome -1 $Sample_R1-paired.fq -2 $Sample_R2-paired.fq -S $OUTPUT-paired.sam &
	
time $hisat --no-spliced-alignment -p $ppn -q -x $ref_genome -U $Sample_R1-unpaired.fq -S $OUTPUT-R1-unpaired.sam &
	
time $hisat --no-spliced-alignment -p $ppn -q -x $ref_genome -U $Sample_R2-unpaired.fq -S $OUTPUT-R2-unpaired.sam &
";
$Combine_comand="time $PICARD MergeSamFiles I=$OUTPUT-paired.sam I=$OUTPUT-R1-unpaired.sam I=$OUTPUT-R2-unpaired.sam O=$OUTPUT.sam"
	}
	if ($aln eq "novoalign")
	{
#	print "Alignment tool is novoalign";
	
$Aln_comand="
time $novoalign -d $novoindex_ref_genome -r None -o Sam -f $Sample_R1-paired.0.fastq $Sample_R2-paired.0.fastq > $OUTPUT-paired.0.sam &
time $novoalign -d $novoindex_ref_genome -r None -o Sam -f $Sample_R1-paired.1.fastq $Sample_R2-paired.1.fastq > $OUTPUT-paired.1.sam &
time $novoalign -d $novoindex_ref_genome -r None -o Sam -f $Sample_R1-paired.2.fastq $Sample_R2-paired.2.fastq > $OUTPUT-paired.2.sam &
time $novoalign -d $novoindex_ref_genome -r None -o Sam -f $Sample_R1-paired.3.fastq $Sample_R2-paired.3.fastq > $OUTPUT-paired.3.sam &
time $novoalign -d $novoindex_ref_genome -r None -o Sam -f $Sample_R1-paired.4.fastq $Sample_R2-paired.4.fastq > $OUTPUT-paired.4.sam &
time $novoalign -d $novoindex_ref_genome -r None -o Sam -f $Sample_R1-paired.5.fastq $Sample_R2-paired.5.fastq > $OUTPUT-paired.5.sam &
time $novoalign -d $novoindex_ref_genome -r None -o Sam -f $Sample_R1-paired.6.fastq $Sample_R2-paired.6.fastq > $OUTPUT-paired.6.sam &
time $novoalign -d $novoindex_ref_genome -r None -o Sam -f $Sample_R1-paired.7.fastq $Sample_R2-paired.7.fastq > $OUTPUT-paired.7.sam &
time $novoalign -d $novoindex_ref_genome -r None -o Sam -f $Sample_R1-unpaired.fq > $OUTPUT-R1-unpaired.sam &
time $novoalign -d $novoindex_ref_genome -r None -o Sam -f $Sample_R2-unpaired.fq > $OUTPUT-R2-unpaired.sam &";

$Combine_comand="time $PICARD MergeSamFiles I=$OUTPUT-paired.0.sam I=$OUTPUT-paired.1.sam I=$OUTPUT-paired.2.sam I=$OUTPUT-paired.3.sam I=$OUTPUT-paired.4.sam I=$OUTPUT-paired.5.sam I=$OUTPUT-paired.6.sam I=$OUTPUT-paired.7.sam I=$OUTPUT-R1-unpaired.sam I=$OUTPUT-R2-unpaired.sam O=$OUTPUT.sam"

	}
open OUT, ">$pbsfile" or die "cannot open file: $!";
my $localtime = localtime();  
print OUT 
"#!/bin/bash	
#PBS -N $aln-$SampleID-$nstr001
#PBS -l nodes=1:ppn=$ppn
#PBS -l vmem=100gb
#PBS -l walltime=$walltime:00:00
#PBS -M $emailaddress
#PBS -m abe
#PBS -j oe

# This pipeline pbs is produced by the perl script:
# perl Make_pipelines-Genome-mapping.pl $ARGV[0] $ARGV[1] $ARGV[2]
# Date and time: $localtime

set +x
module load samtools
module load java
ulimit -s
set -x
cd $work_dir
mkdir $OUTPUT_DIR
mkdir $tmp_DIR
set +x
echo ===============================================================
echo 0. making index files, which should be done before submitting pbs
echo ===============================================================
echo samtools faidx $ref_genome.fasta
echo bwa index $ref_genome.fasta 
echo hisat2-build $ref_genome.fasta $ref_genome 
echo $novoindex $novoindex_ref_genome $ref_genome.fasta 
echo rm $ref_genome.dict
echo $PICARD  CreateSequenceDictionary R=$ref_genome.fasta O=$ref_genome.dict
echo These commands should be executed before submitting this pbs.
echo DO NOT excute these commands repeatedly in the pbs jobs, as it will cause problems when one job is using the index while another job is re-creating the index.
echo ===============================================================
echo 1. After preparing the FASTA file of adapter sequences, trim adapter sequences from sequence reads.
echo ===============================================================

module load java
set -x

time $Trimmomatic PE $Sample_R1.fastq $Sample_R2.fastq $Sample_R1-paired.fq $Sample_R1-unpaired.fq $Sample_R2-paired.fq $Sample_R2-unpaired.fq HEADCROP:3 ILLUMINACLIP:$Adapters:2:30:10:2 SLIDINGWINDOW:4:15 MINLEN:30
set +x
echo ===============================================================
echo 3. Mapping reads to the reference sequence and output bam file.
echo ===============================================================

set -x

$Aln_comand

wait 

set +x
echo ===============================================================
echo 4. Combine the SAM files using Picard.
echo ===============================================================

module load java
set -x

$Combine_comand

set +x
echo ===============================================================
echo 5. Convert the SAM file to the BAM file.
echo ===============================================================

set -x
	
time $samtools view -bS $OUTPUT.sam > $OUTPUT.bam
set +x
echo ===============================================================
echo 6. Sort the BAM file using Picard.
echo ===============================================================

module load java
set -x
	
time $PICARD SortSam INPUT=$OUTPUT.bam OUTPUT=$OUTPUT-Sorted.bam SORT_ORDER=coordinate
set +x

echo ===============================================================
echo 7. Add read groups to the sorted BAM file.
echo ===============================================================
module load java
set -x

time $PICARD AddOrReplaceReadGroups INPUT=$OUTPUT-Sorted.bam OUTPUT=$OUTPUT-RG_Sorted.bam RGID=Daphnia RGLB=bar RGPL=illumina RGSM=$Sample RGPU=6
set +x
echo ===============================================================
echo 8. Mark duplicate reads.
echo ===============================================================

module load java
set -x
time $PICARD MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=100 MAX_RECORDS_IN_RAM=5000000 VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true INPUT=$OUTPUT-RG_Sorted.bam OUTPUT=$OUTPUT-RG_Sorted_dedup.bam METRICS_FILE=$OUTPUT-metrics.txt
set +x
echo ===============================================================
echo 9. Index the BAM file using Picard.
echo ===============================================================

module load java
set -x
time $PICARD BuildBamIndex INPUT=$OUTPUT-RG_Sorted_dedup.bam
set +x
echo ===============================================================
echo 10. Define intervals to target for the local realignment.
echo ===============================================================

module load java
set -x
time $GATK -T RealignerTargetCreator -R $ref_genome.fasta -I $OUTPUT-RG_Sorted_dedup.bam -o $OUTPUT.intervals
set +x
echo ===============================================================
echo 11. Locally realign reads around indels.
echo ===============================================================

module load java
set -x

time $GATK -T IndelRealigner -R $ref_genome.fasta -I $OUTPUT-RG_Sorted_dedup.bam -targetIntervals $OUTPUT.intervals -o $OUTPUT-RG_Sorted_dedup_realigned.bam
set +x
echo ===============================================================
echo 12. Clip overlapping read pairs.
echo ===============================================================

set -x

time $bam clipOverlap --in $OUTPUT-RG_Sorted_dedup_realigned.bam --out $OUTPUT-RG_Sorted_dedup_realigned_Clipped.bam
set +x
echo ===============================================================
echo 13. Index the clipped BAM file using Samtools
echo ===============================================================

set -x

time $samtools index $OUTPUT-RG_Sorted_dedup_realigned_Clipped.bam
set +x
echo ===============================================================
echo 14. Make the mpileup file from the BAM file.
echo ===============================================================

set -x
	
time $samtools mpileup -f $ref_genome.fasta $OUTPUT-RG_Sorted_dedup_realigned_Clipped.bam > $OUTPUT.mpileup

set +x
echo ===============================================================
echo =============Task completed.===================
echo ===============================================================
	
"
}
else
	{ 
		#print ", Ops, this file is not found! \n"; 
	} 
}
if ($n1==0)
{
	print "No R1/R2.fq read file is found in $DATA_DIR.\n\n\n";
}
else
{
print "\n
============================================================
$n1 pbs files are produced and saved in: 
	$DATA_DIR/pbs/
  $aln-$SampleID-000.pbs, 
  $aln-$SampleID-001.pbs,
  ... ... 
  $aln-$SampleID-$n.pbs
============================================================
To submit all of the pbs jobs, type the following commands: 
   chmod 755 ./qsub_all_pbs-$aln.sh 
   ./qsub_all_pbs-$aln.sh 
============================================================
Before submitting these pbs files, reference genome index files
must first be made by using the following commands: 
 \$ samtools faidx $ref_genome.fasta 
 \$ bwa index $ref_genome.fasta 
 \$ hisat2-build $ref_genome.fasta $ref_genome 
 \$ $novoindex $novoindex_ref_genome $ref_genome.fasta 
 \$ rm $ref_genome.dict 
 \$ java -jar /N/soft/rhel6/picard/2.8.1/picard.jar  CreateSequenceDictionary R=PA42.4.1.fasta O=PA42.4.1.dict
============================================================
 PLEASE NOTE: DO NOT excute the indexing commands repeatedly in the jobs,
 as it will cause a problem when one job is using the index files 
 while another job is re-creating them.
============================================================
In these pbs files, $aln-$SampleID-000.pbs is useful for debuging these 
pipelines. The walltime of these pipelines is $walltime hours, 
while the walltime of $aln-$SampleID-000.pbs is only 1.00 hour, 
this will help to identify any problems quickly.
To debug, two small-sized pair-ended fastq read files, named as: 
		$SampleID-000-R1.fq
		$SampleID-000-R2.fq 
should be prepared and saved in the data directory:
		$DATA_DIR
Then, type the following commands: 
		qsub -q debug $DATA_DIR/pbs/$aln-$SampleID-000.pbs			
============================================================

";
}