#! /bin/usr/perl -w

#ref_genome index path and file
$workdir="/N/u/xw63/Carbonate/daphnia/genome_index";
$ref_genome="PA42.4.1";

# The adapter file: an example (Bioo_Adapters.fa) can be found in the same directory ($workdir)
#$Adapters="/PATH/TO/Adapters.fa";
$Adapters="/N/u/xw63/Carbonate/daphnia/Adapters/Bioo_Adapters.fa";

#Save all your raw reads in the DATA_DIR in a sub dir named as SampleID/fastq/
#Name you files like: 
#SampleID-001-R1.fastq
#SampleID-001-R2.fastq
#SampleID-002-R1.fastq
#SampleID-002-R2.fastq
#......
#SampleID-100-R1.fastq
#SampleID-100-R2.fastq

$SampleID="PA2013"; 
$DATA_DIR="/N/dc2/scratch/xw63/".$SampleID;
$tmp_DIR=$DATA_DIR."/tmp";
$MaxNumberofSamples=125;
$emailaddress='ouqd@hotmail.com';

# The paths to the software used in this pipeline
# You must first make sure you have all these software installed and they are all functional
 
$novoalign="/N/soft/rhel6/novoalign/novocraft/novoalign";
$novoindex="/N/soft/rhel6/novoalign/novocraft/novoindex";
$novoindex_ref_genome= "$ref_genome.ndx";
$JAVA="java -XX:ParallelGCThreads=4 -Xmx32g -Xms32g -Djavaio.tmpdir=$tmp_DIR -jar";
$PICARD="$JAVA /N/soft/rhel6/picard/2.8.1/picard.jar";
$samtools="/N/soft/rhel6/samtools/1.3.1/bin/samtools";
$Trimmomatic="$JAVA ~/Trimmomatic-0.36/trimmomatic-0.36.jar";
$GATK="$JAVA /N/soft/rhel6/gatk/3.4-0/GenomeAnalysisTK.jar";
$bam="/N/soft/rhel6/bamUtil/1.0.13/bam";

#Now we find the raw reads and produce pbs files for them

open OUT1, ">./qsub_all_pbs-novoalign.sh" or die "cannot open file: $!";

$n=0;
$n1=0;
while ($n<=$MaxNumberofSamples+1) {
	$n=$n+1;
	$nstr001= sprintf ("%03d", $n-1);
		$Sample=$DATA_DIR."/".$SampleID."-".$nstr001;
		$Sample_R1=$DATA_DIR."/fastq/".$SampleID."-".$nstr001."-R1";
		$Sample_R2=$DATA_DIR."/fastq/".$SampleID."-".$nstr001."-R2";
		$OUTPUT_DIR=$DATA_DIR."/Novoalign";
		$OUTPUT=$OUTPUT_DIR."/".$SampleID."-".$nstr001;
if(-e $Sample_R1.".fastq" && -e $Sample_R2.".fastq"){ 
	#print ", Okay, this pair-end reads fastq file is found! lets make a pbs file:"; 
	$n1=$n1+1;	

	$pbsfile=$DATA_DIR."/pbs/novoalign-".$SampleID."-".$nstr001.".pbs";
	print $n1.": ";
	print $SampleID."-".$nstr001."-R1/R2.fastq  -->  ";
	print $pbsfile."\n";
	print OUT1 "\nqsub ".$pbsfile;			

open OUT, ">$pbsfile" or die "cannot open file: $!";
print OUT 
"#!/bin/bash	
#PBS -N Nov-$SampleID-$nstr001
#PBS -l nodes=1:ppn=8
#PBS -l vmem=100gb
#PBS -l walltime=24:00:00
#PBS -M $emailaddress
#PBS -m abe
#PBS -j oe

set +x
module load samtools
module load java
ulimit -s
set -x
cd $workdir
mkdir $OUTPUT_DIR
mkdir $tmp_DIR
set +x
echo ===============================================================
echo 0. making index files using the following commands
echo ===============================================================
echo samtools faidx $ref_genome.fasta
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
echo # 2. Split the FASTQ files into eight pieces.
echo ===============================================================
set -x	
ngsutilsj fastq-split -n 8 $Sample_R1-paired.fq &
ngsutilsj fastq-split -n 8 $Sample_R2-paired.fq &
wait
set +x
echo ===============================================================
echo 3. Mapping reads to the reference sequence and output bam file.
echo ===============================================================
set -x
time $novoalign -d $novoindex_ref_genome -r None -o Sam -f $Sample_R1-paired.0.fastq $Sample_R2-paired.0.fastq > $OUTPUT-paired.0.sam &
time $novoalign -d $novoindex_ref_genome -r None -o Sam -f $Sample_R1-paired.1.fastq $Sample_R2-paired.1.fastq > $OUTPUT-paired.1.sam &
time $novoalign -d $novoindex_ref_genome -r None -o Sam -f $Sample_R1-paired.2.fastq $Sample_R2-paired.2.fastq > $OUTPUT-paired.2.sam &
time $novoalign -d $novoindex_ref_genome -r None -o Sam -f $Sample_R1-paired.3.fastq $Sample_R2-paired.3.fastq > $OUTPUT-paired.3.sam &
time $novoalign -d $novoindex_ref_genome -r None -o Sam -f $Sample_R1-paired.4.fastq $Sample_R2-paired.4.fastq > $OUTPUT-paired.4.sam &
time $novoalign -d $novoindex_ref_genome -r None -o Sam -f $Sample_R1-paired.5.fastq $Sample_R2-paired.5.fastq > $OUTPUT-paired.5.sam &
time $novoalign -d $novoindex_ref_genome -r None -o Sam -f $Sample_R1-paired.6.fastq $Sample_R2-paired.6.fastq > $OUTPUT-paired.6.sam &
time $novoalign -d $novoindex_ref_genome -r None -o Sam -f $Sample_R1-paired.7.fastq $Sample_R2-paired.7.fastq > $OUTPUT-paired.7.sam &
time $novoalign -d $novoindex_ref_genome -r None -o Sam -f $Sample_R1-unpaired.fq > $OUTPUT-R1-unpaired.sam &
time $novoalign -d $novoindex_ref_genome -r None -o Sam -f $Sample_R2-unpaired.fq > $OUTPUT-R2-unpaired.sam &
wait 

set +x
echo ===============================================================
echo 4. Combine the SAM files using Picard.
echo ===============================================================
module load java
set -x
time $PICARD MergeSamFiles I=$OUTPUT-paired.0.sam I=$OUTPUT-paired.1.sam I=$OUTPUT-paired.2.sam I=$OUTPUT-paired.3.sam I=$OUTPUT-paired.4.sam I=$OUTPUT-paired.5.sam I=$OUTPUT-paired.6.sam I=$OUTPUT-paired.7.sam I=$OUTPUT-R1-unpaired.sam I=$OUTPUT-R2-unpaired.sam O=$OUTPUT.sam

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
$n1 pbs files are produced and saved in: $DATA_DIR/pbs/
  novoalign-$SampleID-000.pbs, 
  novoalign-$SampleID-001.pbs,
  ... ... 
  novoalign-$SampleID-$n.pbs
============================================================ 
In these pbs files, bwa-$SampleID-000.pbs is useful for debuging this 
pipeline, two small-sized pair-ended fastq read files, named as: 
		$SampleID-000-R1.fq
		$SampleID-000-R2.fq 
should be prepared and copied to the data directory:
		$DATA_DIR
Then, type the following commands: 
		qsub $DATA_DIR/pbs/novoalign-$SampleID-000.pbs
This will help to identify any problems quickly.
============================================================
To submit all of the pbs jobs, type the following commands: 
   chmod 755 ./qsub_all_pbs-novoalign.sh 
   ./qsub_all_pbs-novoalign.sh 
============================================================
Before submitting these pbs files, one must first make\n reference genome index files by using the following commands: 
  samtools faidx $ref_genome.fasta 
  $novoindex $novoindex_ref_genome $ref_genome.fasta 
  rm $ref_genome.dict 
  java -jar /N/soft/rhel6/picard/2.8.1/picard.jar  CreateSequenceDictionary R=PA42.4.1.fasta O=PA42.4.1.dict
  
DO NOT excute these commands repeatedly in the pbs jobs, as it will cause problems when one job is using the index files while another job is re-creating the index.
============================================================\n\n";
}