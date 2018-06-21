## Pipelines for NGS genome-mapping and variant-calling 
### Lynch Lab, CME, Biodesign, ASU 
#### Curated and updated by Xiaolong Wang <ouqd@hotmail.com>
#### Initialized May 28, 2018


		======================================================
		   How to make and run the genome mapping pipelines  
		======================================================

This program is useful for genome mapping in large scale. If you have a large number (up to hundreds) of NGS sequencing reads to map onto a reference genome, following these steps for mapping the reads:
1.Save all Your raw reads in DATA_DIR in a sub dir ./fastq/, name you files in this way: 

	SampleID-001-R1.fastq
	SampleID-001-R2.fastq
	SampleID-002-R1.fastq
	SampleID-002-R2.fastq
			...
	SampleID-100-R1.fastq
	SampleID-100-R2.fastq

2.Indexing the reference genome using the following command:

	============================================================
		export ref_genome="PA42.4.0"
		
		bwa index $ref_genome.fasta 
		
		novoindex $ref_genome.ndx $ref_genome.fasta 
		
		hisat2-build $ref_genome.fasta $ref_genome

		samtools faidx $ref_genome.fasta
		
		PICARD  CreateSequenceDictionary R=$ref_genome.fasta O=$ref_genome.dict
	============================================================
	
3.Prepare the adapter files, and save it as:
	
	/PATH/To/Your/Adapters.fa

4.Edit MGMP.pl, change the settings according to your data:

	#Your reference genome: 
		$ref_genome="PA42.4.1"
		
	#Your adapter sequence file name: 
		$Adapters="/PATH/To/Your/Adapters.fa";
		
	#Your maximum number of file numbering: 
		$MaxNumberofSamples=150;
	
	#Your email address: 
		$emailaddress='xxx@xxx.xxx'
	
4.Make pipeline pbs files, to make and run the genome mapping/mapgd pipelines in large scale, simply run:

 ./Make_pipelines.sh

 Your screen will show a menu:

  =================================================================
  #                                                               #
  #     To make genome mapping pipelines, execute the             #
  # following commands, and follow the instructions they give:    #
  #                                                               #
  #     To use bwa:              				  #
  #                                                               #
		perl MGMP.pl bwa <DATA_DIR> <sampleID>
  #                                                               #
  #     To use novalign:          				  #
  #                                                               #
		perl MGMP.pl novoalign <DATA_DIR> <sampleID>
  #                                                               #
  #     To use hisat2:            				  #
  #                                                               #
		perl MGMP.pl hisat <DATA_DIR> <sampleID>
  #                                                               #
  =================================================================
Please copy and execute one of the above commands with appropriate args.
	
pipeline .pbs files will be produced for each pair of reads, and the pbs files are saved in the data directory in a sub directory pbs:

For bwa:

=================================================================

	DATA_DIR/pbs/bwa-SampleID-001.pbs
	DATA_DIR/pbs/bwa-SampleID-002.pbs
			...
	DATA_DIR/pbs/bwa-SampleID-100.pbs
	
  =================================================================

For novoalign:
	
  =================================================================
  
	DATA_DIR/pbs/novoalign-SampleID-001.pbs
	DATA_DIR/pbs/novoalign-SampleID-002.pbs
			...
	DATA_DIR/pbs/novoalign-SampleID-100.pbs

  =================================================================

For hisat:
	
  =================================================================
  
	DATA_DIR/pbs/hisat-SampleID-001.pbs
	DATA_DIR/pbs/hisat-SampleID-002.pbs
			...
	DATA_DIR/pbs/hisat-SampleID-100.pbs
	
  =================================================================
	

Run the  perl script MGMP.pl with appropriate args for the aligner you choose to make the pipelines and then, 

5. The perl script MGMP.pl also produced a shell script to submit these pbs jobs.
To submit the produced pipeline pbs jobs to the system for computing:
	
	For bwa:

	======================================
		chmod 755 qsub_all_pbs-bwa.sh
		./qsub_all_pbs-bwa.sh
	======================================

	For hisat:
	
	======================================
		chmod 755 qsub_all_pbs-hisat.sh
		./qsub_all_pbs-hisat.sh
	======================================

	For novoalign:

	======================================
		chmod 755 qsub_all_pbs-novoalign.sh
		./qsub_all_pbs-novoalign.sh
	======================================

==========End=======
