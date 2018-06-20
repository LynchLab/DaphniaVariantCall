		==============================================================
		#	How to make and run the parallel mapgd pipeline  #
		==============================================================
		
1. After mapping all the reads to the reference genome, you have got a number of .mpileup files: 

	========================================================
		SampleID-001.mpileup
		SampleID-002.mpileup
				...
		SampleID-096.mpileup
	========================================================

2. Make proview files: 
	========================================================
	
        perl Make_pipeline-mapgd-proview.pl
        qsub ./mapgd-parallel-proview.pbs
		
	========================================================
Note: 

In the original pipeline, this step is implemented by the following command:
	=========================================================================
 
		mapgd proview -i *.mpileup -H $HeaderFile > output.pro.txt 
		
 	=========================================================================

However, this is very slow because it is not fully parallelized. This step  will takes up to 50-100 hours for a population with 96 clones. To reduce the computation time, this step is replaced with the pipeline, mapgd-parallel-proview.pbs (produced by Make_pipeline-mapgd-proview.pl). 

In mapgd-parallel-proview.pbs, the proview files are generated for each of the 96 clones:
 
	=========================================================================
	mapgd proview -i $SampleID-001.mpileup -H $HeaderFile > $SampleID-001.proview &
	mapgd proview -i $SampleID-002.mpileup -H $HeaderFile > $SampleID-002.proview &
	mapgd proview -i $SampleID-003.mpileup -H $HeaderFile > $SampleID-003.proview &
	
			... ...
			
	mapgd proview -i $SampleID-096.mpileup -H $HeaderFile > $SampleID-096.proview &
	
	=========================================================================

In this way, mapgd proview file is produced for each of the 96 clones.  Because all processes can be run simutaneously in independent threads, the computation time is greatly deduced:
Then, these mapgd proview files can be combined by the following command using a homemade java program (CombineProview.java):
	=============================================================
	
		java -cp ./ CombineProview <DATA_DIR> <output>
		
	=============================================================
	
3. Submit the parallel mapgd pipeline:
	===============================================

		qsub mapgd-parallel.pbs
		
	===============================================
This pipeline will combine all mapgd proview files into one using the abovemendtioned java program (CombineProview.java), and do the rest of the mapgd pipeline same as the original pipeline.

That's all, thank you!
