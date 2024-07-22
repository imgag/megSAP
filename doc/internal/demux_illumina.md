# Demultiplexing of NovaSeq6000

1. Check basic run statistics with Illumina Sequence Analysis Viewer

    Error-rate < 1%  
    Cluster density  
    Check whether number of lanes is set correctly in NGSD  
    ...  

2. Set the run quality in the NGSD 

	Use comment field to describe problems.

3. Demultiplexing

	The demultiplexing splits the data of the lanes according to project/sample.  
	It reads BCL and writes FASTQ data.
	First we need to create a sample sheet:
	``` bash
	cd /mnt/storage3/raw_data/[instrument]/[run]/
	```
	```bash
	php /mnt/storage2/megSAP/pipeline/src/Tools/export_samplesheet.php -out SampleSheet_bcl2fastq.csv -run [ngsd_run_name]
	```	
	*Note*: The index (==barcode) should be given even if only one sample is on a lane (to avoid contamination).  
	*Note*: If index 2 sequences are reverse-complement (e.g in NovaSeq6000 or MiSeq dual-indexing) of the sequences given in the sample sheet, re-export using the '-mid2_no_rc' flag.  
	*Note*: Restrict the MID lenghs using the parameters '-mid1_len' and '-mid2_len'.      
	*Note*: If there are some samples with only MID_1 and no MID_2. Start demultiplexing without them and check in the unclaimed MIDs what sequence is their MID_2. Stop the demux and enter this found MID_2 into the sample sheet to demultiplex all samples together.
	
	Then we execute the actual demultiplexing based on the sample sheet:
	```bash
	ionice -c 3 /mnt/storage1/share/opt/bcl2fastq2_2.19.1/bcl2fastq -o Unaligned -p 10 --sample-sheet SampleSheet_bcl2fastq.csv --barcode-mismatches 1 >demux.log 2>&1
	```
	```bash
	tail demux.log
	```
 	The demultiplexing of a run with samples that contain UMI sequences (e.g. RNA samples) needs a modification of the given base-mask as well as additional arguments to make sure  the UMI reads are not filtered as to short:
	```bash
  	ionice -c 3 /mnt/storage1/share/opt/bcl2fastq2_2.19.1/bcl2fastq -o Unaligned -p 10 --sample-sheet SampleSheet_bcl2fastq.csv --barcode-mismatches 1 --use-bases-mask Y105,I8Y11,I8,Y105 --mask-short-adapter-reads 0 --minimum-trimmed-read-length=8  > demux.log 2>&1
   	``` 
	```bash
	tail demux.log
	```
 
	*Note*: The tail command must return "Processing completed with 0 errors and 0/1 warnings."  
	*Note*: If the MID lengths of the samples differ or you want to demultiplex different lengths etc. you have to provide the correct recipe using `--use-bases-mask`   
	*Note*: If you have a low barcode distance (Warning of `export_samplesheet.php`) set the parameter `--barcode-mismatches` to `0`   
	*Note*: The parameters '--tiles s_[n]' can be used to demultiplex one/several lanes only.  
	*Note*: To keep index cycles/short reads specified in the basemask argument, use '--minimum-trimmed-read-length=8 --mask-short-adapter-reads=8' to prevent obtaining NNNNNNNN-reads only.  
	*Note*: To use index as a read, write Y instead of I in the base mask argument. Those reads will be called as separate fastq files.  
	*Note*: If a sample has no second index, try AGATCTCGGT. For unknown bases 9/10 of first index read try AT.

	
4. Check that for all samples there is data and that no sample was missing in the sample sheet:
	```bash
	du -sh Unaligned/*/* Unaligned/Undetermined* | egrep -v "Reports|Stats"
	```

5. Check for additional MIDs in run data:
	```bash
	FastqMidParser -in Unaligned/Undetermined_S0_L00?_R?_001.fastq.gz
	```
6. Copy the FASTQ files to the project folders and queue the analysis using:
	```bash
	php /mnt/storage2/megSAP/pipeline/src/NGS/copy_sample.php
	```
	```bash
	ionice -c 3 make all
	```

7. Backup run using backup tool:
	```bash
	sudo -u archive-gs php /mnt/storage2/megSAP/pipeline/src/Tools/backup_queue.php -mode run -in [run] -email [email]
	```

8. Delete the run raw data (when all samples are analyzed with passed QC):
	```bash
	rm -rf [run]
	```

# Demultiplexing of NovaSeqXPlus
1. Make sure the DRAGEN analysis on the NovaSeq X is completed and all files are copied to the storage.   
	A `CopyComplete.txt`file is created in the run folder when the copy is done.    
	A `CopyComplete.txt`file is created in `Analysis/[1-9]` when the DRAGEN analysis is done.   

2. Check basic run statistics with Illumina Sequence Analysis Viewer

    Error-rate < 1%  
    Cluster density  
    Check whether number of lanes is set correctly in NGSD  
    ...  

3. Set the run quality in the NGSD 

	Use comment field to describe problems.

4. Demultiplexing

	Demultiplexing is done on the instrument already. We only have to check the output for additional MIDs:
	
	```bash
	python /mnt/storage2/megSAP/pipeline/src/Tools/extract_demux_stats.py
	```
	Check the output files `DemuxStats.html` and `DemuxStats.png` for additional MIDs.

5. Copy ORA/BAM to the project folders and analyze them using:
	```bash
	php /mnt/storage2/megSAP/pipeline/src/NGS/copy_sample.php
	```
	```bash
	ionice -c 3 make all
	```

6. Backup run using backup tool:
	```bash
	sudo -u archive-gs php /mnt/storage2/megSAP/pipeline/src/Tools/backup_queue.php -mode run -in [run] -email [email]
	```

7. Delete the run raw data (when all samples are analyzed and passed QC):
	```bash
	rm -rf [run]
	```

## Troubleshooting: manual copy of runs
1. ssh into NovaSeq X
2. copy run (meta) data:
```bash
rsync -v -r --progress /usr/local/illumina/mnt/runs/[run] /mnt/storage3/raw_data/LUNA/ # data      (no '/' at the end of [run])
rsync -v -r --progress /usr/local/illumina/runs/[run] /mnt/storage3/raw_data/LUNA/     # meta data (no '/' at the end of [run])
```  

# APPENDIX: How to install bcl2fastq2

1. Download and unzip the RPM file.

2. Install using the command:
	```bash
	sudo alien -i bcl2fastq2-v2.18.0.12-Linux-x86_64.rpm
	```

3. Copy the executable to the opt folder:
	```bash
	mkdir /mnt/storage1/share/opt/bcl2fastq2_1.18/
	cp /usr/local/bin/bcl2fastq /mnt/storage1/share/opt/bcl2fastq2_1.18/
	```
