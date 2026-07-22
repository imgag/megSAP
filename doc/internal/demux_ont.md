# ONT PromethION Data Transfer and Post-Run Process

1. Check run quality using MinKNOW UI software

    Basecalled bases >90Gb (good), 75 to 90Gb (medium), <75GB (bad) (per sample)   
    Approximate N50 >10kb (human WGS)   

2. Set the run quality in the NGSD 

	Use comment field to describe problems.

3. Data transfer

    Raw data is stored on the PromethION in `/data/[MINERVA|VULCAN]/[PROJECT|BATCH NAME]`.
   
    Connect to the device via ssh:
    ```bash
    ssh minerva
    ```
    Check if copy is completed:
   ```bash
   systemctl --user status rsync-seq-data.service
   ```
   Example output for completed copy:
   ```bash
   #systemctl --user status rsync-seq-data.service
	○ rsync-seq-data.service - rsync data from /data to storage3 buffer
	     Loaded: loaded (/home/prom/.config/systemd/user/rsync-seq-data.service; static)
	     Active: inactive (dead) since Sun 2026-06-07 09:49:50 CEST; 51min ago
	   Duration: 12.299s
	TriggeredBy: ● rsync-seq-data.timer
	    Process: 825938 ExecStart=rsync -aLi --no-owner --no-group --no-perms --omit-dir-times --exclude=fastq_fail/ --exclude=fastq_pass/ /data/VULCAN/ /mnt/storage3/raw_data/VULCAN/ (code=exited, status=0/SUCCESS)
	   Main PID: 825938 (code=exited, status=0/SUCCESS)
	        CPU: 141ms
	
	Jun 07 09:49:37 PCA100348 systemd[25641]: Started rsync-seq-data.service - rsync data from /data to storage3 buffer.
   ```

    
    
5. Copy raw data to buffer (diagnostic data):
	```bash
	cd /mnt/storage3/raw_data/[MINERVA/VULCAN]/
 	# diagnostic samples
 	ionice -c 3 rsync -ah --info=stats2,progress2 --omit-dir-times [BATCH_NAME] /mnt/storage3b/promethion_rawdata_buffer/[Minerva|Vulcan]
 	# research
 	ionice -c 3 rsync -ah --info=stats2,progress2 --omit-dir-times [BATCH_NAME] /mnt/storage3/promethion_rawdata_buffer/[YYYYmmdd]
 	 ```

6. Copy runs and start optional basecalling and analysis

    Using `copy_ont_data.php`, the flow cell data is copied to the
    default sample folder location, based on the information entered in
    GSvar.

    ```bash
    ulimit -n 10000
    php /mnt/storage2/megSAP/pipeline/src/IMGAG/copy_ont_run.php \
        -run_dir [run directory] \
        -queue_sample \
        [-queue_basecalling \
        -basecall_model [sup/hac]] \
        [-email [email]]

    # example for diagnostic:
    [cd /mnt/storage3/raw_data/[MINERVA/VULCAN]/[BATCH]]
    php /mnt/storage2/megSAP/pipeline/src/IMGAG/copy_ont_run.php [-email email@test.com] -basecall_model sup -queue_basecalling -queue_sample -run_dir DNA0000000_00000
    
    ```

    The result file is `[processed sample].mod.unmapped.bam` in the corresponding project/sample directory and the sample analysis is automatically queued. If `-queue_basecalling` is set, the basecalled BAM is also stored in the run folder.

    The script aborts if the flow cell directory contains a `pod5_skip`
    subdirectory, indicating a partially basecalled flow cell (see
    below) and no basecalling is queued.

7. Backup raw run data using backup script (**after basecalling is done!!!**):   
	By default no POD5 data is included in the backup, only the BAM files! 
	```bash
	sudo -u archive-gs php /mnt/storage2/megSAP/pipeline/src/IMGAG/backup_queue.php -mode run -in [run] -email [email]
	```

9. Delete the run raw data (when all samples are analyzed with passed QC, raw data is copied to buffer and backup is done):
	- delete from PromethION
 	- delete from /mnt/storage3/raw_data
  	- after 3 months: delete from buffer
	```bash
	rm -rf [run]
	```

## Basecalling

In case of software issues, the device may store non-basecalled data in
a `pod5_skip` subdirectory. This data can be basecalled on a GPU server
(e.g. SRV010) using `dorado`, and the basecalled data is stored in an
extra BAM file in `bam_pass` directory.

```bash
/mnt/storage2/megSAP/tools/dorado-0.6.0-linux-x64/bin/dorado \
    basecaller \
    hac,5mCG_5hmCG \
    <run directory>/<flow cell directory>/pod5_skip \
    --min-qscore 9 \
    > <run directory>/<flow cell directory>/bam_pass/pod5_skip_basecalled.bam
```

Rename the `pod5_skip` directory to continue with the `copy_ont_run.php`
script.

```bash
cd <run directory>/<flow cell directory>
mv pod5_skip pod5_skip_done
```

## Manual data copy:
Create a directory with the run start date:
```bash
mkdir /mnt/storage3/promethion_rawdata_buffer/<YYYYmmdd>
```

Copy the data using `rsync`:
```bash
rsync \
	-ah \
	--info=stats2,progress2 \
	--no-owner \
	--no-group \
	--no-perms \
	--omit-dir-times \
	--exclude="*.fastq.gz" \
	/data/<project name> /mnt/storage3/promethion_rawdata_buffer/<YYYYmmdd>
```
In case the FASTQ output files are required (cDNA sequencing or
custom protocols), the exclude parameter (line 8) needs to be
removed. Sequencing data is scattered in many small files which are
concatenated, the limit for open files needs to be increased with
`ulimit` beforehand.




