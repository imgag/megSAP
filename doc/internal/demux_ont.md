# ONT PromethION Data Transfer and Post-Run Process

1. Check run quality using MinKNOW UI software

    Basecalled bases >50Gb (good), 30 to 50Gb (medium), <30GB (bad)   
    Approximate N50 >10kb (human WGS)   

2. Set the run quality in the NGSD 

	Use comment field to describe problems.

3. Data transfer

    Raw data is stored on the PromethION in `/data/<project name>`.
   
    Connect to the device via ssh:
    ```bash
    ssh minerva
    ```

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

4. Create sample raw data

    Using `copy_ont_data.php`, the flow cell data is copied to the
    default sample folder location, based on the information entered in
    GSvar.

    ```bash
    ulimit -n 10000
    php /mnt/storage2/megSAP/src/NGS/copy_ont_run.php \
        -run_name <run name> \
        -run_dir <run directory> \
        -bam \
        -ignore_aligned \
        -queue_sample
    ```

    The result file is `<processed sample>.mod.unmapped.bam` in the
    corresponding project/sample directory and the sample analysis is
    automatically queued.

    The script aborts if the flow cell directory contains a `pod5_skip`
    subdirectory, indicating a partially basecalled flow cell (see
    below).

5. Backup raw run data using backup script:
	```bash
	sudo -u archive-gs php /mnt/storage2/megSAP/pipeline/src/Tools/backup_queue.php -mode run -in [run] -email [email]
	```

6. Delete the run raw data (when all samples are analyzed with passed QC):
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
    --min-qsore 9 \
    > <run directory>/<flow cell directory>/bam_pass/pod5_skip_basecalled.bam
```

Rename the `pod5_skip` directory to continue with the `copy_ont_run.php`
script.

```bash
cd <run directory>/<flow cell directory>
mv pod5_skip pod5_skip_done
```