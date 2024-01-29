## megSAP and DRAGEN

The megSAP pipeline also supports mapping using the Illumina DRAGEN server. This requires additional software on the DRAGEN server and some modifications to the settings.

(For general setup instructions of DRAGEN use the [Illumina DRAGEN documentation](https://emea.support.illumina.com/sequencing/sequencing_software/dragen-bio-it-platform/documentation.html))

### Required software

To use the DRAGEN server in the megSAP pipeline some additional software has to be installed on the DRAGEN server:

* `php` - is required to run the mapping script which manages the file transfer to and from the DRAGEN server and starts the mapping. Additionally the modules `PDO` and `mysql` are required.
```
[sudo] yum install php php-pdo php-mysqlnd
```

* `SunGridEngine` - is required to queue and execute mapping jobs from other servers. For that a separate queue with 1 slot (on each DRAGEN server) has to be configured. This prevents the DRAGEN server to run more than one mapping job at a time.

* `megSAP` - has to be installed to the same path as it is installed on the server which runs the analysis. megSAP can be installed on a network share which is available on all analysis servers and DRAGEN.

### settings.ini

After the required software is installed, some parameters in the megSAP `settings.ini` have to be adapted: 

* `dragen_user` - User which is used to run the analysis on the DRAGEN server. It has to be the same user who started the complete analysis and has to have read and write access to the folders defined below.

* `dragen_in`/`dragen_out` - Transfer folders which have to be accessible from both the server which performs the main analysis and the DRAGEN server. These two folders are used to transfer data to the DRAGEN server, i.e. FASTQ files, and transfer data from the DRAGEN server to the analysis server, i.e. BAM files.

* `dragen_data` - Temporary folder on the DRAGEN server in which the mapping is performed. This folder should be located on the fast SSD storage of the DRAGEN server (usually: `/staging/...`) and is created for each mapping and deleted after the mapped data has been moved to the transfer folder.

* `dragen_genomes` - Path to the genome reference hash tables. Should also be stored on the DRAGEN SSD storage and has to have the following structure: A folder for each reference containing the FASTA file and its index. Additionally this folder has to contain a subfolder named `dragen` which contains the actual hash tables. A example of the folder structure for `GRCh38`: 
```
├── GRCh38
│   ├── dragen
│   │   ├── hash_table.cfg
│   │   ├── hash_table.cfg.bin
│   │   ├── hash_table.cmp
│   │   ├── hash_table_stats.txt
│   │   ├── reference.bin
│   │   ├── ref_index.bin
│   │   ├── repeat_mask.bin
│   │   ├── replay.json
│   │   ├── streaming_log.csv
│   │   └── str_table.bin
│   ├── GRCh38.fa
│   └── GRCh38.fa.fai
```
To create the hash table run the following command:
```
dragen --build-hash-table true --ht-reference /staging/genomes/GRCh38/GRCh38.fa --output-dir /staging/genomes/GRCh38/dragen/ --ht-suppress-mask yes --ht-suppress-decoys yes --enable-cnv true
```

* `dragen_log` - Folder to store STDOUT and STDERR of the queued DRAGEN mapping jobs to determine if a finished job has ended successfully.

* `queues_dragen` - Queue name(s) to which the DRAGEN mapping jobs are submitted. These queues must have 1 slot.
