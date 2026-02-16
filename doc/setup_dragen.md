## megSAP and DRAGEN

The megSAP pipeline also supports mapping using the Illumina DRAGEN server. This requires additional software on the DRAGEN server and some modifications to the settings.

(For general setup instructions of DRAGEN use the [Illumina DRAGEN documentation](https://emea.support.illumina.com/sequencing/sequencing_software/dragen-bio-it-platform.html))

**Note:** Currently on DRAGEN 4.3 is supported. Support for DRAGEN 4.4 is coming soon.

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

* `dragen_version` - Version of DRAGEN to use, e.g. `4.3.17`.

* `dragen_in`/`dragen_out` - Transfer folders which have to be accessible from both the server which performs the main analysis and the DRAGEN server. These two folders are used to transfer data to the DRAGEN server, i.e. FASTQ files, and transfer data from the DRAGEN server to the analysis server, i.e. BAM files.

* `dragen_data` - Folder used as working directory on the DRAGEN server. This folder should be located on `/staging/. A sub-directory is created for each DRAGEN analysis and deleted after the data has been moved to the transfer folder.

* `dragen_genome` - Path to the DRAGEN genome reference hash tables to use for the analysis. Should also be located on `/staging/`.

To create the hash tables in the folder `/staging/genomes/GRCh38/dragen/`, copy the GRCh38 reference genome to `/staging/genomes/GRCh38/GRCh38.fa` and run the following command:
```
dragen --build-hash-table true --ht-reference /staging/genomes/GRCh38/GRCh38.fa --output-dir /staging/genomes/GRCh38/dragen/ --enable-cnv true --ht-num-threads 40
```

* `dragen_log` - Folder to store STDOUT and STDERR of the queued DRAGEN mapping jobs to determine if a finished job has ended successfully.

* `queues_dragen` - Queue name(s) to which the DRAGEN mapping jobs are submitted. These queues must have 1 slot.


### Running DRAGEN analysis

The DRAGEN analysis is executed before the megSAP analysis.

If you run the DRAGEN analysis on NovaSeqX+ or on a AWS DRAGEN instance, the results of the DRAGEN analysis have to be copied into a sub-folder of the sample folder called `dragen`.  
Then the megSAP germline analysis is run with all steps but the mapping step.

If you use a on-premise DRAGEN via a SGE queue, the DRAGEN analysis is performed with `src/Pipelines/analyze_dragen.php`.  
It automatically starts the megSAP analysis, unless the parameter `-no_queuing` is used.
