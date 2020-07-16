## megSAP and DRAGEN
The megSAP pipeline also supports mapping using the illumina DRAGEN server. For that some additional software is required and some settings has to be adapted.

### Required software

To use the DRAGEN server in the megSAP pipeline some aditional software has to be installed on the DRAGEN server:

* `php` - is required to run the mapping script which manages the file transfer to and from the DRAGEN server and starts the mapping. Additionally the modules `PDO` and `mysql` are required.

* `Sun GridEngine` - is required to queue mapping jobs from other servers. For that a seperate queue with 1 slot (on each DRAGEN server) has to be configured. This prevents the DRAGEN server to run more than one mapping job at a time.

* `megSAP` - has to be installed to the same path as it is installed on the server which runs the analysis. Alternatively megSAP can be installed on a network share which is available DRAGEN and analysis server.

### settings.ini

After the required software is installed and configured probably some parameters in the megSAP `settings.ini` has to be adapted: 

* `dragen_user` - User which is used to run the analysis on the DRAGEN server. It has to be the same user who started the complete analysis and has to have read and write access to the folders defined below.

* `dragen_in`/`dragen_out` - Transfer folders which have to be accessible from both the server which performs the analysis and the DRAGEN server. These two foloders are used to transfer data to the DRAGEN server (e. g. FastQ files) and transfer data from the DRAGEN server to the analysis server (e. g. BAM files).

* `dragen_data` - Temporary folder on the DRAGEN server in which the mapping is performed. This folder should be located on the fast SSD storage of the DRAGEN server (usually: `/staging/...`) and is created for each mapping and deleted after the mapped data has been moved to the tranfer folder.

* `dragen_genomes` - Path to the genome reference hash tables. Should also be stored on the DRAGEN SSD storage.

* `queues_dragen` - Queue(s) with 1 slot where the DRAGEN mapping jobs are submitted to and are run on the DRAGEN server(s).