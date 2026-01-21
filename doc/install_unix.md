# Building megSAP from sources (Linux)

Currently only Linux is supported!  

*Note: Indexing the genome with bwa-mem2 requires about 60BG of RAM. Do not try to install megSAP on a machine with less than 64GB of RAM!*

## Base dependencies

We are providing instructions for the latest Ubuntu LTS distibution (Ubuntu 24.04).  
If you are using other Linux distributions, you have to adapt them yourself.

	> sudo apt-get update
	> sudo apt-get install -y rsync zlib1g bzip2 php8.3-cli php8.3-xml php8.3-mysql make unzip wget git gnumeric pigz ghostscript

Install Apptainer:

	> sudo add-apt-repository -y ppa:apptainer/ppa
	> sudo apt update
	> sudo apt install -y apptainer

## Downloading

Clone the last release of megSAP:

	> git clone --branch 2025_10 https://github.com/imgag/megSAP.git

### Resolving proxy issues with git

If you are behind a proxy that block the standard git port, you see something like this:

    $ git clone https://github.com/imgag/megSAP.git
    Cloning into 'megSAP'...
    fatal: Unable to look up github.com (port 9418) (Name or service not known)

Then you have to adapt your ~/.gitconfig file like that:

    [http]
    proxy = http://[user]:[password]@[host]:[port]

## Initial setup

To install the required tools and databases you will need to execute some installation scripts in the order described here.

First, we make sure the privileges of the installation scripts are correct:

	> cd megSAP/data
	> chmod 755 *.sh

Next, we install a few tools and download apptainer containers for the rest of the tools:

	> ./download_tools.sh
	> ./download_container.sh

Next, we need to download and index the reference genome:
	
	> ./download_GRCh38.sh

Finally, we need to download and convert some open-source databases for annotations:

	> ./download_dbs.sh
	> php ../src/Install/db_download.php # DB downloads that require apptainer containers

**Note:** OMIM, HGMD and COSMIC are not downloaded automatically because of license issues. If you have the license for those databases, download/convert them according to the commented sections in the `download_dbs.sh` script.

## NGSD initialization

If you want to use the NGSD and it is not initialized already, perform the following steps:

1) Install the database server and perform initial import of gene, transcript and disease data as described [here](https://github.com/imgag/ngs-bits/blob/master/doc/install_ngsd.md).
2) Add the NGSD credentials to the megSAP `settings.ini` file.

**Note:** To annotate variants with NGSD in-house counts, classifications, etc., NGSD data has to be exported regularly. To do so, adapt the file `data\dbs\NGSD\Makefile` and execute `make export` once a week using a cronjob.


## Settings

Changing the settings is not absolutely necessary as most entries have defaults set.  
If you want to change the settings, copy `settings.ini.default` to `settings.ini` and adapt the settings to your needs.  

Settings entries are described [here](settings.md)

## Execution

Now the pipelines with all required tools and data are installed. They can be found within the `src/Pipelines` folder. Go to the [documentation](../README.md) for further details.

