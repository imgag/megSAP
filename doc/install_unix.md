# Building megSAP from sources (Linux)

Currently only Linux is supported!  

## Base dependencies

We are providing instructions for the latest Ubuntu LTS distibutions here.  
If you are using other Linux distributions, you have to adapt them yourself.
    
Ubuntu 22.04

	> sudo apt-get update
	> sudo apt-get install -y rsync zlib1g bzip2 php8.1-cli php8.1-xml php8.1-mysql make unzip wget git gnumeric pigz ghostscript

Ubuntu 24.04

	> sudo apt-get update
	> sudo apt-get install -y rsync zlib1g bzip2 php8.3-cli php8.3-xml php8.3-mysql make unzip wget git gnumeric pigz ghostscript

Install Apptainer:

	> sudo add-apt-repository -y ppa:apptainer/ppa
	> sudo apt update
	> sudo apt install -y apptainer

## Downloading

Clone megSAP repository with containerized tools:

	> git clone https://github.com/imgag/megSAP.git

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

1) Install a MariaDB server
2) Create a database and a associated user in the SQL database.
3) Add the NGSD information to the megSAP `settings.ini` file.
4) Create tables using the following script:

	> php ../src/Install/db_init.php

5) Import genomics base data (genes, transcripts, phenotypes, gene-phenotype associations, ...) using the following tools from ngs-bits:

	> NGSDImportQC --help  
	> NGSDImportHGNC --help  
	> NGSDImportEnsembl --help  
	> NGSDImportHPO --help  
	> NGSDImportGeneInfo --help  
	> NGSDImportOMIM --help  
	> NGSDImportORPHA --help  

**Note:** To call ngs-bits tools, you have to call the apptainer container like that `apptainer exec data/tools/apptainer_container/ngs-bits_[version].sif [tool] [parameters]`.

**Note:** To annotate variants with NGSD in-house counts, classifications, etc., NGSD data has to be exported regularly. To do so, adapt the file `data\dbs\NGSD\Makefile` and execute `make export` once a week using a cronjob.


## Settings

Changing the settings is not absolutely necessary as most entries have defaults set.  
If you want to change the settings, copy `settings.ini.default` to `settings.ini` and adapt the settings to your needs.  

Settings entries are described [here](settings.md)

## Execution

Now the pipelines with all required tools and data are installed. They can be found within the `src/Pipelines` folder. Go to the [documentation](../README.md) for further details.
