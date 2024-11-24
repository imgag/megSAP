# Building megSAP from sources (Linux)

Currently only Linux is supported!  

## Dependencies

We are providing instructions for Ubuntu 20.04 and RHEL 8.3 here. However this should be reasonably easy to port to any other Linux distribution.

### Base dependencies for different operating systems:

Ubuntu 20.04

	> sudo apt-get install -y rsync zlib1g-dev bzip2 build-essential php7.4-cli php7.4-xml php7.4-mysql make unzip wget git libssl-dev gnumeric tabix numdiff libdb-dev libgit2-dev libsqlite3-dev libxml2-utils pybind11-dev
    
Ubuntu 22.04

	> //TODO Kilian
    
Ubuntu 24.04

	> //TODO Kilian
    
## Install Apptainer

	> sudo add-apt-repository -y ppa:apptainer/ppa
	> sudo apt update
	> sudo apt install -y apptainer

## Downloading

Clone megSAP repository with containerized tools:

	> git clone -b containerize_tools https://github.com/imgag/megSAP.git

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
	> php ../src/Tools/db_download.php # DB downloads that require apptainer containers

**Note:** OMIM, HGMD and COSMIC are not downloaded automatically because of license issues. If you have the license for those databases, download/convert them according to the commented sections in the download script.

## NGSD initialization

//TODO Marc

**Note:** To annotate variants with NGSD in-house counts, classifications, etc., NGSD data has to be exported regularly. Adapt the file `data\dbs\NGSD\Makefile` and execute `make export` once a week using a cronjob.


## Settings

Changing the settings is not absolutely necessary as most entries have defaults set.  
If you want to change the settings, copy `settings.ini.default` to `settings.ini` and adapt the settings to your needs.  

Settings entries are described [here](settings.md)

## Execution

Now the pipelines with all required tools and data are installed. They can be found within the `src/Pipelines` folder. Go to the [documentation](../README.md) for further details.
