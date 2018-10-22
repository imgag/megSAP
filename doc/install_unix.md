# Building megSAP from sources (Linux)

Currently only Linux is supported!  

It is theoretically feasible to get this running on OSX, but we are not testing this regularily.  
If you do want to work on adding OSX support please contact [Lennard Berger](https://github.com/Fohlen).

## Dependencies

We are providing instructions for Ubuntu 16.04 and 18.04 officially. However this should be reasonably easy to port to any other Linux distribution.

### Base dependencies (Ubuntu 16.04)

```
    > sudo apt-get update && apt-get install -y
    bzip2 \
    default-jre \
    perl-base \
    php7.0-cli \ 
    php7.0-xml \ 
    php7.0-mysql \
    python-matplotlib \ 
    python-numpy \
    python-pysam \
    tabix \
    unzip \
    wget \
    build-essential \ 
    cmake \ 
    cpanminus \
    git \ 
    libbz2-dev \ 
    liblzma-dev \ 
    libncurses5-dev \ 
    libpng-dev \ 
    libqt5sql5-mysql \ 
    libqt5xmlpatterns5-dev \ 
    libssl-dev \ 
    qt5-default \ 
    qt5-qmake \ 
    qtbase5-dev 
```

### Base dependencies (Ubuntu 18.04)

```
    > sudo apt-get update && apt-get install -y \
    bzip2 \
    default-jre \
    perl-base \ 
    php7.2-cli \ 
    php7.2-xml \ 
    php7.2-mysql \ 
    python-matplotlib \ 
    python-numpy \
    python-pysam \
    tabix \ 
    unzip \ 
    wget \
    build-essential \ 
    cmake \ 
    cpanminus \
    git \ 
    libbz2-dev \ 
    liblzma-dev \ 
    libncurses5-dev \ 
    libqt5sql5-mysql \ 
    libpng-dev \
    libqt5xmlpatterns5-dev \ 
    libssl-dev \
    qt5-default \ 
    qt5-qmake \ 
    qtbase5-dev
```

### R dependencies (for somatic CNV calling)

	> sudo apt install r-base-core r-cran-optparse r-cran-robustbase r-cran-foreach r-cran-doparallel r-cran-mass r-cran-data.table

### Perl dependencies (for RNA pipeline)

	> sudo -E cpanm Set::IntervalTree URI::Escape DB_File Carp::Assert JSON::XS PerlIO::gzip

## Downloading

If you haven't already, check out the repository with

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

To install the required tools and data you will need to execute custom script delivered with the repository.
If you work in a security critical environment it is advised that you use a [chroot environment](https://help.ubuntu.com/community/BasicChroot) as the scripts will attempt to install software with administrative privileges.

First, we make sure the privileges of the installation scripts are correct:

	> cd megSAP/data
	> chmod 755 *.sh

Next, we install all required tools

	> ./download_tools.sh
	> ./download_tools_somatic.sh #only needed for tumor analysis
	> ./download_tools_rna.sh #only needed for RNA analysis

Next, we need to download and index the reference genome:
	
	> ./download_GRCh37.sh


Finally, we need to download and convert some open-source databases for annotations:

	> ./download_dbs.sh
	> ./download_dbs_rna.sh #only needed for RNA analysis

**Note:** OMIM, HGMD and COSMIC are not downloaded automatically because of license issues. If you have the license for those databases, download/convert them according to the commented sections in the download script.

## Execution

Now the pipelines with all required tools and data are installed. They can be found within the `src/Pipelines` folder. Go to the [documentation](../README.md) for further details.






