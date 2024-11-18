# Building megSAP from sources (Linux)

Currently only Linux is supported!  

## Dependencies

We are providing instructions for Ubuntu 20.04 and RHEL 8.3 here. However this should be reasonably easy to port to any other Linux distribution.

### Base dependencies for different operating systems:

Ubuntu 20.04

	> sudo apt-get install -y rsync bzip2 default-jre bioperl libset-intervaltree-perl libjson-xs-perl libcarp-assert-perl libgd-dev libdb-dev libxml2-dev libxml2-utils php7.4-cli php7.4-xml php7.4-mysql tabix unzip wget build-essential cmake cpanminus git libbz2-dev liblzma-dev libncurses5-dev libqt5sql5-mysql libpng-dev libqt5xmlpatterns5-dev libssl-dev qt5-default qt5-qmake qtbase5-dev libcurl4-openssl-dev libhts-dev libtabixpp-dev libtabixpp0 meson ninja-build gnumeric numdiff libpcre2-dev libreadline-dev libffi-dev libharfbuzz-dev libfribidi-dev libgit2-dev pybind11-dev libsqlite3-dev openjdk-17-jdk openjdk-17-jre gfortran
    
Ubuntu 22.04

	> sudo apt-get install -y rsync bzip2 default-jre bioperl libset-intervaltree-perl libjson-xs-perl libcarp-assert-perl libgd-dev libdb-dev libxml2-dev libxml2-utils php8.1-cli php8.1-xml php8.1-mysql tabix unzip wget build-essential cmake cpanminus git libbz2-dev liblzma-dev libncurses5-dev libqt5sql5-mysql libpng-dev libqt5xmlpatterns5-dev libssl-dev qtbase5-dev qt5-qmake qtbase5-dev libhts-dev libtabixpp-dev libtabixpp0 meson ninja-build gnumeric numdiff libpcre2-dev libreadline-dev libffi-dev libharfbuzz-dev libfribidi-dev libgit2-dev pybind11-dev libsqlite3-dev openjdk-17-jdk openjdk-17-jre gfortran
    
RHEL 8.3

	Add LANGUAGE and LC_ALL to /etc/locale.conf
    
	Example:
	LANGUAGE="en_US.UTF-8"
	LANG="en_US.UTF-8"
	LC_ALL="en_US.UTF-8"

	> subscription-manager repos --enable codeready-builder-for-rhel-8-x86_64-rpms
	> yum groupinstall "Development Tools" -y
	> yum install zlib-devel bzip2-devel xz-devel ncurses-devel libcurl-devel cpan cpanminus gd-devel libdb-devel -y
	> dnf install php-cli php-xml  php-mysqlnd R-core R-core-devel -y
 	> yum install qt5-qtcharts.x86_64 qt5-qtbase-odbc.x86_64 qt5-qtbase-mysql.x86_64 qt5-qtxmlpatterns.x86_64 libcurl-devel.x86_64 qt5-devel.x86_64
    
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

To install the required tools and data you will need to execute custom script delivered with the repository.
If you work in a security critical environment it is advised that you use a [chroot environment](https://help.ubuntu.com/community/BasicChroot) as the scripts will attempt to install software with administrative privileges.

First, we make sure the privileges of the installation scripts are correct:

	> cd megSAP/data
	> chmod 755 *.sh

Next, we install all required tools and download apptainer containers of required tools:

	> ./download_tools.sh
	> ./download_container.sh

Next, we need to download and index the reference genome:
	
	> ./download_GRCh38.sh

Finally, we need to download and convert some open-source databases for annotations:

	> ./download_dbs.sh
	> ./download_dbs_rna.sh #only needed for RNA analysis
	
	(Some downloads require specific apptainer containers and are executed via the db_download.php script)
	> cd ..
	> php /src/Tools/db_download.php

**Note:** OMIM and HGMD are not downloaded automatically because of license issues. If you have the license for those databases, download/convert them according to the commented sections in the download script.

**Note:** The use of the optional NGSD annotation requires an export of the variants to VCF files. These files should be updated on a regular basis. For example code take a look at the NGSD section in the download script.


## Settings

Changing the settings is not absolutely necessary as most entries have defaults set.  
If you want to change the settings, copy `settings.ini.default` to `settings.ini` and adapt the settings to your needs.  

Settings entries are described [here](settings.md)

## Execution

Now the pipelines with all required tools and data are installed. They can be found within the `src/Pipelines` folder. Go to the [documentation](../README.md) for further details.
