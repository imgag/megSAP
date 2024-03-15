# Building megSAP from sources (Linux)

Currently only Linux is supported!  

## Dependencies

We are providing instructions for Ubuntu 20.04 and RHEL 8.3 here. However this should be reasonably easy to port to any other Linux distribution.

### Base dependencies

Ubuntu 20.04

	> sudo apt-get install -y rsync bzip2 default-jre bioperl libset-intervaltree-perl libjson-xs-perl libcarp-assert-perl libgd-dev libdb-dev libxml2-dev libxml2-utils php7.4-cli php7.4-xml php7.4-mysql tabix unzip wget build-essential cmake cpanminus git libbz2-dev liblzma-dev libncurses5-dev libqt5sql5-mysql libpng-dev libqt5xmlpatterns5-dev libssl-dev qt5-default qt5-qmake qtbase5-dev libcurl4-openssl-dev libhts-dev libtabixpp-dev libtabixpp0 meson ninja-build gnumeric numdiff libpcre2-dev libreadline-dev libffi-dev libharfbuzz-dev libfribidi-dev libgit2-dev pybind11-dev libsqlite3-dev
    
RHEL 8.3

	> add LANGUAGE and LC_ALL to /etc/locale.conf
    
	Example:
	LANGUAGE="en_US.UTF-8"
	LANG="en_US.UTF-8"
	LC_ALL="en_US.UTF-8"

	> subscription-manager repos --enable codeready-builder-for-rhel-8-x86_64-rpms
	> yum groupinstall "Development Tools" -y
	> yum install zlib-devel bzip2-devel xz-devel ncurses-devel libcurl-devel cpan cpanminus gd-devel libdb-devel -y
	> dnf install php-cli php-xml  php-mysqlnd R-core R-core-devel -y
 	
  	# dependencies for ngs-bits:
 	> yum install qt5-qtcharts.x86_64 qt5-qtbase-odbc.x86_64 qt5-qtbase-mysql.x86_64 qt5-qtxmlpatterns.x86_64 libcurl-devel.x86_64 qt5-devel.x86_64
    
## Downloading

Clone the latest release of megSAP:

	> git clone -b 2023_11 https://github.com/imgag/megSAP.git

Or, if you want to test the current development version:

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
	> ./download_tools_somatic.sh #only needed for somatic analysis
	> ./download_tools_rna.sh #only needed for RNA analysis

Next, we need to download and index the reference genome:
	
	> ./download_GRCh38.sh


Finally, we need to download and convert some open-source databases for annotations:

	> ./download_dbs.sh
	> ./download_dbs_rna.sh #only needed for RNA analysis

**Note:** OMIM and HGMD are not downloaded automatically because of license issues. If you have the license for those databases, download/convert them according to the commented sections in the download script.

**Note:** The use of the optional NGSD annotation requires an export of the variants to VCF files. These files should be updated on a regular basis. For example code take a look at the NGSD section in the download script.

## Execution

Now the pipelines with all required tools and data are installed. They can be found within the `src/Pipelines` folder. Go to the [documentation](../README.md) for further details.
