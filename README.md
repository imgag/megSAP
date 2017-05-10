# *megSAP* - a Medical Genetics Sequence Analysis Pipeline

The goal of the megSAP proejct is to develop an NGS data analysis pipeline for medical genetics, that is

 * state-of-the-art in terms of performance,
 * fast
 * and usable for diagnostics:
 	* we use no tools that require a license for diagnostics e.g. [GATK](https://software.broadinstitute.org/gatk/)
    * extensive logging (tools, versions, parameters) ensures reproducability of results
	* extensive testing before adding/updating tools or databases makes sure the results are valid

megSAP is developed by the [Institute of Medical Genetics and Applied Genomics](http://www.uni-tuebingen.de/Klinische_Genetik/start.html) and several collaborators from academia and industry. If you are interested to join the effort, please contact [Marc Sturm](https://github.com/marc-sturm).

A poster about megSap which was presented at NGS'17 in barcelona can be downloaded 
[here](doc/Poster_April_2017.pdf)

## Download

There are no releases yet.  
Please use git to download the most recent development version:

    git clone --recursive https://github.com/imgag/megSAP.git

### Resolving proxy issues with git

If you are behind a proxy that block the standard git port, you see something like this:

    $ git clone --recursive https://github.com/imgag/megSAP.git
    Cloning into 'megSAP'...
    fatal: Unable to look up github.com (port 9418) (Name or service not known)

Then you have to adapt your ~/.gitconfig file like that:

    [http]
    proxy = http://[user]:[password]@[host]:[port]

## Dependencies

megSAP (or one of the used tools) depends mainly on the following software to be installed:

* _g++_ (4.5 or higher)
* _qmake_ (Qt 5.5 or higher, including xmlpatterns and mysql package)
* _git_
* _cmake_
* _python_ (including matplotlib)
* _php_
* _tabix_

For example, the installation of the dependencies using Ubuntu 16.04 looks like that:

	> sudo apt-get install -y wget bzip2 unzip make g++ git cmake tabix build-essential qt5-default qt5-qmake qtbase5-dev libqt5sql5-mysql libqt5xmlpatterns5-dev php7.0-cli php7.0-xml php7.0-mysql python python-matplotlib libncurses5-dev

## Initial setup

First, we need to download the tools the pipeline relies on:

	> cd megSAP/data/
	> chmod 755 download_*.sh
	> ./download_tools.sh

Next, we need to download and index the reference genome:
	
	> ./download_GRCh37.sh


Finally, we need to download and convert some open-source databases for annotations:

	> ./download_dbs.sh

**Note:** OMIM, HGMD and COSMIC are not downloaded automatically because of license issues. If you have the license for those databases, download/convert them according to the commented sections in the download script.

## Running an analysis (single sample DNA)

The analysis pipeline assumes that that all data to analyze resides in a sample folder as produced by Illumina's [bcl2fastq](http://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html) tool. If that is the case, the whole analysis is performed with one command.  
For example like this:

	php megSAP/src/Pipelines/analyze.php -folder Sample_NA12878_01 -name NA12878_01 -system hpHBOCv5.ini -steps ma,vc,an

In the example above, the configuration of the pipeline is done using the `hpHBOCv5.ini` file, which contains information about adapter sequences, target region, data type, etc.

Example data with can be analyzed using the command above can be downloaded from [here](https://medgen.medizin.uni-tuebingen.de/NGS-downloads/NA12878_01.zip).

After the data analysis, the sample folder contains BAM and VCF (gzipped) files as expected. Additionally, several qcML files that contain QC data are generated (open with a browser).

## Documentation

Detailed documentation of the analysis pipeline can be found here:

* [DNA analysis (single sample)](doc/dna_single_sample.md)
* DNA analysis (trios and families) - coming soon
* DNA analysis (tumor-normal pair) - coming soon
* RNA analysis (expression) - coming soon
* RNA analysis (variant calling)  - coming soon


## Support

Please report any issues or questions to the [megSAP issue 
tracker](https://github.com/imgag/megSAP/issues).











