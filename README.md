# *megSAP* - a Medical Genetics Sequence Analysis Pipeline

megSAP is a NGS data analysis pipeline for medical genetics, which is developed by the [Institute of Medical Genetics and Applied Genomics](http://www.uni-tuebingen.de/Klinische_Genetik/start.html) and several collaborators from academia and industry. Since December 2016 the project is publicly available on GitHub, however closed-source development started already in 2012.  

The design goals of the project are:

 * state-of-the-art in terms of sensitivity/specificity,
 * fast
 * and usable for diagnostics:
 	* we use no tools that require a license for diagnostics
    * extensive logging (tools, versions, parameters) ensures reproducability of results
	* extensive testing before adding/updating tools or databases makes sure the results are valid

If you are interested to join the effort, please contact [Marc Sturm](https://github.com/marc-sturm).

## ChangeLog

* TODO.2018: Using Ensembl VEP for variant annotation now (was SnpEff)
* 07.08.2018: Removed most annotation from the 'filter' column and moved the functionality to GSvar.
* 20.07.2018: Updated ABRA2 version (attention: this changes indel positions - see ABRA2 2.06 changelog) 
* 17.07.2018: Added '##PIPELINE' header line to GSvar files to keep track of the megSAP version the file was created with.
* 11.07.2018: Added UPD detection for trios.
* 22.04.2018: NGSD import of germline variants restricted to variants with AF<5% to improve performance.
* 20.04.2018: NGSD now handles analysis job queuing and execution on SGE.
* 13.03.2018: Refactoring of trio analysis: it is now based on multi-sample pipeline, produces a GSvar file with three sample columns, and calls off-target variants.
* 18.01.2018: Added b-allele frequency files for visualization in IGV. 
* 08.01.2018: Updated tools (BWA, samtools, snpEff) and databases (ClinVar, gnomAD, HGMD).
* 15.12.2017: Added runs-of-homozygosity detection to the germline pipeline.

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
* _python_ (including matplotlib)
* _php_
* _tabix_
* _java_ 

For example, the installation of the dependencies using Ubuntu 16.04 looks like that:

	> sudo apt-get install -y wget bzip2 unzip make cmake g++ git tabix build-essential qt5-default qt5-qmake qtbase5-dev libqt5sql5-mysql libqt5xmlpatterns5-dev php7.0-cli php7.0-xml php7.0-mysql python python-matplotlib libncurses5-dev bzip2 libbz2-dev liblzma-dev default-jre libssl-dev libpng-dev perl-base curl mysql-client 

For example, the installation of the dependencies using Ubuntu 18.04 looks like that:

	> sudo apt-get install -y wget bzip2 unzip make cmake g++ git tabix build-essential qt5-default qt5-qmake qtbase5-dev libqt5sql5-mysql libqt5xmlpatterns5-dev php7.2-cli php7.2-xml php7.2-mysql python python-matplotlib libncurses5-dev bzip2 libbz2-dev liblzma-dev default-jre libssl-dev libpng-dev perl-base curl libmysqlclient-dev cpanminus

For molecular barcode processing, several python dependencies are required. They can be installed with ``pip``:

	> sudo apt-get install python-pip
	> pip install -r data/python_requirements.txt

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

## Documentation

Documentation about the different analysis pipelines can be found here:

* [DNA analysis (single sample)](doc/dna_single_sample.md)
* [DNA analysis (multi-sample and trio)](doc/dna_multi_sample.md)
* [DNA analysis (tumor-normal pair)](doc/dna_tumor-normal_pair.md)
* [RNA analysis (expression)](doc/rna_expression.md)
* RNA analysis (variant calling)  - coming soon


## Support

Please report any issues or questions to the [megSAP issue 
tracker](https://github.com/imgag/megSAP/issues).


