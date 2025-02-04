# *megSAP* - Medical Genetics Sequence Analysis Pipelines

megSAP is a NGS data analysis pipeline for medical genetics, which is developed by the [Institute of Medical Genetics and Applied Genomics, University Hospital, Tübingen](http://www.uni-tuebingen.de/Klinische_Genetik/start.html) and several collaborators from academia and industry.

The design goals of the project are:

 * state-of-the-art in terms of sensitivity/specificity,
 * fast
 * and usable for diagnostics:
 	* we use no tools that require a license for diagnostics
    * extensive logging (tools, versions, parameters) ensures reproducability of results
	* extensive testing before adding/updating tools or databases makes sure the results are valid

If you are interested to join the effort, please contact [Marc Sturm](https://github.com/marc-sturm).

## Obtaining megSAP

General remarks:
- Only GRCh38 is supported. There is a [branch](https://github.com/imgag/megSAP/tree/GRCh37) for GRCh37, but is is not updated or maintained since December 2021 anymore.
- The installation of megSAP is quite time-consuming because large databases for annotation of variants need to be downloaded and converted.

The default way of using megSAP is cloning the megSAP repository and calling the analysis pipelines from there.
Installation instructions for this way can be found [here](doc/install_unix.md).  

Alternatively, there is a *beta* version of the megSAP pipeline in a container.  
You still need to download databases and tools, but the container version may be more convenient in some scenarios like AWS.  
The installation instructions for the container version of megSAP can be found [here](doc/install_unix_container_version.md).

## Documentation

Documentation about the different **Illumina short-read pipelines** can be found here:

* [DNA germline analysis (single sample)](doc/dna_single_sample.md)
* [DNA germline analysis (multi-sample and trio)](doc/dna_multi_sample.md)
* [DNA somatic analysis (tumor-normal pair)](doc/dna_tumor-normal_pair.md)
* [DNA somatic analysis(tumor only)](doc/dna_tumor_only.md)
* [RNA analysis (expression, fusions)](doc/rna_expression.md)

Documentation about the different **ONT long-read pipelines** can be found here:

* [DNA germline analysis (single sample)](doc/dna_longread_single_sample.md)

## Support

Please report any issues or questions to the [megSAP issue tracker](https://github.com/imgag/megSAP/issues).

## Citing

You can cite megSAP using Zenodo DOIs:

* 2024_09: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13744183.svg)](https://doi.org/10.5281/zenodo.13744183)

## ChangeLog

Major changes of master since last release:

* tools are not longer installed, but used from pre-built apptainer containers.
* germline DNA: improved speed (C++ tool) and IO (gzipped coverage profiles) of CNV reference cohort selection.

For older changes see [releases](https://github.com/imgag/megSAP/releases).
