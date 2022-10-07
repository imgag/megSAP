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

## Obtaining megSAP

The installation of megSAP is quite time-consuming because many tools and big databases need to installed.  
The installation instructions can be found [here](doc/install_unix.md).

*Note: Only GRCh38 is supported. There is a [branch](https://github.com/imgag/megSAP/tree/GRCh37) for GRCh37, but is is not updated or maintained since December 2021 anymore.*

## Documentation

Documentation about the different analysis pipelines can be found here:

* [DNA germline analysis (single sample)](doc/dna_single_sample.md)
* [DNA germline analysis (multi-sample and trio)](doc/dna_multi_sample.md)
* [DNA somatic analysis (tumor-normal pair)](doc/dna_tumor-normal_pair.md)
* [DNA somatic analysis (tumor only)](doc/dna_tumor_only.md)
* [RNA analysis (expression)](doc/rna_expression.md)


## Support

Please report any issues or questions to the [megSAP issue 
tracker](https://github.com/imgag/megSAP/issues).


## ChangeLog

Changes of master since last release:

* fixed false duplications in GRCh38 (see [details](https://github.com/imgag/megSAP/blob/master/doc/internal/grch38_false_duplications.md)).
* germline DNA pipeline: Improved sensitivity of non-diploid variant calling (mitochondrial, mosaic).
* germline DNA pipeline: Added calling of mitochondrial variants for shallow WGS analysis.
* germline DNA pipeline: Added small variant and structural variant calling on Dragen when '-use_dragen' is enabled (see [details](https://github.com/imgag/megSAP/blob/master/doc/internal/vc_with_dragen.md)).
* updated VEP to version 107.
* updated ngs-bits to version 2022_10 (improved/added multithreading in several tools).
* updates ClinCNV to version 1.18.0 (improved handling of PAR region).
* minor fixes and updates.

For older changes see [releases](https://github.com/imgag/megSAP/releases).
