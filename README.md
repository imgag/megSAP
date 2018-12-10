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
The installaton of megSAP from sources is quite time-consumung.  
Therefore, a ready-to-use version of _megSAP_ is available via a Docker image:

- From **sources** for [Linux](doc/install_unix.md)
- **Docker image** for [Linux](doc/install_docker.md) (beta)


## Documentation

Documentation about the different analysis pipelines can be found here:

* [DNA analysis (single sample)](doc/dna_single_sample.md)
* [DNA analysis (multi-sample and trio)](doc/dna_multi_sample.md)
* [DNA analysis (tumor-normal pair)](doc/dna_tumor-normal_pair.md)
* [RNA analysis (expression)](doc/rna_expression.md)


## Support

Please report any issues or questions to the [megSAP issue 
tracker](https://github.com/imgag/megSAP/issues).


## ChangeLog

* 22.11.2018: Retrieving CGI cancer type from new NGSD entry (somatic pipeline)
* 16.11.2018: Updated to freebayes 1.2.0 (parallelized to compensate for the increased runtime)
* 07.11.2018: Updated VEP from version 93.2 to 94.5
* 19.10.2018: Reverted to freebayes 1.1.0 because version 1.2.0 is too slow.
* 10.10.2018: Updated all tools and databases to the latest releases.
* 04.10.2018: Using Ensembl VEP for variant annotation now (was SnpEff)
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



