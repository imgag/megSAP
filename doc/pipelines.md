# Pipelines

To run a pipeline the following prerequisites must be met:

- properly downloaded and installed _megSAP_ pipeline
- downloaded the reference genome and databases
- use a `settings.ini` so megSAP finds all of these.

**Note:** _megSAP_ uses a fallback `settings.ini` that is optimised to work with the pre-shipped version. Have a look on the [default settings](../settings.ini.default) if you need to customize.

Please refer to the [installation](./install_docker.md) if any of these does not yet apply.

Below is an example how to run the `RNA` pipeline

```
docker run -v /Sample_X_folder:/Sample_X_01/ -v $(pwd)/genomes:/home/ubuntu/megSAP/data/genomes -v $(pwd)/dbs:/home/ubuntu/megSAP/data/dbs 
  -it imgag/megSAP php src/Pipelines/analyze_rna.php \
  -folder Sample_X_01 -name X_01 \
  -system truseq.ini -steps ma,rc,an
```

This will

- mount the Sample_X_folder to the path supplied to the pipeline
- mount the reference databases to `/megSAP/data/genomes`
- finally run the `analyze_rna` pipeline on Sample_X_folder

## Documentation

Documentation about the different analysis pipelines can be found here:

* [DNA analysis (single sample)](doc/dna_single_sample.md)
* [DNA analysis (multi-sample and trio)](doc/dna_multi_sample.md)
* [DNA analysis (tumor-normal pair)](doc/dna_tumor-normal_pair.md)
* [RNA analysis (expression)](doc/rna_expression.md)
* RNA analysis (variant calling)  - coming soon
