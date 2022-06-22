# megSAP - RNA Expression Analysis

### Basics

Single sample RNA expression analysis is performed using the `analyze_rna.php`
script. Please have a look at the help using:

```shell
php megSAP/src/Pipelines/analyze_rna.php --help
```

The main parameters that you have to provide are:

* `folder` - The output folder containing all result files.
* `name` - Basename/prefix for all output files.

In addition, you may want to specify:

* `steps` - Analysis steps to perform. Use `ma,rc,an,fu` to perform
   mapping, read counting, annotation and fusion detection.
* `system` - The [processing system INI file](processing_system_ini_file.md).

### Running an analysis

If all data to analyze resides in a sample folder as produced by Illumina's
[bcl2fastq](http://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html)
tool, the whole analysis is performed with one command, for example like this:

```shell
php megSAP/src/Pipelines/analyze_rna.php \
  -folder Sample_X_01 -name X_01 \
  -system truseq.ini -steps ma,rc,an
```

In the example above, the configuration of the pipeline is done using the
`truseq.ini` file, which contains all necessary information (see [processing
system INI file](processing_system_ini_file.md)).


### Tools used in this analysis pipeline

The tools used in this analysis pipline can are described [here](rna_expression_pipeline.md).

A complete list of all tools and databases used in megSAP and when they were last updated can be found [here](update_overview.md).

### Output

After the analysis, these files are created in the output folder:

1. mapped reads in BAM format
2. raw read counts, in [featureCounts](http://bioinf.wehi.edu.au/featureCounts/)
   tabular output format
3. normalized read counts, annotated with gene symbols
4. QC data in [qcML format](https://www.ncbi.nlm.nih.gov/pubmed/24760958), which
   can be opened with a web browser


### Using other genomes

To use the RNA expression pipeline with other genomes, you need to provide

* the genome FASTA file, e.g. `megSAP/data/genomes/CustomGenome.fa`
* the STAR genome index, e.g. `megSAP/data/genomes/STAR/CustomGenome/`
* the gene annotation file in Ensembl-like GTF format, e.g.
  `megSAP/data/dbs/gene_annotations/CustomGenome.gtf`

The genome can by specified in the [processing system INI
file](processing_system_ini_file.md).

[back to the start page](../README.md)

