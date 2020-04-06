# megSAP - DNA analysis (tumor samples, tumor/normal pairs)

### Introduction

Somatic analysis of tumor/normal DNA sample pairs and of tumor-only DNA samples is
performed with `somatic_dna.php`. It requires samples already aligned and present in
BAM format. Please have a look at the help:

    > php megSAP/src/Pipelines/somatic_dna.php --help

The main parameters that you have to provide are:

* `t_bam` - the tumor DNA sample BAM file
* `n_bam` - the normal DNA sample BAM file (only required for tumor/normal analysis)
* `-out_folder` - the output folder

Alignment of the individual DNA sample(s) can be performed with `analyze.php` using `-somatic`
and `-steps ma` options (also see [single sample DNA analysis](dna_single_sample.md)).


### Output

The main output files are listed below.

**Variant calling step (`vc`, `an`):**

* `somatic_var.vcf.gz` - SNVs and short indels called with Strelka2 (sample pairs) or FreeBayes (tumor sample only)
* `somatic_var_annotated.vcf.gz` and `somatic.GSvar` - annotated variants
* `somatic_var_structural.vcf.gz` and `somatic_var_structural.tsv` - structural variants called with manta
* `somatic_bafs.igv` - B-allele frequencies in IGV-compatible tabular format

**Copy-number variation calling step (`cn`):**

* `somatic_clincnv.tsv` - copy-number variations called with ClinCNV

**Microsatellite analysis step (`msi`):**

* `somatic_msi.tsv` - microsatellite instabilities called with MANTIS

### Additional annotation steps
**RNA annotation (`an_rna`):**
If you have data from the [RNA pipeline](../doc/rna_expression.md), this data can be annotated in the `an_rna` step. This will annotate read depth, allele frequency and TPM from the RNA sample to each variant contained in the GSvar file. Provide the bam file using the parameter `-t_rna_bam` and the RNA sample name using the parameter `-rna_id`. This step also annotates reference gene expression from The Human Protein Atlas. Therefore, you have to specify a reference tissue type. It can be passed via the parameter `-rna_ref_tissue`.

**CGI annotation (`ci`):**
This step will send the SNVs and CNVs from `somatic.GSvar` and `somatic_clincnv.tsv` to CancerGenomeInterpreter.org. The results will be downloaded into the out folder and annotated to both files. Provide a tumor type using the parameter `-cancer_type` for more accurate results, acronyms can be found on CancerGenomeInterpreter.org.


### Single Sample Tumor DNA Analysis

In tumor-only DNA analysis, variants are called with FreeBayes. Identification of somatic
variants out of all the variant calls is realized by applying filters on variant allele
frequencies from public databases (present in the annotated variant lists). This can be
achieved interactively in GSvar (filter panel) or by using `VariantFilterAnnotations`
(from `ngs-bits`) with the following filter specification file (for e.g. 1% filter):

```
Allele frequency    max_af=1
```

[back to the start page](../README.md)
