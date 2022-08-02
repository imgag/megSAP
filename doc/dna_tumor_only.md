# megSAP - DNA analysis (tumor only)

## Analysis pipeline

Somatic analysis of tumor-only DNA samples is performed with `somatic_dna.php`.

### Mapping tumor sample

The somatic pipeline requires requires the tumor sample BAM file.  
Thus, the tumor sample needs to be mapped first using the `analysis.php` with the parameter `-steps ma -somatic`.  
For more information, see the documentation at [single sample DNA analysis](dna_single_sample.md).

### Somatic variant calling

After mapping, the tumor BAM file is used to perform somatic variant calling.  

The main parameters that you have to provide are:

* `t_bam` - the tumor DNA sample BAM file
* `-out_folder` - the output folder

Please have a look at the help:

    > php megSAP/src/Pipelines/somatic_dna.php --help

Somatic variant calling without normal sample leads to a lot of false-positive somatic variant calls because artefacts cannot be removed efficiently.

To remove false-positives, filtering on variant allele frequencies from public databases can be done.  
This can be done interactively in GSvar or using `VariantFilterAnnotations` from `ngs-bits`.  
For 1% AF filter specify this line in the filter definition:

```
Allele frequency    max_af=1
```

### Tools used in this analysis pipeline

The tools used in this analysis can be found here: [Somatic tools](dna_somatic_tools.md).

## Output

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
If you have data from the [RNA pipeline](rna_expression.md), this data can be annotated in the `an_rna` step. This will annotate read depth, allele frequency and TPM from the RNA sample to each variant contained in the GSvar file. Provide the bam file using the parameter `-t_rna_bam` and the RNA sample name using the parameter `-rna_id`. This step also annotates reference gene expression from The Human Protein Atlas. Therefore, you have to specify a reference tissue type. It can be passed via the parameter `-rna_ref_tissue`.


[back to the start page](../README.md)
