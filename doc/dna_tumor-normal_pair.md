# megSAP - DNA analysis (tumor/normal pairs)

## Analysis pipeline

Somatic analysis of tumor/normal DNA sample pairs is performed with `somatic_dna.php`.

### Mapping single samples

The somatic pipeline requires requires tumor and normal sample BAM files.  
Thus, they need to be mapped first using the `analysis.php` with the parameter `-steps ma -somatic`.  
For more information, see the documentation at [single sample DNA analysis](dna_single_sample.md).

### Somatic variant calling (tumor/normal)

After mapping, the tumor and normal BAM files are used to perform somatic variant calling.  

The main parameters that you have to provide are:

* `t_bam` - the tumor DNA sample BAM file
* `n_bam` - the normal DNA sample BAM file (only required for tumor/normal analysis)
* `-out_folder` - the output folder

Please have a look at the help:

    > php megSAP/src/Pipelines/somatic_dna.php --help

### Tools used in this analysis pipeline

The tools used in this analysis can be found here: [Somatic tools](dna_somatic_tools.md).

### Performance

Performance benchmarks of the the megSAP pipeline can be found [here](performance_somatic.md).

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
