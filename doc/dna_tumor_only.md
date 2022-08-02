# megSAP - DNA analysis (tumor only)

## Analysis pipeline

Somatic analysis of tumor-only DNA samples is performed with `somatic_dna.php`.

### Mapping tumor sample

The somatic pipeline requires requires the tumor sample BAM.
Thus, it needs to be mapped first using the `analysis.php` with the parameter `-steps ma -somatic`.  
For more information, see the documentation at [single sample DNA analysis](dna_single_sample.md).

### Somatic variant calling

Somatic variant calling without normal sample leads to a lot of false-positive somatic variant calls because artefacts cannot be removed efficiently.  
To use the tumor-only pipeline, simply skip the `n_bam` parameter of `somatic_dna.php`.

To remove false-positives, filtering on variant allele frequencies from public databases can be done.  
This can be done interactively in GSvar or using `VariantFilterAnnotations` from `ngs-bits`.  
For 1% AF filter specify this line in the filter definition:

```
Allele frequency    max_af=1
```


### Tools used in this analysis pipeline

The following tools are used for mapping and calling of small variants and annotation of small variants:

| step                                           | tool                     | version              | comments                                         |
|------------------------------------------------|--------------------------|----------------------|--------------------------------------------------|
| QC - sample gender verification                | SampleGender             | ngs-bits latest      |                                                  |
| QC - low coverage statistics                   | BedLowCoverage           | ngs-bits latest      |                                                  |
| QC - low coverage annotation                   | BedAnnotateGenes         | ngs-bits latest      |                                                  |
| Variant calling - hla                          | HLA-genotyper            | 2022_05              |   https://github.com/axelgschwind/hla-genotyper  |
| Variant calling - SNVs and Indels              | Varscan2                 | 2.4.4                |                                                  |
| Variant calling - b-allele frequency (t/n)     | VariantAnnotateFrequency | ngs-bits latest      |                                                  |
| Variant calling - left-normalization of InDels | VcfLeftNormalize         | ngs-bits latest      |                                                  |
| Annotation                                     | VEP                      | 104.3                |                                                  |

CNV calling and annotation is performed using these tools:

| step                                               | tool                 | version              | comments                                            |
|----------------------------------------------------|----------------------|----------------------|-----------------------------------------------------|
| CNV calling                                        | ClinCNV              | 1.17.2               |                                                     |
| annotation - general                               | BedAnnotateFromBed   | ngs-bits latest      | Several data sources are annotated using this tool. |
| annotation - gene information                      | CnvGeneAnnotation    | ngs-bits latest      |                                                     |
| annotation - overlapping pathogenic CNVs from NGSD | NGSDAnnotateCNV      | ngs-bits latest      |                                                     |

SV calling and annotation is performed using these tools:

| step                                      | tool                            | version              | comments                                            |
|-------------------------------------------|---------------------------------|----------------------|-----------------------------------------------------|
| SV calling                                | Manta                           | 1.6.0                |                                                     |
| annotation - gene information             | BedpeGeneAnnotation             | ngs-bits latest      |                                                     |

Viral load:

| step                                           | tool                     | version              | comments                                            |
|------------------------------------------------|--------------------------|----------------------|-----------------------------------------------------|
| Filtering reads                                | samtools                 | 1.15.1               |                                                     |
| Mapping                                        | bwa mem2                 | 2.2.1                |                                                     |
| Variant calling                                | freebayes                | 1.3.6                |                                                     |
| Statistics - Coverage                          | BedCoverage              | ngs-bits latest      |                                                     |
| Statistics - BedReadCount                      | BedReadCount             | ngs-bits latest      |                                                     |
| Annotation                                     | BedAnnotateFromBed       | ngs-bits latest      | Several data sources are annotated using this tool. |

Microsatellite instability (MSI):

| step                                           | tool                     | version              | comments                                            |
|------------------------------------------------|--------------------------|----------------------|-----------------------------------------------------|
| Calling                                        | Mantis                   | 1.0.5                |                                                     |



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
