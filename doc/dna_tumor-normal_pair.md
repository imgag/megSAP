# megSAP - DNA analysis (tumor/normal pairs)

## Analysis pipeline

Somatic analysis of tumor/normal DNA sample pairs is performed with `somatic_dna.php`.

### Mapping single samples

The somatic pipeline requires requires tumor and normal sample BAM files.  
Thus, they need to be mapped first using the `analysis.php` with the parameter `-steps ma -somatic`.  
For more information, see the documentation at [single sample DNA analysis](dna_single_sample.md).

### Somatic variant calling (tumor/normal)

After mapping, the tumor and normal BAM files are used to perform somatic variant calling withformat.  

The main parameters that you have to provide are:

* `t_bam` - the tumor DNA sample BAM file
* `n_bam` - the normal DNA sample BAM file (only required for tumor/normal analysis)
* `-out_folder` - the output folder

Please have a look at the help:

    > php megSAP/src/Pipelines/somatic_dna.php --help

### Tools used in this analysis pipeline

The following tools are used for mapping and calling of small variants and annotation of small variants:

| step                                           | tool                     | version              | comments                                         |
|------------------------------------------------|--------------------------|----------------------|--------------------------------------------------|
| QC - sample correlation verification           | SampleSimilarity         | ngs-bits latest      |                                                  |
| QC - sample gender verification                | SampleGender             | ngs-bits latest      |                                                  |
| QC - low coverage statistics                   | BedLowCoverage           | ngs-bits latest      |                                                  |
| QC - low coverage annotation                   | BedAnnotateGenes         | ngs-bits latest      |                                                  |
| Variant calling - hla                          | HLA-genotyper            | 2022_05              |   https://github.com/axelgschwind/hla-genotyper  |
| Variant calling - SNVs and Indels              | Strelka2                 | 2.9.9                |                                                  |
| Variant calling - b-allele frequency (t/n)     | VariantAnnotateFrequency | ngs-bits latest      |                                                  |
| Variant calling - left-normalization of InDels | VcfLeftNormalize         | ngs-bits latest      |                                                  |
| Annotation                                     | VEP                      | 104.3                |                                                  |

CNV calling and annotation is performed using these tools:

| step                                               | tool                 | version              | comments                                            |
|----------------------------------------------------|----------------------|----------------------|-----------------------------------------------------|
| CNV calling                                        | ClinCNV              | 1.17.2               |                                                     |
| annotation - gene information                      | CnvGeneAnnotation    | ngs-bits latest      |                                                     |
| annotation - overlapping pathogenic CNVs from NGSD | NGSDAnnotateCNV      | ngs-bits latest      |                                                     |

SV calling and annotation is performed using these tools:

| step                                      | tool                            | version              | comments                                            |
|-------------------------------------------|---------------------------------|----------------------|-----------------------------------------------------|
| SV calling                                | Manta                           | 1.6.0                |                                                     |
| annotation - gene information             | BedpeGeneAnnotation             | ngs-bits latest      |                                                     |

Determination of viral load is performed using these tools:

| step                                           | tool                     | version              | comments                                            |
|------------------------------------------------|--------------------------|----------------------|-----------------------------------------------------|
| Filtering reads                                | samtools                 | 1.15.1               |                                                     |
| Mapping                                        | bwa mem2                 | 2.2.1                |                                                     |
| Variant calling                                | freebayes                | 1.3.6                |                                                     |
| Statistics - Coverage                          | BedCoverage              | ngs-bits latest      |                                                     |
| Statistics - BedReadCount                      | BedReadCount             | ngs-bits latest      |                                                     |
| Annotation                                     | BedAnnotateFromBed       | ngs-bits latest      | Several data sources are annotated using this tool. |

Determination of microsatellite instability (MSI) is performed using this tool:

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

**CGI annotation (`ci`):**
This step will send the SNVs and CNVs from `somatic.GSvar` and `somatic_clincnv.tsv` to CancerGenomeInterpreter.org. The results will be downloaded into the out folder and annotated to both files. Provide a tumor type using the parameter `-cancer_type` for more accurate results, acronyms can be found on CancerGenomeInterpreter.org.


[back to the start page](../README.md)
