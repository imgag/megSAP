# megSAP - Combined DNA/RNA Somatic Analysis

The pipeline `somatic_dna_rna.php` implements tumor--normal somatic analysis
for DNA and optionally additional RNA tumor/normal samples.

## Basics

The main parameters are used to specify at least the DNA samples:

* `p_folder` - The project folder containing all of the specified samples.
* `t_dna_id` - Tumor DNA sample identifier.
* `n_dna_id` - Normal DNA sample identifier.

Additional RNA samples are specified with the following parameters:

* `t_rna_id` - Tumor RNA sample identifier.
* `n_rna_id` - Normal RNA sample identifier.

If the pipeline is not used with a database, processing systems need to be
specified manually with parameters `t_dna_sys`, `n_dna_sys`, `t_rna_sys` and
`n_rna_sys`.

## Processing

The pipeline implements several processing steps, which can be specified with
the `process` parameter:

* **dna**: DNA pair analysis, by default somatic variant calling (vc) and
  annotation of somatic variants (an), see `somatic_dna.php`
* **rna**: RNA pair analysis, by default expression differences (rc) and somatic
  fusions (fu), see `somatic_rna.php`
* **germline**: germline analysis of DNA normal sample
* **co**: combines DNA and RNA results
* **igv**: creates IGV session file

*Note: If single DNA and RNA samples are not yet analyzed (e.g. mapping not calculated or no expression data available), add the step `ma` to `steps_dna` or `steps_rna`.*

### Germline Analysis

There are two presets for the germline variant calling (based on the DNA normal
sample):

1. `-germline_preset default`: Runs a default germline analysis, target region
   and naming suffix can be specified optionally with `germline_target` and
   `germline_suffix`.
2. `-germline_preset nearby`: Runs a germline analysis using positions of
   somatic variants and their nearby regions as the target regions, somatic
   variants are annotated with information on the presence of nearby somatic
   variants.

## Output

The following output files are generated in the output directory (`out_folder`):

* somatic variants, as annotated VCF and GSvar files, including RNA observation
* structural variants
* CNVs
* expression differences
* somatic fusions
* germline variants
* IGV session

## Examples

### Project eMed-HCC

For samples which are not yet analyzed:

```
php megSAP/src/Pipelines/somatic_dna_rna.php \
        -p_folder ... \
        -t_dna_id t_dna \
        -n_dna_id n_dna \
        -t_rna_id t_rna \
        -n_rna_id n_rna \
        -steps_dna ma,vc,an \
        -steps_rna ma,rc,fu \
        -filter_set not-coding-splicing \
        -germline_suffix _adme \
        -germline_target .../ADME.bed
```

If single sample analysis are already performed:

```
php megSAP/src/Pipelines/somatic_dna_rna.php \
        -p_folder ... \
        -t_dna_id t_dna \
        -n_dna_id n_dna \
        -t_rna_id t_rna \
        -n_rna_id n_rna \
        -filter_set not-coding-splicing \
        -germline_suffix _adme \
        -germline_target .../ADME.bed
```
### Project IVAC-ALL

For samples which are not yet analyzed:

```
php megSAP/src/Pipelines/somatic_dna_rna.php \
        -p_folder ... \
        -t_dna_id t_dna \
        -n_dna_id n_dna \
        -t_rna_id t_rna \
        -n_rna_id n_rna \
        -steps_dna ma,vc,an \
        -steps_rna ma,rc,fu \
        -filter_set not-coding-splicing,synonymous \
        -germline_preset nearby
```

If single sample analysis are already performed:

```
php megSAP/src/Pipelines/somatic_dna_rna.php \
        -p_folder ... \
        -t_dna_id t_dna \
        -n_dna_id n_dna \
        -t_rna_id t_rna \
        -n_rna_id n_rna \
        -steps_dna ma,vc,an \
        -steps_rna ma,rc,fu \
        -filter_set not-coding-splicing,synonymous \
        -germline_preset nearby
```

[back to the start page](../README.md)
