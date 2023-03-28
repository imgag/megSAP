## Tools used in the somatic analysis:

The following tools are used for calling of small variants and annotation of small variants:

| step                                           | tool                     | version              | comments                                         |
|------------------------------------------------|--------------------------|----------------------|--------------------------------------------------|
| Variant calling - hla                          | HLA-genotyper            | 2022_05              | https://github.com/axelgschwind/hla-genotyper    |
| Variant calling - SNVs and Indels              | Strelka2                 | 2.9.10               | (tumor-normal)                                   |
| Variant calling - SNVs and Indels              | Varscan2                 | 2.4.5                | (tumor-only)                                     |
| Variant calling - left-normalization of InDels | VcfLeftNormalize         | ngs-bits latest      |                                                  |
| Annotation - b-allele frequency (t/n)          | VariantAnnotateFrequency | ngs-bits latest      |                                                  |
| Annotation - general                           | VEP                      | 109.3                |                                                  |


CNV calling and annotation is performed using these tools:

| step                                               | tool                 | version              | comments                                            |
|----------------------------------------------------|----------------------|----------------------|-----------------------------------------------------|
| CNV calling                                        | ClinCNV              | 1.18.3               |                                                     |
| annotation - general                               | BedAnnotateFromBed   | ngs-bits latest      | Several data sources are annotated using this tool. (tumor-only) |
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
| Filtering reads                                | samtools                 | 1.17                 |                                                     |
| Mapping                                        | bwa mem2                 | 2.2.1                |                                                     |
| Variant calling                                | freebayes                | 1.3.6                |                                                     |
| Statistics - Coverage                          | BedCoverage              | ngs-bits latest      |                                                     |
| Statistics - BedReadCount                      | BedReadCount             | ngs-bits latest      |                                                     |
| Annotation                                     | BedAnnotateFromBed       | ngs-bits latest      | Several data sources are annotated using this tool. |

Microsatellite instability (MSI):

| step                                           | tool                     | version              | comments                                            |
|------------------------------------------------|--------------------------|----------------------|-----------------------------------------------------|
| Calling                                        | Mantis                   | 1.0.5                |                                                     |


[to tumor-only analysis page](dna_tumor_only.md)

[to tumor-normal analysis page](dna_tumor-normal_pair.md)

[back to the start page](../README.md)

