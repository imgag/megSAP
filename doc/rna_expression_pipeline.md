# RNA-seq Expression Pipeline

This document details all steps of the RNA-seq expression pipeline for paired-end reads (script `analyze_rna.php`).

### Pre-processing

Quality control and adapter removal is performed using [SeqPurge](https://github.com/imgag/ngs-bits) 
([Sturm et al. 2016](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1069-7)).

No read trimming using base quality is performed, for further details see
[Williams et al. 2016](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0956-2).

_Output: read QC data in qcML format_

### Mapping (step `ma`)

Read alignments are computed by [STAR](https://github.com/alexdobin/STAR). On
request. [MappingQC](https://github.com/imgag/ngs-bits) is run on the final
BAM file for quality control and statistics.

_Output: alignment in BAM format, alignment QC data in qcML format_

### Read counting (steps `rc`, `an`)

Mapped read summarization is performed by
[featureCounts](http://bioinf.wehi.edu.au/featureCounts/). Normalization of read
counts by counts per million mapped reads (cpm) method or by fragments per
kilobase of exon per million reads mapped (fpkm) method is calculated by a
script (`rc_normalize.php`). Annotation of read counts based on the provided GTF
file is implemented in `rc_annotate.php`.

_Output: raw read counts in tabular featureCounts format, normalized (and
annotated) read counts in tabular format_

### Fusion detection (step `fu`)

[Arriba](https://github.com/suhrig/arriba) can be used to to identify fusion transcripts.

_Output: candidate fusions in tabular format_
