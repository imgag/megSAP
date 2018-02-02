# megSAP - RNA-Seq Analysis

## Introduction

The RNA-seq pipeline implemented in megSAP is tailored for expression
analyses. It starts from reads in FASTQ format and produces spliced read
alignments against the genome, expression values and several other
results.

Installation procedure is outlined in [RNA-seq pipeline
installation](installation_rna_pipeline.md).

## Running an analysis

The pipeline analyzes a single sample and can be run as follows:

```shell
php megSAP/src/Pipelines/analyze_rna.php -folder Sample_A1 -name A1-system sys.ini
```

The main parameters that you have to provide are:

* `folder` - The folder containing FASTQ files and in which all result
  files will be stored.
* `name` - The basename/prefix for all output files.
* `system` - The [processing system INI file](processing_system_ini_file.md).
  It specifies sequencing adapter reads and the genome build.

In addition, you can specify other options including:

* `steps` - Analysis steps to perform, see below.
* `library_type` - The library type (strandedness of reads), possible options
  are ``reverse`` (default), ``forward`` or ``unstranded``.
* `no_splicing` - Disables spliced alignment.

Sample folders are expected to be structured following [bcl2fastq]
default output, i.e. read 1 files have to match `*_R1_001.fastq.gz` and
read 2 files `*_R2_001.fastq.gz` respectively.

[bcl2fastq]: https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html

## Output and Analysis Details

RNA-seq reads are processed by the following steps:

* [Read quality control and adapter trimming](#read-quality-control-and-adapter-trimming)
  - preprocesses FASTQ reads for subsequent mapping
* [Mapping](#mapping)
  - performs alignment against genome (``-steps`` argument ``ma``)
* [Read counting](#read-counting)
  - summarizes and counts aligned reads for genomic features (``-steps``
	  arguments ``rc`` and ``an``)
* [Fusion detection](#fusion-detection)
  - identifies fusion genes (``-steps`` argument ``fu``)
* [Transcript quantification](#transcript-quantification)
  - produces transcript-level expression values (``-steps`` argument ``tx``)

All steps produce associated quality control values to ensure proper
sample and analysis quality and to facilitate troubleshooting.

### Read quality control and adapter trimming

Quality control and adapter removal is performed by [SeqPurge] ([Sturm
et al. 2016]) for paired-end reads, and by [ReadQC] and [skewer] ([Jiang
et al. 2014]) for single-end reads. Both tools are configured to remove
adapters only, reads are not trimmed by base qualities (for further
details see [Williams et al. 2016]).

[SeqPurge]: https://github.com/imgag/ngs-bits
[ReadQC]: https://github.com/imgag/ngs-bits
[skewer]: https://github.com/relipmoc/skewer
[Sturm et al. 2016]: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1069-7
[Jiang et al. 2014]: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-182
[Williams et al. 2016]: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0956-2

**Output files:**

* ``Prefix_stats_fastq.qcML`` - SeqPurge (or ReadQC) store quality
  control values, including read count, number of bases, quality scores,
  etc. in [qcML format](https://www.ncbi.nlm.nih.gov/pubmed/24760958).

Trimmed FASTQ are stored as temporary files only.

### Mapping

Read alignments are computed by
[STAR](https://github.com/alexdobin/STAR). On request,
[abra2](https://github.com/mozack/abra2) is used to perform an indel
realignment. [MappingQC](https://github.com/imgag/ngs-bits) is run on
the final BAM file for quality control and statistics.

**Output files:**

* ``Prefix.bam`` - alignment in BAM format.
* ``Prefix_splicing.tsv`` - TODO
* ``Prefix_chimeric.tsv`` - TODO
* ``Prefix_stats_map.qcMl`` - statistics on mapping results in qcML
  format.

### Read counting

Mapped read summarization on gene-level is performed by
[featureCounts](http://bioinf.wehi.edu.au/featureCounts/).  Reads are
normalized by counts per million mapped reads method (cpm), fragments
per kilobase of exon per million reads mapped method (fpkm) and
transcripts per million method (tpm). Annotation of read counts with
gene symbols and gene biotypes is based on the provided GTF file.

**Output files:**

* ``Prefix_counts_raw.tsv`` - raw read counts in featureCounts native
  output format (tab-separated)
* ``Prefix_counts.tsv`` - normalized and annotated read counts
  (tab-separated).
* ``Prefix_stats_rc.tsv`` - read assignment and expression quality
  control values.

### Fusion detection

[STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion) is used to
process the STAR alignment output to identify gene fusions.

**Output files:**

* ``Prefix_fusions.tsv``

### Transcript quantification

Transcript quantification is calculated by [salmon] and does not require
the primary STAR output as it directly uses the pre-processed FASTQ
files.

[salmon]: https://github.com/COMBINE-lab/salmon

**Output files:**

* ``Prefix_transcript_quant.tsv``

## Quality Control

## Downstream Analysis

---

## FAQ

### Using other genomes

To use the RNA expression pipeline with custom genomes, you need to
provide:

* the genome FASTA file, e.g. `megSAP/data/genomes/CustomGenome.fa`, and
* the STAR genome index, e.g. `megSAP/data/genomes/STAR/CustomGenome/`
* the gene annotation file in Ensembl-like GTF format, e.g.
  `megSAP/data/dbs/gene_annotations/CustomGenome.gtf` (for read counting
  and annotation)
* the STAR-Fusion database, e.g.
  `megSAP/data/genomes/STAR-Fusion/CustomGenome/` (for gene fusion
  detection)
* the salmon cDNA database, e.g.
  `megSAP/data/genomes/salmon/CustomGenome_cDNA` (for transcript
  quantification)