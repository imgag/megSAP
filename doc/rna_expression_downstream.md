# RNA-seq Expression -- Downstream Analysis

### Expression matrix

The following R code snippets help to create a single expression matrix from
multiple raw read count files as produced by the
[RNA-seq Expression Pipeline](rna_expression.md).

```r
library(tools)
library(data.table)

# create unique sample identifier from filename
extract.sampleid <- function(x) file_path_sans_ext(basename(x))

# read a single-sample featureCounts output file
read.featurecounts <- function(f) {
  df <- fread(f, sep="\t", skip=1, header=TRUE, drop=2:6, data.table=FALSE)
  row.names(df) <- df$Geneid
  df <- df[,-1, drop=FALSE]
  df <- df[,sort(colnames(df)), drop=FALSE]
  colnames(df) <- extract.sampleid(colnames(df))
  return(df)
}

# read multiple single-sample featureCounts output files and construct expression matrix
read.many.featurecounts <- function(files) {
  all.tables <- sapply(files,
                       read.featurecounts,
                       simplify=FALSE, USE.NAMES=TRUE)
  genes.ref <- row.names(all.tables[[1]])
  check.ids <- all(sapply(all.tables, function(t) all(row.names(t) == genes.ref)))
  stopifnot(check.ids)

  df <- Reduce(cbind, all.tables)

  return(df)
}
```

Using the default project/sample directory structure, one way to create the
expression matrix is then simply performed by:
```r
fs <- list.files(pattern = "_counts_raw.tsv$", full.names = TRUE, recursive = TRUE)
counts <- read.many.featurecounts(fs)
head(counts)
```

```
                SampleA SampleB SampleC
ENSG00000223972       0       0       0
ENSG00000227232       0       0       0
ENSG00000243485       0       0       0
ENSG00000237613       0       0       0
ENSG00000268020       0       0       0
ENSG00000240361       0       0       0
```
