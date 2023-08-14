# illumina NovaSeq X Plus

## create custom genome

- can only be done via BaseSpace (see https://knowledge.illumina.com/software/on-instrument-analysis-software/software-on-instrument-analysis-software-troubleshooting-list/000007055)

- create new project on BaseSpace
- upload *.fa file
- open `DRAGEN Reference Builder` in Apps (use version matching the Dragen version on the sequencer)
- select the created project as output
- choose uploaded Fasta file as input
- leave `Liftover Validation` checked and `Mask BED File` at `Autodetect`
- activate the includes for RNA, CNV and HLA (do not include RNA in 4.2.4)
- run application 


## Errors in the documentation from Illumina (SampleSheet_v2)
[illumnia documentation](https://support-docs.illumina.com/SHARE/SampleSheetv2/Content/SHARE/SampleSheetv2/Settings_fNV_mX.htm)

### BCLConvert
- `BarcodeMismatchesIndex1` and `BarcodeMismatchesIndex2` have to be columns in the data section
- Of Index2 the first n bases has to be cut: OverrideCycles: `Y151;I8N2;N2I8;Y151`

### DRAGEN Enrichment/Germline
- Parameter `Bedfile` has to be named `BedFile`
- Parameter `KeepFastQ` has to be named `KeepFastq`