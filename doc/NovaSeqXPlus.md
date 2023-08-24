# illumina NovaSeq X Plus

## create custom genome

- can only be done via BaseSpace (see https://knowledge.illumina.com/software/on-instrument-analysis-software/software-on-instrument-analysis-software-troubleshooting-list/000007055)

- create new project on BaseSpace
- upload *.fa file
- open `Reference Builder (Instruments) v1.1.0` in Apps 
- select the created project as output
- select the correct Sequencer (`NovaSeq X Series`)
- enter name. This has to be exactly the name the genome should be called! It cannot be changed afterwards.
- enter organization and species
- give a link/description for the source of the *.fa file
(- set `Mask BED File` to `hg38` (worked once for me))
- leave `SAM Liftover File`at `None`
(- GTF Annotation and Methylation was left empty/unchecked)
- run application 
(the job crashes since 2 weeks, simply retry it several times)


## Uploading enrichment files (BED) to the sequencer

- In the upload dialog set the file type to `BedFile`(not `bed`), otherwise it will not be found in the sample sheet 


## Re-analyse samples

- go to completed runs
- select the run and click on requeue secondary analysis
- upload SampleSheet
 

## Errors in the documentation from Illumina (SampleSheet_v2)
[illumnia documentation](https://support-docs.illumina.com/SHARE/SampleSheetv2/Content/SHARE/SampleSheetv2/Settings_fNV_mX.htm)

### BCLConvert
- `BarcodeMismatchesIndex1` and `BarcodeMismatchesIndex2` have to be columns in the data section
- Of Index2 the first n bases has to be cut: OverrideCycles: `Y151;I8N2;N2I8;Y151`

### DRAGEN Enrichment/Germline
- Parameter `Bedfile` has to be named `BedFile`
- Parameter `KeepFastQ` has to be named `KeepFastq`
