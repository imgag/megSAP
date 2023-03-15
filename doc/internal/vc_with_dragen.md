# Variant calling with Dragen

Small variant calling and structural variant calling with Dragen is possible now.

## Small variants comparison

|WES (NA12878_58)                                                     |freebayes               |dragen v3.8.4           |dragen v4.0.3           |
|---------------------------------------------------------------------|------------------------|------------------------|------------------------|
|all                                                                  |128106                  |140852                  |140724                  |
|AF<1%                                                                |5573                    |6277                    |6179                    |
|AF<1%, NGSD<20                                                       |1400                    |2919                    |2857                    |
|AF<1%, NGSD<20, high/moderate/low impact                             |356                     |410                     |392                     |
|recessive stringent                                                  |38                      |62                      |55                      |
|recessive relaxed                                                    |69                      |111                     |99                      |
|dominant stringent                                                   |10                      |9                       |9                       |
|dominant relaxed                                                     |162                     |185                     |174                     |
|dominant relaxed + pheno                                             |72                      |198                     |197                     |
|rank1:                                                               |chr11:6617154 C>T (8.00)|chr11:6617154 C>T (8.00)|chr11:6617154 C>T (8.00)|
|                                                                     |                        |                        |                        |
|overlap with freebayes calls (Agilent v7 ROI)                        |                        |97.05% + 1.27% GT diff  |96.94% + 1.26% GT diff  |
|overlap with freebayes calls  (Agilent v7 ROI, only high conf region)|                        |99.31% + 0.18% GT diff  |99.27% + 0.17%  GT diff |


|WGS (NA12878_45)                                                     |freebayes               |dragen v4.0.3           |
|---------------------------------------------------------------------|------------------------|------------------------|
|all                                                                  |365620                  |373779                  |
|AF<1%                                                                |278628                  |284172                  |
|AF<1%, NGSD<20                                                       |43214                   |92845                   |
|AF<1%, NGSD<20, high/moderate/low impact                             |406                     |478                     |
|recessive stringent                                                  |46                      |100                     |
|recessive relaxed                                                    |89                      |167                     |
|dominant stringent                                                   |26                      |100                     |
|dominant relaxed                                                     |189                     |229                     |
|dominant relaxed + pheno                                             |1140                    |1896                    |
|rank1:                                                               |chr11:6617154 C>T (8.50)|chr11:6617154 C>T (8.50)|
|                                                                     |                        |                        |
|overlap with freebayes calls (Agilent v7 ROI)                        |                        |95.21% + 1.81%  GT diff |
|overlap with freebayes calls  (Agilent v7 ROI, only high conf region)|                        |98.05% + 0.91%  GT diff |

**Conclusion: Small variant calls are pretty similar, except for artefacts in low complexity regions.**

## Small variants benchmark

Bechmarks on GIAB sample NA12878 based on the Agilent V7 exome target region:

|SNVs                                                                 |sensitivity             |ppv                     |genotyping|
|---------------------------------------------------------------------|------------------------|------------------------|----------|
|WES default                                                          |99,87                   |98,92                   |99,96     |
|WES dragen v4.0.3                                                    |99,89                   |99,52                   |99,98     |
|WGS default                                                          |99,95                   |98,22                   |99,60     |
|WGS dragen v4.0.3                                                    |99,91                   |99,19                   |99,95     |

|InDels                                                               |sensitivity             |ppv                     |genotyping|
|---------------------------------------------------------------------|------------------------|------------------------|----------|
|WES default                                                          |95,29                   |93,72                   |96,21     |
|WES dragen v4.0.3                                                    |97,64                   |95,93                   |98,49     |
|WGS default                                                          |98,24                   |97,85                   |98,67     |
|WGS dragen v4.0.3                                                    |99,74                   |99,15                   |99,79     |


## Structural variants comparison

Structural variant calling with Dragen is performed only if `dragen_sv_calling` is `true` in the settings.

|WES (NA12878_58)                                                      |freebayes               |dragen v3.8.4           |dragen v4.0.3 |
|----------------------------------------------------------------------|------------------------|------------------------|--------------|
|all                                                                   |201                     |198                     |238           |
|nospecial, filter PASS                                                |82                      |75                      |86            |
|nospecial, filter PASS, qual>=100                                     |72                      |72                      |86            |
|nospecial, filter PASS, qual>=100, PE read depth>=5                   |50                      |46                      |44            |
|nospecial, filter PASS, qual>=100, PE read depth>=5, NGSD<=30         |5                       |27                      |27            |
|nospecial, filter PASS, qual>=100, PE read depth>=5, NGSD<=30, bpd<=20|0                       |7                       |4             |
|exome/wgs stringent                                                   |5                       |27                      |27            |
|exome/wgs relaxed                                                     |15                      |97                      |129           |
|                                                                      |                        |                        |              |
|SVs MANTA (qual>100):                                                 |                        |146                     |146           |
|SVs DRAGEN (qual>100):                                                |                        |172                     |201           |
|overlap (pos):                                                        |                        |93                      |92            |
|overlap (pos+type):                                                   |                        |(type diff) 63          |(type diff) 62|
|overlap (pos+type+GT):                                                |                        |63                      |60            |

|WGS (NA12878_45)                                                      |freebayes               |dragen v4.0.3           |
|----------------------------------------------------------------------|------------------------|------------------------|
|all                                                                   |10992                   |13699                   |
|nospecial, filter PASS                                                |8236                    |12298                   |
|nospecial, filter PASS, qual>=100                                     |7734                    |12074                   |
|nospecial, filter PASS, qual>=100, PE read depth>=5                   |3725                    |4908                    |
|nospecial, filter PASS, qual>=100, PE read depth>=5, NGSD<=30         |65                      |579                     |
|nospecial, filter PASS, qual>=100, PE read depth>=5, NGSD<=30, bpd<=20|43                      |151                     |
|exome/wgs stringent                                                   |116                     |621                     |
|exome/wgs relaxed                                                     |765                     |3272                    |
|                                                                      |                        |                        |
|SVs MANTA (qual>100):                                                 |                        |8585                    |
|SVs DRAGEN (qual>100):                                                |                        |12643                   |
|overlap (pos):                                                        |                        |7197                    |
|overlap (pos+type):                                                   |                        |6941                    |
|overlap (pos+type+GT):                                                |                        |6797                    |


**Conclusion: Structural variant calls differ quite a bit. Dragen calls more variants, but is more conservative when not sure about the type (reports BND instead of DEL/DUP from Manta). It also reports small tandem duplications (< 1000 bases) as insertions!**

## Internal documentation

GutHub issue: https://github.com/imgag/megSAP/issues/132  
Benchmarks were performed in the folder: /mnt/storage3/users/ahsturm1/Sandbox/2022\_09\_07\_dragen\_vc/

