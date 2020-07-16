### Comparison megSAP pipline using DRAGEN mapping vs. bwa mapping

All runs are performed on NA12878 with the standard parameters of megSAP and compared to the GIAB reference for NA12878. This means for all analysis an indel realignment with Abra2 was done. 

The NA12878 WGS sample was sequenced with the illumina TruSeq DNA PCR-Free kit and the NA12878 WES sample was sequences using the Agilent SureSelectXT Human All Exon V7 kit.

| Type <br> | Mapping <br> | <br> recall/sensitivity | SNVS <br> precision/ppv | <br> genotyping_accuracy | <br> recall/sensitivity | INDEL <br> precision/ppv | <br> genotyping_accuracy |
| --------- | ------------ | ----------------------: | ----------------------: | -----------------------: | ----------------------: | -----------------------: | -----------------------: |
| WGS       | DRAGEN       |                  99.96% |                  99.74% |                   99.98% |                  97.71% |                   99.52% |                   98.36% |
| WGS       | BWA          |                  99.96% |                  99.72% |                   99.98% |                  97.73% |                   99.52% |                   98.34% |
|           |              |                         |                         |                          |                         |                          |                          |
| WES       | DRAGEN       |                  99.61% |                  99.59% |                   99.97% |                  93.23% |                   94.96% |                   96.45% |
| WES       | BWA          |                  99.62% |                  99.53% |                   99.97% |                  93.13% |                   95.05% |                   96.50% |

