### Comparison megSAP pipline using DRAGEN mapping vs. bwa mapping

All runs are performed with the standard parameters of megSAP (e.g. with indel realignment using ABRA2) and compared to the GIAB reference for NA12878.

| Sample <br> | Type <br> | Mapping <br> | <br> recall/sensitivity | Combined <br> precision/ppv | <br> genotyping_accuracy | <br> recall/sensitivity | SNVS <br> precision/ppv | <br> genotyping_accuracy | <br> recall/sensitivity | INDEL <br> precision/ppv | <br> genotyping_accuracy |
| ----------- | --------- | ------------ | ----------------------: | --------------------------: | -----------------------: | ----------------------: | ----------------------: | -----------------------: | ----------------------: | -----------------------: | -----------------------: |
| NA12878_45  | WGS       | DRAGEN       |                  99.65% |                      99.71% |                   99.76% |                  99.96% |                  99.74% |                   99.98% |                  97.71% |                   99.52% |                   98.36% |
| NA12878_45  | WGS       | BWA          |                  99.65% |                      99.70% |                   99.76% |                  99.96% |                  99.72% |                   99.98% |                  97.73% |                   99.52% |                   98.34% |
|             |           |              |                         |                             |                          |                         |                         |                          |                         |                          |                          |
| NA12878_58  | WES       | DRAGEN       |                  99.17% |                      99.28% |                   99.74% |                  99.61% |                  99.59% |                   99.97% |                  93.23% |                   94.96% |                   96.45% |
| NA12878_58  | WES       | BWA          |                  99.18% |                      99.22% |                   99.75% |                  99.62% |                  99.53% |                   99.97% |                  93.13% |                   95.05% |                   96.50% |

