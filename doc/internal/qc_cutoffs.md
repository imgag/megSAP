# Sample QC cutoffs


### lrGS

| QC value                 | good           | medium          | bad             |
|--------------------------|----------------|-----------------|-----------------|
| target region read depth | > 30x          | > 25x           | < 25x           |
| 20x coverage             | > 95%          | > 85%           | < 85%           |
| SNV deviation            | < 5            |                 |                 |
| N50                      | > 10kb         |  < 10kb         |                 |


### WGS

| QC value                 | good           | medium           | bad             |
|--------------------------|----------------|------------------|-----------------|
| target region read depth | > 33x          | > 30x            | < 30x           |
| 20x coverage             | > 98.5%        | 90-98.5%         | < 90%           |
| SNV deviation            | < 3%           | 3-6%             | > 6%            |
| CNV count                | 1000-2500      | < 1000 \| > 2500 |                |
| AT dropout               | < 3%           | < 8%             | > 8%           |

- high SNV deviation (>5) indicates contanmination!


### WES normal

| QC value                 | good           | medium          | bad             |
|--------------------------|----------------|-----------------|-----------------|
| target region read depth | > 100x         | > 80x           | < 80x           |
| 20x coverage             | > 95%          | > 90%           | < 90%           |
| SNV deviation            | < 4%           | > 4%            | > 5%            |

### WES tumor

| QC value                 | good           | medium          | bad             |
|--------------------------|----------------|-----------------|-----------------|
| target region read depth | > 200x         | > 150x          | < 150x          |
| 20x coverage             | > 95%          | > 90%           | < 90%           |


### RNA

| QC value                         | good           | medium          | bad             |
|----------------------------------|----------------|-----------------|-----------------|
| target region read depth         | > 8x           | > 8x            | < 8x            |
| house keeping genes 10x coverage | > 50%          | > 27%           | < 27%           |


### cfDNA

| QC value                 | good           | medium          | bad             |
|--------------------------|----------------|-----------------|-----------------|
| target region read depth | > 3000x        | > 1000x         | < 1000x         |
| target read depth 2-fold | > 1000x        | > 500x          | < 500x          |
| monitoring 250x coverage | ~ 100%         |                 |                 |
| cfDNA-tumor correlation  | > 95%          |  > 80%          |                 |

- if cfDNA-tumor correlation is below 95%: check correlation to other cfDNA/normal sample 


### RPGR-Ex15

| QC value                 | good           | medium          | bad             |
|--------------------------|----------------|-----------------|-----------------|
| 20x coverage             | ~ 99%          |                 |                 |
| target region read depth | > 20 000x      |                 |                 |


