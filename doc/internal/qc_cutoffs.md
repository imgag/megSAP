# Sample QC cutoffs


### lrGS germline

| QC value                 | good           | medium          | bad             |
|--------------------------|----------------|-----------------|-----------------|
| target region read depth | > 30x          | 25-30x          | < 25x           |
| 20x coverage             | > 95%          | 85-95%          | < 85%           |
| SNV deviation            | < 4%           | 4-6%            | > 6%            |
| N50                      | >=10kb         |  < 10kb         |                 |

see <https://vswliscurator01.ukt.ad.local/web/0/intra/index.php?art_id=dc_2023_11_15_fc131c46c46463ef86>

### WGS germline

| QC value                 | good           | medium           | bad             |
|--------------------------|----------------|------------------|-----------------|
| target region read depth | > 30x          | 25-30x           | < 25x           |
| 20x coverage             | > 98%          | 95-98%           | < 95%           |
| SNV deviation            | < 3%           | 3-6%             | > 6%            |
| CNV count                | <=2500         | > 2500           |                 |

see <https://vswliscurator01.ukt.ad.local/web/0/intra/index.php?art_id=dc_2021_09_24_c4e4ec2d1b5b118a21>

### WES germline

| QC value                 | good           | medium          | bad             |
|--------------------------|----------------|-----------------|-----------------|
| target region read depth | > 100x         | 80-100x         | < 80x           |
| 20x coverage             | > 95%          | 90-95%          | < 90%           |
| SNV deviation            | < 3%           | 3-6%            | > 6%            |

see <https://vswliscurator01.ukt.ad.local/web/0/intra/index.php?art_id=dc_2021_09_24_c4e4ec2d1b5b118a21>

### WES tumor

| QC value                 | good           | medium          | bad             |
|--------------------------|----------------|-----------------|-----------------|
| target region read depth | > 200x         | 150-200x        | < 150x          |
| 20x coverage             | > 95%          | 95-90%          | < 90%           |

see <https://vswliscurator01.ukt.ad.local/web/0/intra/index.php?art_id=dc_2021_09_24_adb98ae06465357102>


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
| 20x coverage             | > 99%          | 89-99%          | <98%            |
| target region read depth | > 1000x        | > 1000x         |                 |
