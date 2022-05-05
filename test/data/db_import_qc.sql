
INSERT INTO processing_system (name_short, name_manufacturer, shotgun, genome_id) VALUES ('ssX', 'SureSelectXT X-Chromosome', '1', 1);

INSERT INTO device (type, name) VALUES ('GAIIx', 'test');

INSERT INTO sequencing_run (name, fcid, start_date, end_date, device_id, recipe) VALUES ('#002', 'FCID4711', '2013-02-04', '2013-02-04', 1, '76+6+76');

INSERT INTO sender (name) VALUES ('John Doe');

INSERT INTO project (name, type, internal_coordinator_id) VALUES ('X-Chr','research', 1);

INSERT INTO sample (name, sample_type, species_id, gender, tumor, ffpe, sender_id) VALUES ('GS120676', 'DNA', 1, 'male', '0', '0', 1);
INSERT INTO processed_sample (sample_id, process_id, sequencing_run_id, lane, processing_system_id, project_id) VALUES (1,1,1,'1',1,1);

INSERT INTO sample (name, sample_type, species_id, gender, tumor, ffpe, sender_id) VALUES ('GS120677', 'DNA', 1, 'male', '0', '0', 1);
INSERT INTO processed_sample (sample_id, process_id, sequencing_run_id, lane, processing_system_id, project_id) VALUES (2,1,1,'1',1,1);

INSERT INTO sample (name, sample_type, species_id, gender, tumor, ffpe, sender_id) VALUES ('GS120678', 'DNA', 1, 'male', '0', '0', 1);
INSERT INTO processed_sample (sample_id, process_id, sequencing_run_id, lane, processing_system_id, project_id) VALUES (3,1,1,'1',1,1);

INSERT INTO sample (name, sample_type, species_id, gender, tumor, ffpe, sender_id) VALUES ('DXcf120678', 'cfDNA', 1, 'male', '0', '0', 1);
INSERT INTO processed_sample (sample_id, process_id, sequencing_run_id, lane, processing_system_id, project_id) VALUES (4,1,1,'1',1,1);

INSERT INTO sample (name, sample_type, species_id, gender, tumor, ffpe, sender_id) VALUES ('RX120678', 'RNA', 1, 'male', '0', '0', 1);
INSERT INTO processed_sample (sample_id, process_id, sequencing_run_id, lane, processing_system_id, project_id) VALUES (5,1,1,'1',1,1);

INSERT INTO `qc_terms` (`id`, `qcml_id`, `name`, `description`, `type`, `obsolete`) VALUES
(1, 'QC:2000005', 'read count', 'Total number of reads (one cluster in a paired-end experiment generates two reads).', 'int', 0),
(2, 'QC:2000006', 'read length', 'Raw read length of a single read before trimming. Comma-separated list of lengths if several.', 'string', 0),
(3, 'QC:2000007', 'Q20 read percentage', 'The percentage of reads with a mean base quality score greater than Q20.', 'float', 0),
(4, 'QC:2000008', 'Q30 base percentage', 'The percentage of bases with a minimum quality score of Q30.', 'float', 0),
(5, 'QC:2000009', 'no base call percentage', 'The percentage of bases without base call (N).', 'float', 0),
(6, 'QC:2000010', 'gc content percentage', 'The percentage of bases that are called to be G or C.', 'float', 0),
(7, 'QC:2000013', 'variant count', 'Total number of variants in the target region.', 'int', 0),
(8, 'QC:2000014', 'known variants percentage', 'Percentage of variants that are known polymorphisms in the dbSNP database.', 'float', 0),
(9, 'QC:2000015', 'high-impact variants percentage', 'Percentage of variants with high impact on the protein, i.e. stop-gain, stop-loss, frameshift, splice-acceptor or splice-donor variants.', 'float', 0),
(10, 'QC:2000016', 'homozygous variants percentage', 'Percentage of variants that are called as homozygous.', 'float', 0),
(11, 'QC:2000017', 'indel variants percentage', 'Percentage of variants that are insertions/deletions.', 'float', 0),
(12, 'QC:2000018', 'transition/transversion ratio', 'Transition/transversion ratio of single nucleotide variants.', 'float', 0),
(13, 'QC:2000019', 'trimmed base percentage', 'Percentage of bases that were trimmed during to adapter or quality trimming.', 'float', 0),
(14, 'QC:2000020', 'mapped read percentage', 'Percentage of reads that could be mapped to the reference genome.', 'float', 0),
(15, 'QC:2000021', 'on-target read percentage', 'Percentage of reads that could be mapped to the target region.', 'float', 0),
(16, 'QC:2000022', 'properly-paired read percentage', 'Percentage of properly paired reads (for paired-end reads only).', 'float', 0),
(17, 'QC:2000023', 'insert size', 'Average insert size (for paired-end reads only).', 'float', 0),
(18, 'QC:2000024', 'duplicate read percentage', 'Percentage of reads removed because they were duplicates (PCR, optical, etc).', 'float', 0),
(19, 'QC:2000025', 'target region read depth', 'Average sequencing depth in target region.', 'float', 0),
(20, 'QC:2000026', 'target region 10x percentage', 'Percentage of the target region that is covered at least 10-fold.', 'float', 0),
(21, 'QC:2000027', 'target region 20x percentage', 'Percentage of the target region that is covered at least 20-fold.', 'float', 0),
(22, 'QC:2000028', 'target region 30x percentage', 'Percentage of the target region that is covered at least 30-fold.', 'float', 0),
(23, 'QC:2000029', 'target region 50x percentage', 'Percentage of the target region that is covered at least 50-fold.', 'float', 0),
(24, 'QC:2000030', 'target region 100x percentage', 'Percentage of the target region that is covered at least 100-fold.', 'float', 0),
(25, 'QC:2000031', 'target region 200x percentage', 'Percentage of the target region that is covered at least 200-fold.', 'float', 0),
(26, 'QC:2000032', 'target region 500x percentage', 'Percentage of the target region that is covered at least 500-fold.', 'float', 0),
(27, 'QC:2000033', 'error estimation read depth', 'Average read depth on the special target region used for error estimation after mapping.', 'float', 1),
(28, 'QC:2000034', 'error estimation N percentage', 'No base call (N) percentage determined on special target region after mapping.', 'float', 1),
(29, 'QC:2000035', 'error estimation SNV percentage', 'SNV error percentage determined on special target region after mapping.', 'float', 1),
(30, 'QC:2000036', 'error estimation indel percentage', 'indel error percentage determined on special target region after mapping.', 'float', 1),
(31, 'QC:2000039', 'gender check', 'Gender check result: \'n/a\', \'male_passed\', \'male_failed\', \'female_passed\', \'female_failed\'.', 'string', 0),
(32, 'QC:2000040', 'sample correlation', 'Correlation of high-quality SNV genotypes between tumor and normal sample.', 'float', 0),
(33, 'QC:2000041', 'somatic variant count', 'Total number of somatic variants in the target region.', 'int', 0),
(34, 'QC:2000042', 'somatic indel variants percentage', 'Percentage of somatic variants that are insertions/deletions.', 'float', 0),
(35, 'QC:2000043', 'somatic transition/transversion ratio', 'Transition/transversion ratio of somatic single nucleotide variants.', 'float', 0),
(36, 'QC:2000044', 'somatic CNVs count', 'Count of large somatic CNVs.', 'float', 0),
(37, 'QC:2000045', 'known somatic variants percentage', 'Percentage of somatic variants that are listed as germline variants in public datbases (e.g. AF>1% in ExAC).', 'float', 0),
(38, 'QC:2000049', 'bases sequenced (MB)', 'Bases sequenced in total (in megabases).', 'float', 0),
(39, 'QC:2000050', 'bases usable (MB)', 'Bases sequenced that are usable for variant calling (in megabases).', 'float', 0),
(40, 'QC:2000051', 'SNV allele frequency deviation', 'Percentage of common SNPs that deviate from the expected allele frequency (i.e. 0.0, 0.5 or 1.0 for diploid organisms).', 'float', 0),
(41, 'QC:2000052', 'clipped base percentage', 'Percentage of the bases that are soft-clipped or hand-clipped during mapping.', 'float', 0),
(42, 'QC:2000053', 'somatic variant rate', 'Categorized somatic variant rate followed by the somatic variant rate [variants/Mbp] normalized for the target region and exome size and corrected for tumor suppressors.', 'float', 0),
(43, 'QC:2000054', 'tumor content estimate', 'Estimate of tumor content.', 'float', 0),
(44, 'QC:2000057', 'near-target read percentage', 'Percentage of reads that were mapped to the target region or near the target region (at most 250 bases away).', 'float', 0),
(45, 'QC:2000058', 'target region half depth percentage', 'Percentage of the target region that is covered at least with half of the target region average depth. This is a measure of coverage uniformity.', 'float', 0),
(46, 'QC:2000059', 'AT dropout', 'Illumina-style AT dropout metric. Calculated by taking each GC bin independently and calculating (%ref_at_gc - %reads_at_gc) and summing all positive values for GC=[0..50].', 'float', 0),
(47, 'QC:2000060', 'GC dropout', 'Illumina-style GC dropout metric. Calculated by taking each GC bin independently and calculating (%ref_at_gc - %reads_at_gc) and summing all positive values for GC=[50..100].', 'float', 0),
(48, 'QC:2000065', 'target region 1000x percentage', 'Percentage of the target region that is covered at least 1000-fold.', 'float', 0),
(49, 'QC:2000066', 'target region 2500x percentage', 'Percentage of the target region that is covered at least 2500-fold.', 'float', 0),
(50, 'QC:2000067', 'target region 5000x percentage', 'Percentage of the target region that is covered at least 5000-fold.', 'float', 0),
(51, 'QC:2000068', 'target region 7500x percentage', 'Percentage of the target region that is covered at least 7500-fold.', 'float', 0),
(52, 'QC:2000069', 'target region 10000x percentage', 'Percentage of the target region that is covered at least 10000-fold.', 'float', 0),
(53, 'QC:2000070', 'target region 15000x percentage', 'Percentage of the target region that is covered at least 15000-fold.', 'float', 0),
(54, 'QC:2000071', 'target region read depth 2-fold duplication', 'Average coverage with at least 2-fold duplication', 'float', 0),
(55, 'QC:2000072', 'target region read depth 3-fold duplication', 'Average coverage with at least 3-fold duplication', 'float', 0),
(56, 'QC:2000073', 'target region read depth 4-fold duplication', 'Average coverage with at least 4-fold duplication', 'float', 0),
(57, 'QC:2000074', 'raw target region read depth', 'Raw average sequencing depth in target region.', 'float', 0),
(58, 'QC:2000077', 'monitoring variant read depth', 'Average sequencing depth of monitoring variants.', 'float', 0),
(59, 'QC:2000078', 'ID variant read depth', 'Average sequencing depth of ID variants.', 'float', 0),
(60, 'QC:2000079', 'monitoring variant count', 'Total number of monitoring variants.', 'int', 0),
(61, 'QC:2000080', 'monitoring variant 250x percentage', 'Percentage of monitoring variants which have at least 250x coverage.', 'int', 0),
(62, 'QC:2000081', 'ID variant count', 'Total number of ID variants.', 'int', 0),
(63, 'QC:2000082', 'ID variant 250x percentage', 'Percentage of ID variants which have at least 250x coverage.', 'int', 0),
(64, 'QC:2000083', 'cfDNA-tumor correlation', 'Sample correlation of high-quality SNV genotypes between cfDNA and tumor sample.', 'float', 0),
(65, 'QC:2000084', 'cfDNA-cfDNA correlation', 'Sample correlation of high-quality SNV genotypes between related cfDNA samples.', 'string', 0),
(66, 'QC:2000085', 'umiVar error rate 1-fold duplication', 'umiVar error rate for 1-fold duplicated bases.', 'float', 0),
(67, 'QC:2000086', 'umiVar error rate 2-fold duplication', 'umiVar error rate for 2-fold duplicated bases.', 'float', 0),
(68, 'QC:2000087', 'umiVar error rate 3-fold duplication', 'umiVar error rate for 3-fold duplicated bases.', 'float', 0),
(69, 'QC:2000088', 'umiVar error rate 4-fold duplication', 'umiVar error rate for 4-fold duplicated bases.', 'float', 0),
(70, 'QC:2000089', 'raw somatic variant rate', 'Somatic variant rate in variants per Megabase without normalization to TSG/Oncogenes or exome size. SNVs in blacklisted genes were discarded for the calculation.', 'float', 0),
(71, 'QC:2000090', 'somatic custom target 10x percentage', 'Percentage of the somatic custom panel that is covered at least 10-fold.', 'float', 0),
(72, 'QC:2000091', 'somatic custom target 20x percentage', 'Percentage of the somatic custom panel that is covered at least 20-fold.', 'float', 0),
(73, 'QC:2000092', 'somatic custom target 30x percentage', 'Percentage of the somatic custom panel that is covered at least 30-fold.', 'float', 0),
(74, 'QC:2000093', 'somatic custom target 50x percentage', 'Percentage of the somatic custom panel that is covered at least 50-fold.', 'float', 0),
(75, 'QC:2000094', 'somatic custom target 100x percentage', 'Percentage of the somatic custom panel that is covered at least 100-fold.', 'float', 0),
(76, 'QC:2000095', 'somatic custom target 200x percentage', 'Percentage of the somatic custom panel that is covered at least 200-fold.', 'float', 0),
(77, 'QC:2000096', 'somatic custom target 500x percentage', 'Percentage of the somatic custom panel that is covered at least 500-fold.', 'float', 0),
(78, 'QC:2000097', 'somatic custom target region read depth', 'Average sequencing read depth in somatic custom panel region.', 'float', 0),
(79, 'QC:2000098', 'somatic custom target 60x percentage', 'Percentage of the somatic custom panel that is covered at least 60-fold.', 'float', 0),
(80, 'QC:2000099', 'target region 60x percentage', 'Percentage of the target region that is covered at least 60-fold.', 'float', 0),
(81, 'QC:2000100', 'housekeeping genes read percentage', 'Percentage of reads that could be mapped to the exon region of housekeeping genes.', 'float', 0),
(82, 'QC:2000101', 'housekeeping genes read depth', 'Average sequencing depth in exon region of housekeeping genes.', 'float', 0),
(83, 'QC:2000102', 'housekeeping genes 10x percentage', 'Percentage of the exon region of housekeeping genes that is covered at least 10-fold.', 'float', 0),
(84, 'QC:2000103', 'housekeeping genes 20x percentage', 'Percentage of the exon region of housekeeping genes that is covered at least 20-fold.', 'float', 0),
(85, 'QC:2000104', 'housekeeping genes 30x percentage', 'Percentage of the exon region of housekeeping genes that is covered at least 30-fold.', 'float', 0),
(86, 'QC:2000105', 'housekeeping genes 50x percentage', 'Percentage of the exon region of housekeeping genes that is covered at least 50-fold.', 'float', 0),
(87, 'QC:2000106', 'housekeeping genes 100x percentage', 'Percentage of the exon region of housekeeping genes that is covered at least 100-fold.', 'float', 0),
(88, 'QC:2000107', 'housekeeping genes 200x percentage', 'Percentage of the exon region of housekeeping genes that is covered at least 200-fold.', 'float', 0),
(89, 'QC:2000108', 'housekeeping genes 500x percentage', 'Percentage of the exon region of housekeeping genes that is covered at least 500-fold.', 'float', 0),
(90, 'QC:2000109', 'covered gene count', 'Number of covered genes (TPM >= 10.0)', 'int', 0),
(91, 'QC:2000110', 'aberrant spliced gene count', 'Number of aberrant spliced genes (>= 5%)', 'int', 0),
(92, 'QC:2000111', 'outlier gene count', 'Number of outlier genes (zscore >= 3.0)', 'int', 0);
