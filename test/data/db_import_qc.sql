
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

INSERT INTO `qc_terms` (`id`, `qcml_id`, `name`, `description`, `obsolete`) VALUES
(1, 'QC:2000005', 'read count', 'Total number of reads (one cluster in a paired-end experiment generates two reads).', 0),
(2, 'QC:2000006', 'read length', 'Raw read length of a single read before trimming. Comma-separated list of lenghs if several.', 0),
(3, 'QC:2000007', 'Q20 read percentage', 'The percentage of reads with a mean base quality score greater than Q20.', 0),
(4, 'QC:2000008', 'Q30 base percentage', 'The percentage of bases with a minimum quality score of Q30.', 0),
(5, 'QC:2000009', 'no base call percentage', 'The percentage of bases without base call (N).', 0),
(6, 'QC:2000010', 'gc content percentage', 'The percentage of bases that are called to be G or C.', 0),
(7, 'QC:2000013', 'variant count', 'Total number of variants in the target region.', 0),
(8, 'QC:2000014', 'known variants percentage', 'Percentage of variants that are known polymorphisms in the dbSNP database.', 1),
(9, 'QC:2000015', 'high-impact variants percentage', 'Percentage of variants with high impact on the protein, i.e. stop-gain, stop-loss, frameshift, splice-acceptor or splice-donor variants.', 0),
(10, 'QC:2000016', 'homozygous variants percentage', 'Percentage of variants that are called as homozygous.', 0),
(11, 'QC:2000017', 'indel variants percentage', 'Percentage of variants that are insertions/deletions.', 0),
(12, 'QC:2000018', 'transition/transversion ratio', 'Transition/transversion ratio of single nucleotide variants.', 0),
(13, 'QC:2000019', 'trimmed base percentage', 'Percentage of bases that were trimmed during to adapter or quality trimming.', 0),
(14, 'QC:2000020', 'mapped read percentage', 'Percentage of reads that could be mapped to the reference genome.', 0),
(15, 'QC:2000021', 'on-target read percentage', 'Percentage of reads that could be mapped to the target region.', 0),
(16, 'QC:2000022', 'properly-paired read percentage', 'Percentage of properly paired reads (for paired-end reads only).', 1),
(17, 'QC:2000023', 'insert size', 'Average insert size (for paired-end reads only).', 1),
(18, 'QC:2000024', 'duplicate read percentage', 'Percentage of reads removed because they were duplicates (PCR, optical, etc).', 1),
(19, 'QC:2000025', 'target region read depth', 'Average sequencing depth in target region.', 1),
(20, 'QC:2000026', 'target region 10x percentage', 'Percentage the target region that is covered at least 10-fold.', 1),
(21, 'QC:2000027', 'target region 20x percentage', 'Percentage the target region that is covered at least 20-fold.', 1),
(22, 'QC:2000028', 'target region 30x percentage', 'Percentage the target region that is covered at least 30-fold.', 1),
(23, 'QC:2000029', 'target region 50x percentage', 'Percentage the target region that is covered at least 50-fold.', 1),
(24, 'QC:2000030', 'target region 100x percentage', 'Percentage the target region that is covered at least 100-fold.', 1),
(25, 'QC:2000031', 'target region 200x percentage', 'Percentage the target region that is covered at least 200-fold.', 1),
(26, 'QC:2000032', 'target region 500x percentage', 'Percentage the target region that is covered at least 500-fold.', 1),
(27, 'QC:2000033', 'error estimation read depth', 'Average read depth on the special target region used for error estimation after mapping.', 1),
(28, 'QC:2000034', 'error estimation N percentage', 'No base call (N) percentage determined on special target region after mapping.', 1),
(29, 'QC:2000035', 'error estimation SNV percentage', 'SNV error percentage determined on special target region after mapping.', 1),
(30, 'QC:2000036', 'error estimation indel percentage', 'indel error percentage determined on special target region after mapping.', 1),
(31, 'QC:2000039', 'gender check', 'Gender check result: ''n/a'', ''male_passed'', ''male_failed'', ''female_passed'', ''female_failed''.', 1);
