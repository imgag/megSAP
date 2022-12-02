INSERT INTO user (user_id, password, user_role, name, email, created, active) VALUES
('ahtesto1', 'password', 'user', 'Testing Testitus', 'test@med.uni-tuebingen.de', NOW(), 1);

INSERT INTO sender (name) VALUES
('Sending');

INSERT INTO device (type, name) VALUES
('GAIIx', 'test');

INSERT INTO project (name, type, internal_coordinator_id, analysis) VALUES
('SomaticAndTreatment','diagnostic', 2, 'variants');

INSERT INTO processing_system (name_short, name_manufacturer, shotgun, genome_id, type) VALUES
('ssHAEv6', 'SureSelect Human All Exon v6', '1', 1, 'WES');

INSERT INTO sequencing_run (name, fcid, start_date, end_date, device_id, recipe) VALUES
('#00001', 'FCID0001', '2022-11-11', '2022-11-11', 1, '100+8+100');

INSERT INTO sample (name, sample_type, species_id, gender, tumor, ffpe, sender_id, disease_group) VALUES 
('DNA220001', 'DNA', 1, 'male', '1', '1', 1, 'n/a'),
('DNA220002', 'DNA', 1, 'male', '0', '0', 1, 'n/a'),
('DNA220003', 'DNA', 1, 'male', '1', '1', 1, 'n/a'),
('DNA220004', 'DNA', 1, 'male', '0', '0', 1, 'n/a');

INSERT INTO processed_sample (sample_id, process_id, sequencing_run_id, lane, processing_system_id, project_id, normal_id) VALUES
(1,1,1,'1',1,1, NULL),
(2,1,1,'1',1,1, NULL),
(3,1,1,'1',1,1, NULL),
(4,1,1,'1',1,1, NULL);

INSERT INTO somatic_report_configuration (ps_tumor_id, ps_normal_id, created_by, created_date, last_edit_by) VALUES
(3, 4, 1, NOW(), 1);

INSERT INTO somatic_cnv_callset (ps_tumor_id, ps_normal_id, caller, caller_version, call_date, quality) VALUES
(3, 4, "ClinCNV", "test_version", NOW(), "n/a");

INSERT INTO somatic_cnv (somatic_cnv_callset_id, chr, start, end, cn, tumor_cn, tumor_clonality) VALUES
(1, "chr8",	213201,	30544791, 2, 3, 0.5),
(1, "chr8",	30547271,	34321971, 2, 3, 0.5),
(1, "chr8",	34322971,	43966949, 2, 3, 0.5),
(1, "chr8",	45272604,	72830077, 2, 3, 0.5),
(1, "chr8",	72830077,	84773558, 2, 3, 0.5),
(1, "chr8",	84774558,	130053710, 2, 3, 0.5),
(1, "chr8",	130357996,	145138635, 2, 3, 0.5),
(1, "chr6", 34887766, 59819533, 2, 3, 0.5),
(1, "chr6", 81061187, 170583767, 2, 3, 0.5),
(1, "chr9",	79038588,		83753379, 2, 3, 0.5),
(1, "chr9",	129008663,		133939239, 2, 3, 0.5),
(1, "chr10",	42281285,	133624078, 2, 3, 0.5),
(1, "chr11",	180188,		50898143, 2, 3, 0.5),
(1, "chr11",	53492736,	66192465, 2, 3, 0.5),
(1, "chr11",	66909869,	68047592, 2, 3, 0.5),
(1, "chr11",	68047626,	69738009, 2, 3, 0.5),
(1, "chr11",	69773288,	71147346, 2, 3, 0.5),
(1, "chr11",	71241501,	74054691, 2, 3, 0.5),
(1, "chr11",	74057385,	77702515, 2, 3, 0.5),
(1, "chr11",	77725524,	79367484, 2, 3, 0.5),
(1, "chr11",	79367484,	135086621, 2, 3, 0.5),
(1, "chr12",	66862,		15624427, 2, 3, 0.5),
(1, "chr12",	15631421,	16550201, 2, 3, 0.5),
(1, "chr12",	16551201,	22200531, 2, 3, 0.5),
(1, "chr12",	22201531,	24372030, 2, 3, 0.5),
(1, "chr12",	24372030,	25553114, 2, 3, 0.5),
(1, "chr12",	27081947,	30418573, 2, 3, 0.5),
(1, "chr12", 30418573, 31729028, 2, 3, 0.5),
(1, "chr12", 31754699, 31992555, 2, 3, 0.5),
(1, "chr12", 38905416, 47994497, 2, 3, 0.5),
(1, "chr12", 47995234, 48350082, 2, 3, 0.5),
(1, "chr12", 48351082, 51329853, 2, 3, 0.5),
(1, "chr12", 51920738, 52520316, 2, 3, 0.5),
(1, "chr12", 52544511, 60171936, 2, 3, 0.5),
(1, "chr12", 60372912, 63785092, 2, 3, 0.5),
(1, "chr12", 63802070, 133234376, 2, 3, 0.5);	





INSERT INTO somatic_report_configuration_cnv (somatic_report_configuration_id, somatic_cnv_id, exclude_artefact, exclude_low_tumor_content, exclude_low_copy_number, exclude_high_baf_deviation, exclude_other_reason, comment) VALUES
(1, 1, 1, 0, 0, 0, 0, "no comment"),
(1, 2, 1, 0, 0, 0, 0, "no comment"),
(1, 3, 1, 0, 0, 0, 0, "no comment"),
(1, 4, 1, 0, 0, 0, 0, "no comment"),
(1, 5, 1, 0, 0, 0, 0, "no comment"),
(1, 6, 1, 0, 0, 0, 0, "no comment"),
(1, 7, 1, 0, 0, 0, 0, "no comment"),
(1, 8, 1, 0, 0, 0, 0, "no comment"),
(1, 9, 1, 0, 0, 0, 0, "no comment"),
(1, 10, 1, 0, 0, 0, 0, "no comment"),
(1, 11, 1, 0, 0, 0, 0, "no comment"),
(1, 12, 1, 0, 0, 0, 0, "no comment"),
(1, 13, 1, 0, 0, 0, 0, "no comment"),
(1, 14, 1, 0, 0, 0, 0, "no comment"),
(1, 15, 1, 0, 0, 0, 0, "no comment"),
(1, 16, 1, 0, 0, 0, 0, "no comment"),
(1, 17, 1, 0, 0, 0, 0, "no comment"),
(1, 18, 1, 0, 0, 0, 0, "no comment"),
(1, 19, 1, 0, 0, 0, 0, "no comment"),
(1, 20, 1, 0, 0, 0, 0, "no comment"),
(1, 21, 1, 0, 0, 0, 0, "no comment"),
(1, 22, 1, 0, 0, 0, 0, "no comment"),
(1, 23, 1, 0, 0, 0, 0, "no comment"),
(1, 24, 1, 0, 0, 0, 0, "no comment"),
(1, 25, 1, 0, 0, 0, 0, "no comment"),
(1, 26, 1, 0, 0, 0, 0, "no comment"),
(1, 27, 1, 0, 0, 0, 0, "no comment"),
(1, 28, 1, 0, 0, 0, 0, "no comment"),
(1, 29, 1, 0, 0, 0, 0, "no comment"),
(1, 30, 1, 0, 0, 0, 0, "no comment"),
(1, 31, 1, 0, 0, 0, 0, "no comment"),
(1, 32, 1, 0, 0, 0, 0, "no comment"),
(1, 33, 1, 0, 0, 0, 0, "no comment"),
(1, 34, 1, 0, 0, 0, 0, "no comment"),
(1, 35, 1, 0, 0, 0, 0, "no comment"),
(1, 36, 1, 0, 0, 0, 0, "no comment");