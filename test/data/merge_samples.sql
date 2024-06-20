INSERT INTO user (user_id, password, user_role, name, email, created, active) VALUES
('ahtesto1', 'password', 'user', 'Testing Testitus', 'test@med.uni-tuebingen.de', NOW(), 1);

INSERT INTO sender (name) VALUES
('Sending');

INSERT INTO device (type, name) VALUES
('GAIIx', 'test');


ALTER TABLE project MODIFY COLUMN type enum(
	'diagnostic',
	'research',
	'test',
	'external',
	'megSAP-Tests'
)
NOT NULL AFTER aliases;

INSERT INTO qc_terms (qcml_id, name, description, type, obsolete) VALUES
("QC00001", "test_qc_1", "test qc value for merged sample 1", "string", 0),
("QC00002", "test_qc_2", "test qc value for merged sample 2", "string", 0),
("QC:2000005", "read count", "Total number of reads (one cluster in a paired-end...", "int", 0);

INSERT INTO project (name, type, internal_coordinator_id, analysis) VALUES
('merge_samples','megSAP-Tests', 2, 'variants');

INSERT INTO processing_system (name_short, name_manufacturer, shotgun, genome_id, type) VALUES
('ssHAEv6', 'SureSelect Human All Exon v6', '1', 2, 'WES'),
('rna_sys', 'RNA processing system version great', '1', 2, 'RNA');

INSERT INTO sequencing_run (name, fcid, start_date, end_date, device_id, recipe) VALUES
('#00001', 'FCID0001', '2022-11-11', '2022-11-11', 1, '100+8+100');

INSERT INTO sample (name, sample_type, species_id, gender, tumor, ffpe, sender_id, disease_group) VALUES 
('DNA220001', 'DNA', 1, 'male', '1', '1', 1, 'n/a'),
('DNA220002', 'DNA', 1, 'male', '0', '0', 1, 'n/a'),
('RNA220003', 'RNA', 1, 'male', '1', '1', 1, 'n/a'),
('DNA220004', 'DNA', 1, 'male', '0', '0', 1, 'n/a'),
('RNA220005', 'RNA', 1, 'male', '1', '1', 1, 'n/a'),
('DNA220006', 'DNA', 1, 'male', '0', '0', 1, 'n/a'),
('DNA220007', 'DNA', 1, 'male', '1', '1', 1, 'n/a'),
('DNA220008', 'DNA', 1, 'male', '0', '0', 1, 'n/a');

INSERT INTO processed_sample (sample_id, process_id, sequencing_run_id, lane, processing_system_id, project_id, normal_id) VALUES
(1,1,1,'1',1,1, NULL),
(2,1,1,'1',1,1, NULL),
(3,1,1,'1',2,1, NULL),
(4,1,1,'1',1,1, NULL),
(5,1,1,'1',2,1, NULL),
(6,1,1,'1',1,1, NULL),
(7,1,1,'1',1,1, NULL),
(8,1,1,'1',1,1, NULL),
(1,2,1,'1',1,1, NULL),
(2,2,1,'1',1,1, NULL),
(3,2,1,'1',2,1, NULL),
(4,2,1,'1',1,1, NULL),
(5,2,1,'1',2,1, NULL),
(6,2,1,'1',1,1, NULL),
(7,2,1,'1',1,1, NULL),
(8,2,1,'1',1,1, NULL);

UPDATE processed_sample SET scheduled_for_resequencing=1 WHERE id=1;
UPDATE processed_sample SET scheduled_for_resequencing=1 WHERE id=3;
UPDATE processed_sample SET scheduled_for_resequencing=1 WHERE id=4;
UPDATE processed_sample SET scheduled_for_resequencing=1 WHERE id=6;
UPDATE processed_sample SET scheduled_for_resequencing=1 WHERE id=7;
UPDATE processed_sample SET scheduled_for_resequencing=1 WHERE id=8;

INSERT INTO processed_sample_qc (processed_sample_id, qc_terms_id, value) VALUES
(1, 1, "qc 1 value for ps_sample id 1"),
(2, 1, "qc 1 value for ps_sample id 2"),
(3, 1, "qc 1 value for ps_sample id 3"),
(4, 1, "qc 1 value for ps_sample id 4"),
(5, 1, "qc 1 value for ps_sample id 5"),
(6, 1, "qc 1 value for ps_sample id 6"),
(7, 1, "qc 1 value for ps_sample id 7"),
(8, 1, "qc 1 value for ps_sample id 8"),
(9, 1, "qc 1 value for ps_sample id 9"),
(10, 1, "qc 1 value for ps_sample id 10"),
(11, 1, "qc 1 value for ps_sample id 11"),
(12, 1, "qc 1 value for ps_sample id 12"),
(13, 1, "qc 1 value for ps_sample id 13"),
(14, 1, "qc 1 value for ps_sample id 14"),
(15, 1, "qc 1 value for ps_sample id 15"),
(16, 1, "qc 1 value for ps_sample id 16"),
(1, 2, "qc 2 value for ps_sample id 1"),
(2, 2, "qc 2 value for ps_sample id 2"),
(3, 2, "qc 2 value for ps_sample id 3"),
(4, 2, "qc 2 value for ps_sample id 4"),
(5, 2, "qc 2 value for ps_sample id 5"),
(6, 2, "qc 2 value for ps_sample id 6"),
(7, 2, "qc 2 value for ps_sample id 7"),
(8, 2, "qc 2 value for ps_sample id 8"),
(9, 2, "qc 2 value for ps_sample id 9"),
(10, 2, "qc 2 value for ps_sample id 10"),
(11, 2, "qc 2 value for ps_sample id 11"),
(12, 2, "qc 2 value for ps_sample id 12"),
(13, 2, "qc 2 value for ps_sample id 13"),
(14, 2, "qc 2 value for ps_sample id 14"),
(15, 2, "qc 2 value for ps_sample id 15"),
(16, 2, "qc 2 value for ps_sample id 16");

INSERT INTO expression_gene (id, symbol) VALUES
(1, "BRCA1"),
(2, "BRCA2"),
(3, "SRY"),
(4, "SKI");

INSERT INTO expression (processed_sample_id, symbol_id, tpm, raw) VALUES
(3, 1, 1, 35),
(3, 2, 1, 40),
(3, 3, 0, 0),
(5, 1, 4, 40),
(5, 2, 2, 20),
(5, 3, 4, 60),
(5, 4, 1, 5);

INSERT INTO expression_exon (processed_sample_id, chr, start, end, rpb, srpb, raw) VALUES
(3, "chr17", 43125271, 43125364, 0.2, 4, 28),
(3, "chr17", 43124017, 43124115, 0.1, 2, 15),
(3, "chr17", 43115726, 43115779, 0.22, 4.4, 25),
(5, "chr17", 43125271, 43125364, 0.1, 1.4, 20),
(5, "chr17", 43124017, 43124115, 0.3, 6, 30),
(5, "chr17", 43115726, 43115779, 0.2, 4, 25),
(5, "chr17", 43106456, 43106533, 0.1, 2, 17);

INSERT INTO variant (chr, start, end, ref, obs) VALUES
("chr1" , 1, 2, "A", "T"),
("chr2" , 1, 2, "A", "T"),
("chr3" , 1, 2, "A", "T"),
("chr4" , 1, 2, "A", "T"),
("chr5" , 1, 2, "A", "T"),
("chr6" , 1, 2, "A", "T"),
("chr7" , 1, 2, "A", "T"),
("chr8" , 1, 2, "A", "T"),
("chr9" , 1, 2, "A", "T"),
("chr10", 1, 2, "A", "T"),
("chr11", 1, 2, "A", "T"),
("chr12", 1, 2, "A", "T"),
("chr13", 1, 2, "A", "T"),
("chr14", 1, 2, "A", "T"),
("chr15", 1, 2, "A", "T"),
("chr16", 1, 2, "A", "T"),
("chr17", 1, 2, "A", "T"),
("chr18", 1, 2, "A", "T"),
("chr19", 1, 2, "A", "T"),
("chr20", 1, 2, "A", "T"),
("chr21", 1, 2, "A", "T"),
("chr22", 1, 2, "A", "T"),
("chrX" , 1, 2, "A", "T"),
("chrY" , 1, 2, "A", "T");


INSERT INTO somatic_report_configuration (ps_tumor_id, ps_normal_id, created_by, created_date, last_edit_by) VALUES
(7, 8, 1, NOW(), 1);


INSERT INTO detected_somatic_variant (processed_sample_id_tumor, processed_sample_id_normal, variant_id, variant_frequency, depth, quality_snp) VALUES
(1, 2, 1, 0.3, 150, 100),
(1, 2, 2, 0.3, 150, 100),
(1, 2, 3, 0.3, 150, 100),
(3, 4, 4, 0.3, 150, 100),
(3, 4, 5, 0.3, 150, 100),
(3, 4, 6, 0.3, 150, 100),
(3, 4, 7, 0.3, 150, 100),
(5, 6, 8, 0.3, 150, 100),
(5, 6, 9, 0.3, 150, 100),
(5, 6, 10, 0.3, 150, 100),
(5, 6, 11, 0.3, 150, 100),
(5, 6, 12, 0.3, 150, 100),
(7, 8, 13, 0.3, 150, 100),
(7, 8, 14, 0.3, 150, 100),
(7, 8, 15, 0.3, 150, 100),
(7, 8, 16, 0.3, 150, 100),
(7, 8, 17, 0.3, 150, 100),
(7, 8, 18, 0.3, 150, 100);


INSERT INTO somatic_cnv_callset (ps_tumor_id, ps_normal_id, caller, caller_version, call_date, quality) VALUES
(1, 2, "ClinCNV", "test_version", NOW(), "n/a"),
(5, 6, "ClinCNV", "test_version", NOW(), "n/a"),
(7, 8, "ClinCNV", "test_version", NOW(), "n/a");

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
(2, "chr9",	129008663,		133939239, 2, 3, 0.5),
(2, "chr10",	42281285,	133624078, 2, 3, 0.5),
(2, "chr11",	180188,		50898143, 2, 3, 0.5),
(2, "chr11",	53492736,	66192465, 2, 3, 0.5),
(2, "chr11",	66909869,	68047592, 2, 3, 0.5),
(2, "chr11",	68047626,	69738009, 2, 3, 0.5),
(2, "chr11",	69773288,	71147346, 2, 3, 0.5),
(2, "chr11",	71241501,	74054691, 2, 3, 0.5),
(2, "chr11",	74057385,	77702515, 2, 3, 0.5),
(2, "chr11",	77725524,	79367484, 2, 3, 0.5),
(2, "chr11",	79367484,	135086621, 2, 3, 0.5),
(2, "chr12",	66862,		15624427, 2, 3, 0.5),
(3, "chr12",	15631421,	16550201, 2, 3, 0.5),
(3, "chr12",	16551201,	22200531, 2, 3, 0.5),
(3, "chr12",	22201531,	24372030, 2, 3, 0.5),
(3, "chr12",	24372030,	25553114, 2, 3, 0.5),
(3, "chr12",	27081947,	30418573, 2, 3, 0.5),
(3, "chr12", 30418573, 31729028, 2, 3, 0.5),
(3, "chr12", 31754699, 31992555, 2, 3, 0.5),
(3, "chr12", 38905416, 47994497, 2, 3, 0.5),
(3, "chr12", 47995234, 48350082, 2, 3, 0.5),
(3, "chr12", 48351082, 51329853, 2, 3, 0.5),
(3, "chr12", 51920738, 52520316, 2, 3, 0.5),
(3, "chr12", 52544511, 60171936, 2, 3, 0.5),
(3, "chr12", 60372912, 63785092, 2, 3, 0.5),
(3, "chr12", 63802070, 133234376, 2, 3, 0.5);	



INSERT INTO report_configuration (processed_sample_id, created_by, created_date, last_edit_by) VALUES
(6, 1, NOW(), 1);


INSERT INTO detected_variant (processed_sample_id, variant_id, genotype, mosaic) VALUES
(2, 1,  "het", 0),
(2, 2,  "het", 0),
(2, 3,  "het", 0),
(4, 4,  "het", 0),
(4, 5,  "het", 0),
(4, 6,  "het", 0),
(4, 7,  "het", 0),
(6, 8,  "het", 0),
(6, 9,  "het", 0),
(6, 10, "het", 0),
(6, 11, "het", 0),
(6, 12, "het", 0),
(8, 13, "het", 0),
(8, 14, "het", 0),
(8, 15, "het", 0),
(8, 16, "het", 0),
(8, 17, "het", 0),
(8, 18, "het", 0);



INSERT INTO cnv_callset (processed_sample_id, caller, caller_version, call_date, quality) VALUES
(2, "ClinCNV", "test_version", NOW(), "n/a"),
(4, "ClinCNV", "test_version", NOW(), "n/a"),
(6, "ClinCNV", "test_version", NOW(), "n/a");

INSERT INTO cnv (cnv_callset_id, chr, start, end, cn, quality_metrics) VALUES
(1, "chr8",	213201,	30544791, 2, ""),
(1, "chr8",	30547271,	34321971, 2, ""),
(1, "chr8",	34322971,	43966949, 2, ""),
(1, "chr8",	45272604,	72830077, 2, ""),
(1, "chr8",	72830077,	84773558, 2, ""),
(1, "chr8",	84774558,	130053710, 2, ""),
(1, "chr8",	130357996,	145138635, 2, ""),
(1, "chr6", 34887766, 59819533, 2, ""),
(1, "chr6", 81061187, 170583767, 2, ""),
(1, "chr9",	79038588,	83753379, 2, ""),
(2, "chr9",	129008663,	133939239, 2, ""),
(2, "chr10", 42281285,	133624078, 2, ""),
(2, "chr11", 180188,	50898143, 2, ""),
(2, "chr11", 53492736,	66192465, 2, ""),
(2, "chr11", 66909869,	68047592, 2, ""),
(2, "chr11", 68047626,	69738009, 2, ""),
(2, "chr11", 69773288,	71147346, 2, ""),
(2, "chr11", 71241501,	74054691, 2, ""),
(2, "chr11", 74057385,	77702515, 2, ""),
(2, "chr11", 77725524,	79367484, 2, ""),
(2, "chr11", 79367484,	135086621, 2, ""),
(2, "chr12", 66862,		15624427, 2, ""),
(3, "chr12", 15631421,	16550201, 2, ""),
(3, "chr12", 16551201,	22200531, 2, ""),
(3, "chr12", 22201531,	24372030, 2, ""),
(3, "chr12", 24372030,	25553114, 2, ""),
(3, "chr12", 27081947,	30418573, 2, ""),
(3, "chr12", 30418573, 31729028, 2, ""),
(3, "chr12", 31754699, 31992555, 2, ""),
(3, "chr12", 38905416, 47994497, 2, ""),
(3, "chr12", 47995234, 48350082, 2, ""),
(3, "chr12", 48351082, 51329853, 2, ""),
(3, "chr12", 51920738, 52520316, 2, ""),
(3, "chr12", 52544511, 60171936, 2, ""),
(3, "chr12", 60372912, 63785092, 2, ""),
(3, "chr12", 63802070, 133234376, 2, "");

INSERT INTO sv_callset (processed_sample_id, caller, caller_version, call_date) VALUES
(2, "Manta", "test_version", NOW()),
(4, "Manta", "test_version", NOW()),
(6, "Manta", "test_version", NOW());

INSERT INTO sv_deletion (sv_callset_id, chr, start_min, start_max, end_min, end_max) VALUES
(1, "chr1", 1, 300, 700, 800),
(1, "chr2", 1, 300, 700, 800),
(2, "chr3", 1, 300, 700, 800),
(2, "chr4", 1, 300, 700, 800),
(2, "chr5", 1, 300, 700, 800),
(2, "chr6", 1, 300, 700, 800),
(3, "chr7", 1, 300, 700, 800),
(3, "chr8", 1, 300, 700, 800),
(3, "chr9", 1, 300, 700, 800);

INSERT INTO sv_duplication (sv_callset_id, chr, start_min, start_max, end_min, end_max) VALUES
(1, "chr1", 1, 300, 700, 800),
(1, "chr2", 1, 300, 700, 800),
(2, "chr3", 1, 300, 700, 800),
(2, "chr4", 1, 300, 700, 800),
(2, "chr5", 1, 300, 700, 800),
(2, "chr6", 1, 300, 700, 800),
(3, "chr7", 1, 300, 700, 800),
(3, "chr8", 1, 300, 700, 800),
(3, "chr9", 1, 300, 700, 800);

INSERT INTO sv_insertion (sv_callset_id, chr, pos, ci_lower, ci_upper) VALUES
(1, "chr1", 5000, 0, 25),
(1, "chr2", 5000, 0, 25),
(2, "chr3", 5000, 0, 25),
(2, "chr4", 5000, 0, 25),
(2, "chr5", 5000, 0, 25),
(2, "chr6", 5000, 0, 25),
(3, "chr7", 5000, 0, 25),
(3, "chr8", 5000, 0, 25),
(3, "chr9", 5000, 0, 25);

INSERT INTO sv_inversion (sv_callset_id, chr, start_min, start_max, end_min, end_max) VALUES
(1, "chr1", 1, 300, 700, 800),
(1, "chr2", 1, 300, 700, 800),
(2, "chr3", 1, 300, 700, 800),
(2, "chr4", 1, 300, 700, 800),
(2, "chr5", 1, 300, 700, 800),
(2, "chr6", 1, 300, 700, 800),
(3, "chr7", 1, 300, 700, 800),
(3, "chr8", 1, 300, 700, 800),
(3, "chr9", 1, 300, 700, 800);

INSERT INTO sv_translocation (sv_callset_id, chr1, start1, end1, chr2, start2, end2) VALUES
(1, "chr1", 1, 300, "chr1", 700, 800),
(1, "chr2", 1, 300, "chr2", 700, 800),
(2, "chr3", 1, 300, "chr3", 700, 800),
(2, "chr4", 1, 300, "chr4", 700, 800),
(2, "chr5", 1, 300, "chr5", 700, 800),
(2, "chr6", 1, 300, "chr6", 700, 800),
(3, "chr7", 1, 300, "chr7", 700, 800),
(3, "chr8", 1, 300, "chr8", 700, 800),
(3, "chr9", 1, 300, "chr9", 700, 800);