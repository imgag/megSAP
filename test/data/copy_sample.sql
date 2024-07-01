
INSERT INTO user (user_id, password, user_role, name, email, created, active) VALUES
('ahklauo1', 'password', 'user', 'Orang Utan Klaus', 'ouk@med.uni-tuebingen.de', NOW(), 1);

INSERT INTO sender (name) VALUES
('Peppa Pig');

INSERT INTO device (type, name) VALUES
('GAIIx', 'test'),
('NovaSeqXPlus', 'big_iPod');

INSERT INTO project (name, type, internal_coordinator_id, analysis) VALUES
('Exome_Diagnostik','diagnostic', 2, 'variants'),
('IVac','research', 2, 'fastq'),
('SomaticAndTreatment','diagnostic', 2, 'variants');

INSERT INTO processing_system (name_short, name_manufacturer, shotgun, genome_id, type, umi_type) VALUES
('ssHAEv6', 'SureSelect Human All Exon v6', '1', 1, 'WES', 'n/a'),
('TruSeqPCRfree', 'TruSeq DNA PCR-Free', '1', 1, 'WGS', 'n/a'),
('TwistCustomExomeV2', 'Twist Custom Exome V2 IMGAG', '1', 1, 'WES', 'n/a'),
('nebRNAU2_mrna_UMI', 'NEBNext Ultra II Directional RNA mRNA UMI', '1', 1, 'RNA', 'IDT-UDI-UMI');

INSERT INTO sequencing_run (name, fcid, start_date, end_date, device_id, recipe) VALUES
('#00001', 'FCID4711', '2018-02-04', '2018-02-04', 1, '100+8+100'),
('#00002', 'FCID4712', '2018-02-05', '2018-02-05', 1, '100+8+100'),
('#01489', 'FCID4713', '2012-06-27', '2020-06-29', 1, '100+8+100'),
('#00123', 'FCID0001', '2023-08-23', '2023-08-23', 2, '151+10+10+151'),
('#00000', 'old_run', '2014-02-05', '2014-02-05', 1, '100+8+100');

INSERT INTO sample (id, name, sample_type, species_id, gender, tumor, ffpe, sender_id, disease_group) VALUES 
(1, 'DX180049', 'DNA', 1, 'male', '0', '0', 1, 'n/a'),
(2, 'DX180050', 'DNA', 1, 'male', '0', '0', 1, 'n/a'),
(3, 'FO180004', 'DNA', 1, 'male', '1', '1', 1, 'n/a'),
(4, 'FO180005', 'DNA', 1, 'male', '0', '0', 1, 'n/a'),
(5, 'DX181277', 'DNA', 1, 'male', '0', '0', 1, 'n/a'),
(6, 'DX181278', 'DNA', 1, 'male', '1', '1', 1, 'n/a'),
(7, 'DX181279', 'DNA', 1, 'male', '1', '1', 1, 'n/a'),
(8, 'DX181280', 'DNA', 1, 'male', '1', '1', 1, 'n/a'),
(9, 'DX203663', 'DNA', 1, 'male', '0', '0', 1, 'n/a'),
(10, 'DX203664', 'DNA', 1, 'female', '0', '0', 1, 'n/a'),
(11, 'DX203665', 'DNA', 1, 'male', '0', '0', 1, 'Diseases of the respiratory system'),
(12, 'RX123456', 'RNA', 1, 'male', '0', '0', 1, 'n/a');

INSERT INTO processed_sample (sample_id, process_id, sequencing_run_id, lane, processing_system_id, project_id, normal_id) VALUES
(1,1,1,'1',1,1, NULL),
(2,2,1,'1',1,1, NULL),
(3,1,1,'1',1,2, NULL),
(4,1,1,'1',1,2, NULL),
(5,1,1,'1',1,3, NULL),
(6,1,1,'1',1,3, 5),
(7,1,1,'1',1,3, NULL),
(8,1,1,'1',1,3, 5),
(9,2,3,'1',1,1, NULL),
(10,3,3,'1',1,1, NULL),
(11,4,3,'1',1,1, NULL),
(12,3,4,'6,8',4,3, NULL),
(1,5,4,'1,2,3,4',2,1, NULL),
(5,5,4,'5',1,3, NULL),
(6,5,4,'7',1,3, 5),
(1,4,5,'1,2,3,4',2,1, NULL); 

UPDATE processed_sample SET scheduled_for_resequencing=1 WHERE id=16;

INSERT INTO sample_relations (sample1_id, relation, sample2_id) VALUES
(8,'tumor-normal', 5),
(6, 'tumor-normal', 5);
