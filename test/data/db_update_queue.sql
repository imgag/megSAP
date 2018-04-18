
INSERT INTO user (user_id, password, user_role, name, email, created, active) VALUES
('ahklauo1', 'password', 'user', 'Orang Utan Klaus', 'ouk@med.uni-tuebingen.de', NOW(), 1);

INSERT INTO sender (name) VALUES
('Peppa Pig');

INSERT INTO device (type, name) VALUES
('GAIIx', 'test');

INSERT INTO project (name, type, internal_coordinator_id, analysis) VALUES
('Exome_Diagnostik','diagnostic', 2, 'annotation'),
('IVac','research', 2, 'fastq'),
('SomaticAndTreatment','diagnostic', 2, 'annotation');

INSERT INTO processing_system (name_short, name_manufacturer, shotgun, genome_id, type) VALUES
('ssHAEv6', 'SureSelect Human All Exon v6', '1', 1, 'WES'),
('WGS', 'WGS kit', '1', 1, 'WGS'),
('RNA', 'RNA kit', '1', 1, 'RNA');

INSERT INTO sequencing_run (name, fcid, start_date, end_date, device_id, recipe) VALUES
('#00001', 'FCID4711', '2018-02-04', '2018-02-04', 1, '100+8+100');

INSERT INTO sample (name, sample_type, species_id, gender, tumor, ffpe, sender_id) VALUES
('DX150001', 'DNA', 1, 'n/a', '0', '0', 1),
('DX181109', 'DNA', 1, 'n/a', '0', '0', 1),
('DX181110', 'DNA', 1, 'male', '0', '0', 1),
('DX181111', 'DNA', 1, 'female', '0', '0', 1),
('DX181116', 'DNA', 1, 'n/a', '0', '0', 1),
('DX181117', 'RNA', 1, 'n/a', '0', '0', 1);

INSERT INTO processed_sample (sample_id, process_id, sequencing_run_id, lane, processing_system_id, project_id, normal_id) VALUES
(1,1,1,'1',1,1, NULL),
(2,1,1,'1',1,1, NULL),
(3,1,1,'1',1,1, NULL),
(4,1,1,'1',1,1, NULL),
(5,1,1,'1',2,1, NULL), /* WGS */
(6,1,1,'1',3,1, NULL); /* RNA */

INSERT INTO analysis_job (type, `high_priority`, args, sge_id, sge_queue) VALUES
('single sample', 1, '-steps ma,vc,an', '999001', 'priority_srv018'),
('single sample', 1, '', '', ''),
('single sample', 0, '', '999002', 'default_srv018'),
('single sample', 0, '', '999003', 'default_srv018'),
('multi sample',  0, '', '', ''),
('trio',          0, '', '', ''),
('somatic',       0, '', '', ''),
('single sample', 0, '', '', ''), /* WGS */
('single sample', 0, '', '', ''); /* RNA */

INSERT INTO analysis_job_sample (analysis_job_id, processed_sample_id, info) VALUES 
(1, 1, ''),
(2, 2, ''),
(3, 3, ''),
(4, 4, ''),
(5, 2, 'affected'),
(5, 3, 'control'),
(6, 2, 'child'),
(6, 3, 'father'),
(6, 4, 'mother'),
(7, 2, 'tumor'),
(7, 3, 'normal'),
(8, 5, ''),
(9, 6, '');

INSERT INTO analysis_job_history (analysis_job_id, time, user_id, status, output) VALUES 
(1, '2017-01-01T00:00:00', 1, 'finished', ''),
(2, NOW(), 1, 'queued', ''),
(3, NOW(), 1, 'queued', ''),
(3, NOW(), 1, 'started', ''),
(3, NOW(), 1, 'cancel', ''),
(4, NOW(), 1, 'queued', ''),
(4, NOW(), 1, 'started', ''),
(5, NOW(), 1, 'queued', ''),
(6, NOW(), 1, 'queued', ''),
(7, NOW(), 1, 'queued', ''),
(8, NOW(), 1, 'queued', ''),
(9, NOW(), 1, 'queued', '');
