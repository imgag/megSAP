
INSERT INTO user (user_id, password, user_role, name, email, created, active) VALUES
('ahklauo1', 'password', 'user', 'Orang Utan Klaus', 'ouk@med.uni-tuebingen.de', NOW(), 1);

INSERT INTO sender (name) VALUES
('Peppa Pig');

INSERT INTO device (type, name) VALUES
('GAIIx', 'test');

INSERT INTO project (name, type, internal_coordinator_id, analysis) VALUES
('Exome_Diagnostik','diagnostic', 2, 'variants'),
('IVac','research', 2, 'fastq'),
('SomaticAndTreatment','diagnostic', 2, 'variants');

INSERT INTO processing_system (name_short, name_manufacturer, shotgun, genome_id, type) VALUES
('ssHAEv6', 'SureSelect Human All Exon v6', '1', 1, 'WES');

INSERT INTO sequencing_run (name, fcid, start_date, end_date, device_id, recipe) VALUES
('#00001', 'FCID4711', '2018-02-04', '2018-02-04', 1, '100+8+100');

INSERT INTO sample (name, sample_type, species_id, gender, tumor, ffpe, sender_id) VALUES 
('DX181277', 'DNA', 1, 'n/a', '0', '0', 1),
('DX181278', 'DNA', 1, 'male', '0', '0', 1),
('DX181279', 'DNA', 1, 'female', '0', '0', 1);

INSERT INTO processed_sample (sample_id, process_id, sequencing_run_id, lane, processing_system_id, project_id, normal_id) VALUES
(1,1,1,'1',1,1, NULL),
(2,1,1,'1',1,1, NULL),
(3,1,1,'1',1,1, NULL);