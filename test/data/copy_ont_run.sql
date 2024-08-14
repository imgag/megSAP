
INSERT INTO user (user_id, password, user_role, name, email, created, active) VALUES
('ahklauo1', 'password', 'user', 'Orang Utan Klaus', 'ouk@med.uni-tuebingen.de', NOW(), 1);

INSERT INTO sender (name) VALUES
('Peppa Pig');

INSERT INTO device (type, name) VALUES
('Promethion', 'ONT');

INSERT INTO project (name, type, internal_coordinator_id, analysis) VALUES
('LR_Diagnostik_Genome','diagnostic', 2, 'variants');

INSERT INTO processing_system (name_short, name_manufacturer, shotgun, genome_id, type, umi_type) VALUES
('LR-ONT-SQK-LSK114', 'Oxford Nanopore Tech. Ligation Sequencing Kit V14 (SQK-LSK114)', 0, 1, 'lrGS', 'n/a');

INSERT INTO sequencing_run (name, fcid, start_date, end_date, device_id, recipe, status) VALUES
('#01234', 'PAW01234', '2000-01-01', '2000-01-04', 1, '1D', 'run_finished');

INSERT INTO sample (id, name, sample_type, species_id, gender, tumor, ffpe, sender_id, disease_group) VALUES 
(1, '21073LRa277', 'DNA', 1, 'male', '0', '0', 1, 'n/a');

INSERT INTO processed_sample (sample_id, process_id, sequencing_run_id, lane, processing_system_id, project_id, normal_id) VALUES
(1,1,1,'1',1,1, NULL); 
