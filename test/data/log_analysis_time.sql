INSERT INTO processing_system (name_short, name_manufacturer, shotgun, genome_id, type) VALUES
('ssHAEv6', 'SureSelect Human All Exon v6', '1', 2, 'WES');

INSERT INTO device (type, name) VALUES
('GAIIx', 'test');

INSERT INTO sequencing_run (name, fcid, start_date, end_date, device_id, recipe) VALUES
('#00001', 'FCID0001', '2022-11-11', '2022-11-11', 1, '100+8+100');

INSERT INTO sender (name) VALUES
('Sending');

INSERT INTO sample (name, sample_type, species_id, gender, tumor, ffpe, sender_id, disease_group) VALUES 
('NA12878', 'DNA', 1, 'male', '1', '1', 1, 'n/a');

INSERT INTO project (name, type, internal_coordinator_id, analysis) VALUES
('GiaB project','test', 2, 'variants');

INSERT INTO processed_sample (sample_id, process_id, sequencing_run_id, lane, processing_system_id, project_id, normal_id) VALUES
(1,58,1,'1',1,1, NULL);