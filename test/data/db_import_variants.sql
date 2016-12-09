
INSERT INTO processing_system (name_short, name_manufacturer, shotgun, genome_id) VALUES ('ssX', 'SureSelectXT X-Chromosome', '1', 1);

INSERT INTO device (type, name) VALUES ('GAIIx', 'test');

INSERT INTO sequencing_run (name, fcid, start_date, end_date, device_id, recipe) VALUES ('#002', 'FCID4711', '2013-02-04', '2013-02-04', 1, '76+6+76');

INSERT INTO sender (name) VALUES ('John Doe');

INSERT INTO project (name, type, internal_coordinator_id) VALUES ('X-Chr','research', 1);

-- germline samples

INSERT INTO sample (name, sample_type, species_id, gender, tumor, ffpe, sender_id) VALUES ('GS120676', 'DNA', 1, 'male', '0', '0', 1);
INSERT INTO processed_sample (sample_id, process_id, sequencing_run_id, lane, processing_system_id, project_id) VALUES (1,1,1,'1',1,1);

INSERT INTO sample (name, sample_type, species_id, gender, tumor, ffpe, sender_id) VALUES ('GS120677', 'DNA', 1, 'male', '0', '0', 1);
INSERT INTO processed_sample (sample_id, process_id, sequencing_run_id, lane, processing_system_id, project_id) VALUES (2,1,1,'1',1,1);

INSERT INTO sample (name, sample_type, species_id, gender, tumor, ffpe, sender_id) VALUES ('GS120678', 'DNA', 1, 'male', '0', '0', 1);
INSERT INTO processed_sample (sample_id, process_id, sequencing_run_id, lane, processing_system_id, project_id) VALUES (3,1,1,'1',1,1);

-- somatic samples

INSERT INTO sample (name, sample_type, species_id, gender, tumor, ffpe, sender_id) VALUES ('GS130796', 'DNA', 1, 'male', '1', '0', 1);
INSERT INTO processed_sample (sample_id, process_id, sequencing_run_id, lane, processing_system_id, project_id) VALUES (4,1,1,'1',1,1);

INSERT INTO sample (name, sample_type, species_id, gender, tumor, ffpe, sender_id) VALUES ('GS130797', 'DNA', 1, 'male', '0', '0', 1);
INSERT INTO processed_sample (sample_id, process_id, sequencing_run_id, lane, processing_system_id, project_id) VALUES (5,1,1,'1',1,1);

