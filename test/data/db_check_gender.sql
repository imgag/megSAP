
INSERT INTO sender (name) VALUES ('Peppa Pig');

INSERT INTO device (type, name) VALUES ('GAIIx', 'test');

INSERT INTO project (name, type, internal_coordinator_id) VALUES ('X-Chr','research', 1);

INSERT INTO processing_system (name_short, name_manufacturer, shotgun, genome_id, type, target_file) VALUES ('ssHAEv6', 'SureSelect Human All Exon v6', '1', 1, 'Panel', '/mnt/storage2/megSAP/data/enrichment/ssX_2022_03_11.bed');

INSERT INTO sequencing_run (name, fcid, start_date, end_date, device_id, recipe) VALUES ('#002', 'FCID4711', '2013-02-04', '2013-02-04', 1, '76+6+76');

INSERT INTO sample (name, sample_type, species_id, gender, tumor, ffpe, sender_id) VALUES ('GS140127', 'DNA', 1, 'male', '0', '0', 1);

INSERT INTO processed_sample (sample_id, process_id, sequencing_run_id, lane, processing_system_id, project_id) VALUES (1,1,1,'1',1,1);