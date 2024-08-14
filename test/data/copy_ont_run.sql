
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

-- TODO: remove
CREATE  TABLE IF NOT EXISTS `runqc_ont`
(
  `id` INT(11) NOT NULL AUTO_INCREMENT,
  `sequencing_run_id` INT(11) NOT NULL,
  `read_num` INT NOT NULL,
  `yield` INT NOT NULL,
  `passing_filter_perc` FLOAT NOT NULL,
  `fraction_skipped` FLOAT NOT NULL,
  `q30_perc` FLOAT NOT NULL,
  `q20_perc` FLOAT NOT NULL,
  `n50` INT NOT NULL,
  `protocol_id`  VARCHAR(128) NOT NULL,
  `software_args` TEXT NOT NULL,
  `device_firmware_versions` VARCHAR(128) NOT NULL,
  `minknow_version` VARCHAR(16) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE INDEX `fk_sequencing_run_id` (`sequencing_run_id` ASC),
  CONSTRAINT `fk_sequencing_run_id`
    FOREIGN KEY (`sequencing_run_id`)
    REFERENCES `sequencing_run` (`id`)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION
)
ENGINE = InnoDB
DEFAULT CHARACTER SET = utf8;