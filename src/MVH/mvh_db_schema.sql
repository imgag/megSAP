--
-- Tabellenstruktur für Tabelle `case_data`
--

CREATE TABLE `case_data` (
  `id` int(11) NOT NULL,
  `cm_id` varchar(20) NOT NULL,
  `cm_data` mediumtext NOT NULL COMMENT 'case managment data in XML format as proviced by the RedCap API',
  `se_id` text DEFAULT NULL COMMENT 'ID in SE RedCap',
  `se_data` mediumtext DEFAULT NULL COMMENT 'Entries in SE RedCap (several are possible)',
  `rc_data` text DEFAULT NULL COMMENT 'Research consent data from meDIC converted to XML',
  `rc_data_json` text DEFAULT NULL COMMENT 'Research consent data in JSON format as provided by meDIC',
  `sap_id` varchar(20) NOT NULL,
  `ps` varchar(23) DEFAULT NULL COMMENT 'germline sample',
  `ps_t` varchar(23) DEFAULT NULL COMMENT 'tumor sample for tumor-normal'
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;

--
-- Tabellenstruktur für Tabelle `submission_grz`
--

CREATE TABLE `submission_grz` (
  `id` int(11) NOT NULL,
  `case_id` int(11) NOT NULL,
  `date` date NOT NULL,
  `type` enum('test','initial','followup','addition','correction') NOT NULL,
  `tang` varchar(64) NOT NULL,
  `pseudog` varchar(64) NOT NULL DEFAULT '',
  `status` enum('pending','done','failed') NOT NULL,
  `submission_id` text DEFAULT NULL,
  `submission_output` text DEFAULT NULL,
  `metadata` mediumtext DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;

--
-- Tabellenstruktur für Tabelle `submission_kdk_se`
--

CREATE TABLE `submission_kdk_se` (
  `id` int(11) NOT NULL,
  `case_id` int(11) NOT NULL,
  `date` date NOT NULL,
  `type` enum('test','initial','followup','addition','correction') NOT NULL,
  `tank` varchar(64) NOT NULL,
  `pseudok` varchar(64) NOT NULL DEFAULT '',
  `status` enum('pending','done','failed') NOT NULL,
  `submission_id` text DEFAULT NULL,
  `submission_output` text DEFAULT NULL,
  `metadata` mediumtext DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;
