<?php
/** 
	@page mvh_grz_export
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//returns case management data as XML
function get_cm_data($case_id)
{
	global $db_mvh;
	
	$lines = explode("\n", $db_mvh->getValue("SELECT cm_data FROM case_data WHERE id='{$case_id}'", ""));
	$lines_not_repeat = [];
	foreach($lines as $line)
	{
		$line = trim($line);
		if ($line=="") continue;
		if (strpos($line, "<redcap_repeat_instance></redcap_repeat_instance>")===false) continue;

		$lines_not_repeat[] = $line;
	}
	if (count($lines_not_repeat)!=1) trigger_error("Could not parse CM data. ".count($lines_not_repeat)." XML items which are not 'redcap_repeat_instance'!", E_USER_ERROR);
	
	return simplexml_load_string($lines_not_repeat[0], "SimpleXMLElement", LIBXML_NOCDATA);;
}

function not_empty($value, $element)
{
	if (trim($value)=="") trigger_error("Value '{$value}' of element {$element} must not be empty!", E_USER_ERROR);
	
	return $value;
}

function convert_gender($gender)
{
	if ($gender=="n/a") $gender = "unknown";
	return $gender;
}

function convert_tissue($tissue)
{
	if ($tissue=="blood") return "BTO:0000089";
	if ($tissue=="buccal mucosa") return "BTO:0003833";
	if ($tissue=="fibroblast") return "BTO:0000452";
	if ($tissue=="lymphocyte") return "BTO:0000775";
	if ($tissue=="skin") return "BTO:0001253";
	if ($tissue=="muscle") return "BTO:0000887";
		
	trigger_error("Could not convert tissue '{$tissue}' to BTO id!", E_USER_ERROR);
}

function convert_system_type($type)
{
	if ($type=="WES") return "wes";
	if ($type=="WGS") return "wgs";
	if ($type=="lrGS") return "wgs_lr";
		
	trigger_error("Unhandled processing system type '{$type}'!", E_USER_ERROR);
}

function convert_kit_manufacturer($name)
{
	if ($name=="TruSeqPCRfree") return "Illumina";
	if ($name=="twistCustomExomeV2") return "Twist";
	if ($name=="twistCustomExomeV2Covaris") return "Twist";
	
	trigger_error("Unhandled processing system name '{$name}'!", E_USER_ERROR);
}

function convert_sequencer_manufacturer($name)
{
	if ($name=="NovaSeq6000" || $name=="NovaSeqXPlus") return "Illumina";
	if ($name=="PromethION") return "ONT";
	
	trigger_error("Unhandled sequencer name '{$name}'!", E_USER_ERROR);
}

function megsap_version($gsvar)
{
	list($stdout) = exec2("grep '##PIPELINE=' {$gsvar}");
	foreach($stdout as $line)
	{
		$line = trim($line);
		if (starts_with($line, "##PIPELINE="))
		{
			return explode("=", $line, 2)[1];
		}
	}
	
	trigger_error("Could not determine megSAP version from '{$gsvar}'!", E_USER_ERROR);
}

//TODO check that the results are really the same as with the official QC pipeline bofore going live
function run_qc_pipeline($ps, $bam, $fq1, $fq2, $roi, $is_tumor)
{
	global $parser;
	global $seq_mode;
	global $is_somatic;
	global $qc_folder;
	$qc_wf_folder = "/mnt/storage2/MVH/tools/GRZ_QC_Workflow/";
	
	//run mosdepth if necessary
	$mosdepth_folder = "{$qc_folder}/mosdepth_".($is_tumor ? "tumor" : "germline")."/";
	$mosdepth_summary = "{$mosdepth_folder}/output_prefix.mosdepth.summary.txt";
	$mosdepth_regions = "{$mosdepth_folder}/output_prefix.regions.bed.gz";
	if (!file_exists($mosdepth_summary) || !file_exists($mosdepth_regions))
	{
		print "  running mosdepth for {$ps} in folder {$mosdepth_folder} ...\n";
		
		//make sure output folder exits
		exec2("mkdir -p {$mosdepth_folder}");
		
		//run
		$args = [];
		$args[] = "--threads 10";
		$args[] = "--by ".($roi!="" ? realpath($roi) : "{$qc_wf_folder}/assets/default_files/hg38_440_omim_genes.bed"); 
		$args[] = "--fasta /tmp/local_ngs_data_GRCh38/GRCh38.fa ";
		exec2("/mnt/storage2/MVH/tools/mosdepth ".implode(" ", $args)." {$mosdepth_folder}/output_prefix {$bam}");
	}
	else
	{
		print "  note: mosdepth results already exist in folder {$mosdepth_folder} - using it\n";
	}
	
	//run fastp if necessary
	$fastp_folder = "{$qc_folder}/fastp".($is_tumor ? "tumor" : "germline")."/";
	$fastp_json = "{$fastp_folder}/{$ps}.json";
	if (!file_exists($fastp_json))
	{
		print "  running fastp for {$ps} in folder {$fastp_folder} ...\n";
		
		//make sure output folder exits
		exec2("mkdir -p {$fastp_folder}");
		
		//run fastp
		$args = [];
		$args[] = "--thread 10";
		$args[] = "--in1 ".realpath($fq1);
		$args[] = "--in2 ".realpath($fq2);
		$args[] = "--out1 {$fastp_folder}/R1.fastq.gz";
		$args[] = "--out2 {$fastp_folder}/R2.fastq.gz";
		$args[] = "--json {$fastp_json}";
		$args[] = "--html {$fastp_folder}/{$ps}.html";
		$args[] = "--detect_adapter_for_pe";
		exec2("/mnt/storage2/MVH/tools/fastp ".implode(" ", $args));
	}
	else
	{
		print "  note: fastp results already exist in folder {$fastp_folder} - using it\n";
	}
	
	//generate report
	$report = "{$qc_folder}/{$ps}_report.csv";
	$args = [];
	$args[] = "--sample_id '{$ps}'";
	$args[] = "--labDataName 'blood DNA'";
	$args[] = "--libraryType ".strtolower($seq_mode);
	$args[] = "--sequenceSubtype ".($is_tumor ? "somatic" : "germline");
	$args[] = "--genomicStudySubtype ".($is_somatic ? "tumor+germline" : "germline-only");
	$args[] = "--fastp_json {$fastp_json}";
	$args[] = "--mosdepth_global_summary {$mosdepth_summary}";
	$args[] = "--mosdepth_target_regions_bed {$mosdepth_regions}";
	$args[] = " --thresholds {$qc_wf_folder}/assets/default_files/thresholds.json";
	$args[] = "--output {$report}";
	exec2("python3 {$qc_wf_folder}/bin/compare_threshold.py ".implode(" ", $args));
	
	//parse QC report
	$file = file($report);
	$headers = explode(",", trim($file[0]));
	$metrics = explode(",", trim($file[1]));
	return array_combine($headers, $metrics);
}

function get_read_length($ps, $info)
{
	global $db_ngsd;
	
	$sys_type = $info['sys_type'];
	if ($sys_type=="lrGS")
	{
		$read_length = get_processed_sample_qc($db_ngsd, $ps, "QC:2000131");	
	}
	else
	{
		$read_length = get_processed_sample_qc($db_ngsd, $ps, "QC:2000006");
		while (contains($read_length, "-")) $read_length = explode("-", $read_length, 2)[1];
	}
	
	return $read_length;
}


function create_files_json($files_to_submit, $info, $read_length)
{
	$output = [];
	foreach($files_to_submit as $file)
	{
		$file = realpath($file);
		$data = [];
		$data["filePath"] = basename($file);
		if (ends_with($file, ".fastq.gz")) $data["fileType"] = "fastq";
		else if (ends_with($file, ".bam")) $data["fileType"] = "bam";
		else if (ends_with($file, ".vcf")) $data["fileType"] = "vcf";
		else if (ends_with($file, ".bed")) $data["fileType"] = "bed";
		$data["checksumType"] = "sha256";
		$checksum_file = "{$qc_folder}/checksums/".basename($file).".sha256sum";
		if (!file_exists($checksum_file))
		{
			exec2("sha256sum {$file} > {$checksum_file}");
		}
		else
		{
			print "    note: checksum file for ".basename($file)." exists - using it\n";
		}
		$checksum = explode(" ", trim(implode(" ", file($checksum_file))))[0];
		$data["fileChecksum"] = $checksum;
		$data["fileSizeInBytes"] = filesize($file);
		if (ends_with($file, "_R1.fastq.gz"))
		{
			$data["readOrder"] = "R1";
		}
		else if (ends_with($file, "_R2.fastq.gz"))
		{
			$data["readOrder"] = "R2";
		}
		if (ends_with($file, "_R1.fastq.gz") || ends_with($file, "_R2.fastq.gz"))
		{
			$data["readLength"] = (int)$read_length;
			$data["flowcellId"] = $info["run_fcid"];
			$data["laneId"] =  implode(",", $info["ps_lanes"]);
		}
		
		$output[] = $data;
	}
	
	return $toutput;
}

function create_lab_data_json($files, $info, $grz_qc, $is_tumor)
{
	global $megsap_ver;
	global $db_ngsd;
	
	$output = [
				"labDataName" => "DNA ".($is_tumor ? "tumor" : "normal"),
				"tissueOntology" => [
						"name" => "BRENDA tissue ontology",
						"version" => "2021-10-26"
					],
				"tissueTypeId" => convert_tissue($info["tissue"]),
				"tissueTypeName" => $info["tissue"],
				"sampleDate" => not_empty($info["s_received"], "sampleDate"),
				"sampleConservation" => ($info["is_ffpe"] ? "ffpe" : "fresh-tissue"),
				"sequenceType" => "dna",
				"sequenceSubtype" => ($is_tumor ? "somatic" : "germline"),
				"fragmentationMethod" => ($info['sys_name_short']=="twistCustomExomeV2" ? "enzymatic" : "sonication"),
				"libraryType" => convert_system_type($info['sys_type']), 
				"libraryPrepKit" => $info["sys_name"],
				"libraryPrepKitManufacturer" => convert_kit_manufacturer($info['sys_name_short']), 
				"sequencerModel" => $info["device_type"],
				"sequencerManufacturer" => convert_sequencer_manufacturer($info["device_type"]),
				"kitName" => convert_sequencer_manufacturer($info["device_type"])." sequencing kit",
				"kitManufacturer" => convert_sequencer_manufacturer($info["device_type"]),
				"enrichmentKitManufacturer" => convert_kit_manufacturer($info['sys_name_short']),
				"enrichmentKitDescription" => $info["sys_name"],
				"barcode" => $info["sys_adapter1"]."/".$info["sys_adapter2"],
				"sequencingLayout" => "paired-end",
				"sequenceData"=> [
						"bioinformaticsPipelineName" => "megSAP",
						"bioinformaticsPipelineVersion" => $megsap_ver,
						"referenceGenome" => "GRCh38",
						"percentBasesAboveQualityThreshold" =>  [
								"minimumQuality" => (float)($grz_qc["qualityThreshold"]),
								"percent" => (float)($grz_qc["percentBasesAboveQualityThreshold"])
							], 
						"meanDepthOfCoverage" => (float)($grz_qc["meanDepthOfCoverage"]),
						"minCoverage" => (float)($grz_qc["minCoverage"]),
						"targetedRegionsAboveMinCoverage" => (float)(number_format($grz_qc["targetedRegionsAboveMinCoverage"],2)),
						"nonCodingVariants" => "true",
						"callerUsed" => [
								0 => [ //TODO remove hard-coded somatic caller when somatic small variants callset is implemented in NGSD
									"name" => $is_tumor ? "DRAGEN" : $db_ngsd->getValue("SELECT caller FROM small_variants_callset WHERE processed_sample_id='".$info["ps_id"]."'"),
									"version" => $is_tumor ? "4.3.13" : $db_ngsd->getValue("SELECT caller_version FROM small_variants_callset WHERE processed_sample_id='".$info["ps_id"]."'"),
									]
							],
						"files" => $files
					]
		 ];
		 
	if ($is_tumor)
	{
		$tcc = [];
		
		//tumor content - bioinformatics
		$tc = get_processed_sample_qc($db_ngsd, $info['ps_name'], "QC:2000054");	
		if (is_numeric($tc) && $tc>=0 && $tc<=100)
		{
			$tcc[] = ["count" => (float)$tc, "method" => "bioinformatics"];
		}
		
		if (count($tcc)>0)
		{
			$output["tumor cell count"] = $tcc;
		}
	}
	 
	return $output;
}

//parse command line arguments
$parser = new ToolBase("mvh_grz_export", "GRZ export for Modellvorhaben.");
$parser->addInt("case_id", "'id' in 'data_data' of 'MVH' database.", false);
$parser->addFlag("clear", "Clear export and QC folder before running this script.");
$parser->addFlag("test", "Test mode.");
extract($parser->parse($argv));

//init
$db_mvh = DB::getInstance("MVH");
$db_ngsd = DB::getInstance("NGSD");
$mvh_folder = get_path("mvh_folder");

//check that case ID is valid
$id = $db_mvh->getValue("SELECT id FROM case_data WHERE id='{$case_id}'");
if ($id=="") trigger_error("No case with id '{$case}' in MVH database!", E_USER_ERROR);

//clear
$folder = realpath($mvh_folder)."/grz_export/{$case_id}/";
$qc_folder = realpath($mvh_folder)."/grz_qc/{$case_id}/";
if ($clear) exec2("rm -rf {$folder} {$qc_folder}");

//parse case management data
$cm_data = get_cm_data($case_id);

//check network > determine KDK and study_subtype
$network = $cm_data->network_title;
$kdk = "";
$study_subtype = "";
$disease_type = "";
$is_somatic  = false;
if ($network=="Netzwerk Seltene Erkrankungen")
{
	$kdk = "KDKTUE002"; //NSE - Tübingen
	$study_subtype = "germline-only";
	$disease_type = "rare";
}
else if ($network=="Deutsches Netzwerk für Personalisierte Medizin")
{
	$kdk = "KDKTUE005"; //DNPM - Tübingen
	$study_subtype = "tumor+germline";
	$disease_type = "oncological";
	$is_somatic = true;
}
else if ($network=="Deutsches Konsortium Familiärer Brust- und Eierstockkrebs")
{
	$kdk = "KDKL00003"; //DK-FBREK - Leipzig
	$study_subtype = "germline";
	$disease_type = "hereditary";
}
else trigger_error("Unhandled network type '{$network}'!", E_USER_ERROR); 

//check seqencing mode
$seq_mode = $cm_data->seq_mode;
if ($test && $seq_mode=="") $seq_mode="WES"; //TODO
if ($seq_mode!="WGS" && $seq_mode!="WES") trigger_error("Unhandled seq_mode '{$seq_mode}'!", E_USER_ERROR);

print "case: {$case_id} (seq_mode: {$seq_mode} / network: {$network})\n";

//check germline processed sample is ok
$ps = $db_mvh->getValue("SELECT ps FROM case_data WHERE id='{$case_id}'");
$info = get_processed_sample_info($db_ngsd, $ps);
$sys = $info['sys_name_short'];
$patient_id = $info['patient_identifier'];
if ($patient_id=="") trigger_error("No patient identifier set for sample '{$ps}'!", E_USER_ERROR);
print "germline sample: {$ps} (system: {$sys})\n";

//check germline processed sample is ok
$ps_t = "";
$info_t = "";
if ($is_somatic)
{
	$ps_t = $db_mvh->getValue("SELECT ps_t FROM case_data WHERE id='{$case_id}'");
	if ($ps_t=="") trigger_error("Could no tumor sample for case '{$case_id}' in network '{$network}'!", E_USER_ERROR);
	
	//check tumor and germline have same system
	$info_t = get_processed_sample_info($db_ngsd, $ps_t);
	$sys_t = $info_t['sys_name_short'];
	if ($sys!=$sys_t) trigger_error("Mismatching tumor/normal processing systems for case '{$case_id}' in network '{$network}': '{$sys}' vs '{$sys_t}'", E_USER_ERROR);
	
	print "tumor sample: {$ps_t} (system: {$sys})\n";
}

//determine ROI
$roi = "";
if ($sys=="twistCustomExomeV2" || $sys=="twistCustomExomeV2Covaris") $roi = "{$mvh_folder}/rois/twist_exome_core_plus_refseq.bed";
if ($seq_mode!="WGS" && $roi=="") trigger_error("Could not determine target region for sample '{$ps}' with processing system '{$sys}'!", E_USER_ERROR);

//determine tanG==VNg
$sub_ids = $db_mvh->getValues("SELECT id FROM `submission_grz` WHERE status='pending' AND case_id='$case_id}'");
if (count($sub_ids)!=1)  trigger_error(count($sub_ids)." pending GRZ submissions for case {$case_id}. Must be one!", E_USER_ERROR);
$sub_id = $sub_ids[0];
$tang = $db_mvh->getValue("SELECT tang FROM submission_grz WHERE id='{$sub_id}'");

//determine megSAP version from germline GSvar file
$megsap_ver = megsap_version($info['ps_folder']."/{$ps}.GSvar");

//create export folder
exec2("mkdir -p {$folder}");
print "export folder: {$folder}\n";
exec2("mkdir -p {$folder}/files/");
exec2("mkdir -p {$folder}/metadata/");

//determine read length
$read_length = get_read_length($ps, $info);
if ($is_somatic)
{
	$read_length_t = get_read_length($ps_t, $info_t);
}

//TODO add support for lrGS
//TODO? add support for trio and RNA
//create germline raw data (FASTQs + germline VCF)
$n_bam = $info['ps_bam'];
$n_fq1 = "{$folder}/files/{$tang}_normal_R1.fastq.gz";
$n_fq2 = "{$folder}/files/{$tang}_normal_R2.fastq.gz";
$n_vcf = "{$folder}/files/{$tang}_normal.vcf";
if (!file_exists($n_fq1) || !file_exists($n_fq2) || !file_exists($n_vcf))
{
	$parser->execApptainer("ngs-bits", "BamToFastq", "-in {$n_bam} -out1 {$n_fq1} -out2 {$n_fq2} -extend {$read_length}", [$n_bam], ["{$folder}/files/"]);
	exec2("zcat ". $info['ps_folder']."/{$ps}_var.vcf.gz > {$n_vcf}");
}
else
{
	print "  note: FASTQ/VCF for germline {$ps} already exist - using them\n";
}
$files_to_submit = [$n_fq1, $n_fq2, $n_vcf];

//create tumor raw data (FASTQs, somatic VCF)
if ($is_somatic)
{
	$t_bam = $info_t['ps_bam'];
	$t_fq1 = "{$folder}/files/{$tang}_tumor_R1.fastq.gz";
	$t_fq2 = "{$folder}/files/{$tang}_tumor_R2.fastq.gz";
	$tn_vcf = "{$folder}/files/{$tang}_somatic.vcf";
	if (!file_exists($t_fq1) || !file_exists($t_fq2) || !file_exists($tn_vcf))
	{
		$parser->execApptainer("ngs-bits", "BamToFastq", "-in {$t_bam} -out1 {$t_fq1} -out2 {$t_fq2} -extend {$read_length_t}", [$t_bam], ["{$folder}/files/"]);
		exec2("zcat ". $info_t['ps_folder']."/../Somatic_{$ps_t}-{$ps}/{$ps_t}-{$ps}_var.vcf.gz > {$tn_vcf}");
	}
	else
	{
		print "  note: FASTQ/VCF for tumor {$ps_t} already exist - using them\n";
	}
	$files_to_submit_t = [$t_fq1, $t_fq2, $tn_vcf];
}

//copy target region
if ($roi!="")
{
	exec2("cp {$roi} {$folder}/files/target_regions.bed");
	$files_to_submit[] = "{$folder}/files/target_regions.bed";
	if ($is_somatic)
	{
		$files_to_submit_t[] = "{$folder}/files/target_regions.bed";
	}
}

//run QC pipeline for germline sample
exec2("mkdir -p {$qc_folder}/checksums/");
print "QC folder: {$qc_folder}\n";
$grz_qc = run_qc_pipeline($ps, $bam, $n_fq1, $n_fq2, $roi, false);

//run QC pipeline for somatic sample
if ($is_somatic)
{
	$grz_qc_t = run_qc_pipeline($ps_t, $t_bam, $t_fq1, $t_fq2, $roi, true);
}

//prepare file info for meta data JSON
print "  calculating file checksums...\n";
$files = create_files_json($files_to_submit, $info, $read_length);
if ($is_somatic)
{
	$files_t = create_files_json($files_to_submit_t, $info_t, $read_length_t);
}

//prepare lab data for meta data JSON
$lab_data = [];
$lab_data[] = create_lab_data_json($files, $info, $grz_qc, false);
if ($is_somatic)
{
	$lab_data[] = create_lab_data_json($files_t, $info_t, $grz_qc_t, true);
}

//create meta data JSON
print "  creating metadata...\n";
$json = [];
$json['$schema'] = "https://raw.githubusercontent.com/BfArM-MVH/MVGenomseq/refs/tags/v1.1.7/GRZ/grz-schema.json";
$json['submission'] = [
	"submissionDate" => $db_mvh->getValue("SELECT date FROM submission_grz WHERE id='{$sub_id}'"),
	"submissionType" => $db_mvh->getValue("SELECT type FROM submission_grz WHERE id='{$sub_id}'"),
	"tanG" => $tang,
	"localCaseId" => $patient_id,
	"coverageType" => "UNK", //TODO get from GenLab
	"submitterId" => "260840108",
	"genomicDataCenterId" => "GRZTUE002",
	"clinicalDataNodeId" => $kdk,
	"diseaseType" => $disease_type,
	"genomicStudyType" => "single",
	"genomicStudySubtype" => $study_subtype,
	"labName" => "Institute of Medical Genetics and Applied Genomics, Tübingen, Germany",
];
$json['donors'] = [
	0 => [
		"donorPseudonym" => $tang,
		"gender" => convert_gender($info["gender"]),
		"relation" => "index",
		"mvConsent" => [
			"presentationDate" => (string)($cm_data->mvconsentpresenteddate),
			"version" => (string)($cm_data->version_teilnahme),
			"scope" => [
					0 => [
						"type" => ((string)($cm_data->particip_4)=="Ja" ? "permit" : "deny"),
						"date" => (string)($cm_data->datum_teilnahme),
						"domain" => "mvSequencing"
						],
					1 => [
						"type" => ((string)($cm_data->particip_4_1)=="Ja" ? "permit" : "deny"),
						"date" => (string)($cm_data->datum_teilnahme),
						"domain" => "caseIdentification"
						],
					2 => [
						"type" => ((string)($cm_data->particip_4_2)=="Ja" ? "permit" : "deny"),
						"date" => (string)($cm_data->datum_teilnahme),
						"domain" => "reIdentification"
						],
				]
			],
		"researchConsents" => [
				0 => [
					"schemaVersion" => "2025.0.1",
					"presentationDate" => (string)($cm_data->bc_date),
					"scope" => ""
					]
			],
		"labData" => $lab_data		 
		]
	];

//write meta data JSON
file_put_contents("{$folder}/metadata/metadata.json", json_encode($json, JSON_PRETTY_PRINT));


//validate the submission
print "running grz-cli validate...\n";
$output = [];
$exit_code = null;
$stdout = "{$folder}/logs/grz_cli_validate.stdout";
$stderr = "{$folder}/logs/grz_cli_validate.stderr";
exec2("/mnt/storage2/MVH/tools/miniforge3/envs/grz-tools/bin/grz-cli validate --submission-dir {$folder} > {$stdout} 2> {$stderr}", $output, $exit_code); //when using exec2 the process hangs indefinitely sometimes
if ($exit_code!=0)
{
	trigger_error("grz-cli validate failed - see {$folder}/logs/ for output!\n", E_USER_ERROR);
}

//Encrypt the submission
print "running grz-cli encrypt...\n";
$config = ""; //TODO
if ($test) $config = "/mnt/storage2/MVH/config/config_test_phase.txt";
$stdout = "{$folder}/logs/grz_cli_encrypt.stdout";
$stderr = "{$folder}/logs/grz_cli_encrypt.stderr";
$output = [];
$exit_code = null;
exec("/mnt/storage2/MVH/tools/miniforge3/envs/grz-tools/bin/grz-cli encrypt --submission-dir {$folder} --config-file {$config} > {$stdout} 2> {$stderr}", $output, $exit_code); //when using exec2 the process hangs indefinitely sometimes
if ($exit_code!=0)
{
	print_r($output);
	trigger_error("grz-cli encrypt failed - see {$folder}/logs/ for output!\n", E_USER_ERROR);
}

//Upload the submission
print "running grz-cli upload...\n";
$stdout = "{$folder}/logs/grz_cli_upload.stdout";
$stderr = "{$folder}/logs/grz_cli_upload.stderr";
$output = [];
$exit_code = null;
exec("/mnt/storage2/MVH/tools/miniforge3/envs/grz-tools/bin/grz-cli upload --submission-dir {$folder} --config-file {$config} > {$stdout} 2> {$stderr}", $output, $exit_code); //when using exec2 the process hangs indefinitely sometimes
if ($exit_code!=0)
{
	print_r($output);
	trigger_error("grz-cli upload failed - see {$folder}/logs/ for output!\n", E_USER_ERROR);
}

//update submission status in database: MVH/submission_grz
$submission_id = trim(implode("", file($stdout)));
$db_mvh->executeStmt("UPDATE submission_grz SET status='done', submission_id='{$submission_id}' WHERE id='{$sub_id}'");
//TODO how do we see that the script failed and update MVH database? A wrapper around this script would be best...

//clean up export folder if successfull
exec2("rm -rf {$folder}");

/*
//TODO:
- run on SRV005 as bioinf when updated to Ubuntu 24.04
- store schema of MVH DB

#Test samples
- WGS: NA12878x3_20 - 99999
- WES T/N: NA12878x3_44 / NA12877_32 - 99998
- lrGS: 24067LRa008_01

#Installation of mosdepth

> wget https://github.com/brentp/mosdepth/releases/download/v0.3.11/mosdepth --no-check-certificate

#Installation of fastp

> wget http://opengene.org/fastp/fastp
> chmod a+x ./fastp

#Installation notes GRZ-CLI (see https://github.com/BfArM-MVH/grz-cli):

	- Install miniforge at /mnt/storage2/megSAP/tools/miniforge3/
		> curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
		> bash Miniforge3-$(uname)-$(uname -m).sh
	- Install GRZ-CLI
		> /mnt/storage2/megSAP/tools/miniforge3/bin/conda create -n grz-tools -c conda-forge -c bioconda "grz-cli"
		> /mnt/storage2/megSAP/tools/miniforge3/bin/conda activate grz-tools
	- Updates with:
		> /mnt/storage2/megSAP/tools/miniforge3/bin/conda update -n grz-tools "grz-cli"
	- Activate:
		> /mnt/storage2/megSAP/tools/miniforge3/bin/conda activate grz-tools
*/

?>
