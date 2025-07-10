<?php
/** 
	@page mvh_grz_export
*/

require_once("mvh_functions.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

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

function run_qc_pipeline($ps, $bam, $fq1, $fq2, $roi, $is_tumor)
{
	global $parser;
	global $seq_mode;
	global $is_somatic;
	global $qc_folder;
	global $is_lrgs;
	global $patient_id;
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
		print "  note: mosdepth results already exist in folder {$mosdepth_folder} - using them\n";
	}
	
	//run fastp if necessary
	$fastp_folder = "{$qc_folder}/fastp_".($is_tumor ? "tumor" : "germline")."/";
	$fastp_json = "{$fastp_folder}/{$ps}.json";
	if (!file_exists($fastp_json))
	{
		print "  running fastp for {$ps} in folder {$fastp_folder} ...\n";
		
		//make sure output folder exits
		exec2("mkdir -p {$fastp_folder}");
		
		//run fastp
		$args = [];
		$args[] = "--thread 10";
		if ($is_lrgs)
		{
			$args[] = "--in ".realpath($fq1);
			$args[] = "--out {$fastp_folder}/R1.fastq.gz";
		}
		else
		{
			$args[] = "--in1 ".realpath($fq1);
			$args[] = "--in2 ".realpath($fq2);
			$args[] = "--out1 {$fastp_folder}/R1.fastq.gz";
			$args[] = "--out2 {$fastp_folder}/R2.fastq.gz";
			$args[] = "--detect_adapter_for_pe";
		}
		$args[] = "--json {$fastp_json}";
		$args[] = "--html {$fastp_folder}/{$ps}.html";
		exec2("/mnt/storage2/MVH/tools/fastp".($is_lrgs ? "long" : "")." ".implode(" ", $args));
	}
	else
	{
		print "  note: fastp results already exist in folder {$fastp_folder} - using them\n";
	}
	
	//generate report
	$report = "{$qc_folder}/{$ps}_report.csv";
	$args = [];
	$args[] = "--donorPseudonym '{$patient_id}'";
	$args[] = "--sample_id '{$ps}'";
	$args[] = "--labDataName 'blood DNA'";
	$args[] = "--libraryType ".($is_lrgs ? "wgs_lr" : strtolower($seq_mode));
	$args[] = "--sequenceSubtype ".($is_tumor ? "somatic" : "germline");
	$args[] = "--genomicStudySubtype ".($is_somatic ? "tumor+germline" : "germline-only");
	$args[] = "--fastp_json {$fastp_json}";
	$args[] = "--mosdepth_global_summary {$mosdepth_summary}";
	$args[] = "--mosdepth_target_regions_bed {$mosdepth_regions}";
	$args[] = "--output {$report}";
	$args[] = "--meanDepthOfCoverage 1"; //dummy value, not used on our case
	$args[] = "--targetedRegionsAboveMinCoverage 1"; //dummy value, not used on our case
	$args[] = "--percentBasesAboveQualityThreshold 1"; //dummy value, not used on our case
	exec2("/mnt/storage2/MVH/tools/python3/bin/python3 {$qc_wf_folder}/bin/compare_threshold.py ".implode(" ", $args));
	
	//parse QC report
	$file = file($report);
	$headers = explode(",", trim($file[0]));
	$metrics = explode(",", trim($file[1]));
	return array_combine($headers, $metrics);
}

function get_read_length($ps, $sys_type)
{
	global $db_ngsd;
	
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
	global $qc_folder;
	
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
		if (ends_with($file, "_R1.fastq.gz") || ends_with($file, "_R2.fastq.gz") || ends_with($file, ".bam"))
		{
			$data["readLength"] = (int)$read_length;
			$data["flowcellId"] = $info["run_fcid"];
			$data["laneId"] =  implode(",", $info["ps_lanes"]);
		}
		
		$output[] = $data;
	}
	
	return $output;
}

function create_lab_data_json($files, $info, $grz_qc, $is_tumor)
{
	global $megsap_ver;
	global $db_ngsd;
	global $test;
	global $is_lrgs;
	
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
				"enrichmentKitManufacturer" => convert_kit_manufacturer($info['sys_name_short'], false),
				"enrichmentKitDescription" => $info["sys_name"],
				"barcode" => $info["sys_adapter1"]."/".$info["sys_adapter2"],
				"sequencingLayout" => ($is_lrgs ? "single-end" : "paired-end"),
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
		
		if($test && count($tcc)==0)
		{
			$tcc[] = ["count" => 47.11, "method" => "bioinformatics"];
		}
		
		if (count($tcc)>0)
		{
			$output["tumorCellCount"] = $tcc;
		}
	}
	 
	return $output;
}

//parse command line arguments
$parser = new ToolBase("mvh_grz_export", "GRZ export for Modellvorhaben.");
$parser->addInt("case_id", "'id' in 'case_data' of 'MVH' database.", false);
$parser->addFlag("clear", "Clear export and QC folder before running this script.");
$parser->addFlag("test", "Test mode.");
extract($parser->parse($argv));

//init
$db_mvh = DB::getInstance("MVH");
$db_ngsd = DB::getInstance("NGSD");
$mvh_folder = get_path("mvh_folder");

//check that case ID is valid
$id = $db_mvh->getValue("SELECT id FROM case_data WHERE id='{$case_id}'");
if ($id=="") trigger_error("No case with id '{$case_id}' in MVH database!", E_USER_ERROR);
$cm_id = $db_mvh->getValue("SELECT cm_id FROM case_data WHERE id='{$case_id}'");

//clear
$folder = realpath($mvh_folder)."/grz_export/{$case_id}/";
$qc_folder = realpath($mvh_folder)."/grz_qc/{$case_id}/";
if ($clear) exec2("rm -rf {$folder} {$qc_folder}");

//get data we need from MVH database
$cm_data = get_cm_data($db_mvh, $case_id);
$gl_data = get_gl_data($db_mvh, $case_id);

//check network > determine KDK and study_subtype
$network = xml_str($cm_data->network_title);
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
$seq_mode = xml_str($cm_data->seq_mode);
if ($seq_mode!="WGS" && $seq_mode!="WES") trigger_error("Unhandled seq_mode '{$seq_mode}'!", E_USER_ERROR);

//get patient identifer (pseudonym from case management) - this is the ID that is used to identify submissions from the same case by GRZ/KDK
$patient_id = xml_str($cm_data->psn); //TODO change to CM Fallnummer and pseudonymize via meDIC
if ($patient_id=="") trigger_error("No patient identifier set for sample '{$ps}'!", E_USER_ERROR);

print "MVH DB id: {$case_id} (CM ID: {$cm_id} / CM pseudonym: {$patient_id} / seq_mode: {$seq_mode} / network: {$network})\n";

//check germline processed sample is ok
$ps = $db_mvh->getValue("SELECT ps FROM case_data WHERE id='{$case_id}'");
$info = get_processed_sample_info($db_ngsd, $ps);
$sys = $info['sys_name_short'];
$sys_type = $info['sys_type'];
$is_lrgs = $sys_type=="lrGS";
print "germline sample: {$ps} (system: {$sys}, type: {$sys_type})\n";

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
$sub_ids = $db_mvh->getValues("SELECT id FROM `submission_grz` WHERE status='pending' AND case_id='{$case_id}'");
if (count($sub_ids)!=1)  trigger_error(count($sub_ids)." pending GRZ submissions for case {$case_id}. Must be one!", E_USER_ERROR);
$sub_id = $sub_ids[0];
print "ID in submission_grz table: {$sub_id}\n";
$tan_g = $db_mvh->getValue("SELECT tang FROM submission_grz WHERE id='{$sub_id}'");
print "TAN: {$tan_g}\n";

//determine megSAP version from germline GSvar file
$megsap_ver = megsap_version($info['ps_folder']."/{$ps}.GSvar");

//create export folder
exec2("mkdir -p {$folder}");
print "export folder: {$folder}\n";
exec2("mkdir -p {$folder}/files/");
exec2("mkdir -p {$folder}/metadata/");
exec2("mkdir -p {$folder}/logs/");

//determine read length
$read_length = get_read_length($ps, $sys_type);
if ($is_somatic)
{
	$read_length_t = get_read_length($ps_t, $sys_type);
}

//TODO add support for trio
//TODO add support for RNA?
//create germline raw data (FASTQs + germline VCF)
$n_bam = $info['ps_bam'];
$n_fq1 = "{$folder}/files/{$tan_g}_normal_R1.fastq.gz";
$n_fq2 = $is_lrgs ? "" : "{$folder}/files/{$tan_g}_normal_R2.fastq.gz";
$n_vcf = "{$folder}/files/{$tan_g}_normal.vcf";
$lrgs_bam = "{$folder}/files/{$tan_g}_normal.bam";
if ($is_lrgs && !file_exists($lrgs_bam)) //for lrGS we submit BAM: convert CRAM to BAM
{
	print "  generating BAM file for germline sample {$ps}...\n";
	$parser->execTool("Tools/cram_to_bam.php", "-cram {$n_bam} -bam {$lrgs_bam} -threads 10");
}
if (!$is_lrgs && (!file_exists($n_fq1) || !file_exists($n_fq2)))
{
	print "  generating FASTQ files for germline sample {$ps}...\n";
	$parser->execApptainer("ngs-bits", "BamToFastq", "-in {$n_bam} -out1 {$n_fq1} -out2 {$n_fq2} -extend {$read_length}", [$n_bam], ["{$folder}/files/"]);
}
if ($is_lrgs && !file_exists($n_fq1))
{
	print "  generating FASTQ file for germline sample {$ps} (needed for QC only)...\n";
	$parser->execApptainer("ngs-bits", "BamToFastq", "-in {$lrgs_bam} -out1 {$n_fq1}", [$lrgs_bam], ["{$folder}/files/"]);
}
if (!file_exists($n_vcf))
{
	print "  generating VCF file for germline sample {$ps}...\n";
	exec2("zcat ". $info['ps_folder']."/{$ps}_var.vcf.gz > {$n_vcf}");
}
$files_to_submit = $is_lrgs ? [$lrgs_bam, $n_vcf] : [$n_fq1, $n_fq2, $n_vcf];

//create tumor raw data (FASTQs, somatic VCF)
if ($is_somatic)
{
	$t_bam = $info_t['ps_bam'];
	$t_fq1 = "{$folder}/files/{$tan_g}_tumor_R1.fastq.gz";
	$t_fq2 = "{$folder}/files/{$tan_g}_tumor_R2.fastq.gz";
	$tn_vcf = "{$folder}/files/{$tan_g}_somatic.vcf";
	if (!file_exists($t_fq1) || !file_exists($t_fq2))
	{
		print "  generating FASTQ files for tumor sample {$ps_t}...\n";
		$parser->execApptainer("ngs-bits", "BamToFastq", "-in {$t_bam} -out1 {$t_fq1} -out2 {$t_fq2} -extend {$read_length_t}", [$t_bam], ["{$folder}/files/"]);
	}
	if (!file_exists($tn_vcf))
	{
		print "  generating VCF file for tumor sample {$ps_t}...\n";
		exec2("zcat ". $info_t['ps_folder']."/../Somatic_{$ps_t}-{$ps}/{$ps_t}-{$ps}_var.vcf.gz > {$tn_vcf}");
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
$grz_qc = run_qc_pipeline($ps, $is_lrgs ? $lrgs_bam : $n_bam, $n_fq1, $n_fq2, $roi, false);

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

//prepare research consent data - for format see https://www.medizininformatik-initiative.de/Kerndatensatz/KDS_Consent_V2025/MII-IG-Modul-Consent-TechnischeImplementierung-FHIRProfile-Consent.html
$active_consent_count = 0;
$research_use_allowed = false;
$date = "";
$start = "";
$end = "";
$rc_data = get_rc_data($db_mvh, $case_id);
foreach($rc_data->consent as $consent)
{
	if ($consent->status!="active") continue;
	++$active_consent_count;
	
	$date = xml_str($consent->date);
	$start = xml_str($consent->start);
	$end = xml_str($consent->end);
	
	foreach($consent->permit as $permit)
	{
		if ($permit->code=="2.16.840.1.113883.3.1937.777.24.5.3.8" || $permit->code=="2.16.840.1.113883.3.1937.777.24.5.3.1")
		{
			$research_use_allowed = true;
		}
	}
}
if ($active_consent_count==0) trigger_error("No active consent found in MVH data:\n{$rc_data}", E_USER_ERROR);
if ($active_consent_count>1) trigger_error("More than one active consent found in MVH data:\n{$rc_data}", E_USER_ERROR);
$research_consent = [
	"status" => "active",
	"scope" => [
		"coding" => [
				0 => [
					"system" => "http://terminology.hl7.org/CodeSystem/consentscope",
					"code" => "research"					
				]
			]
		],
	"category" => [
			[
				"coding" => [
						[
							"system" => "http://loinc.org",
							"code" => "57016-8"					
						]
					]
			],
			[
				"coding" => [
						[
							"system" => "https://www.medizininformatik-initiative.de/fhir/modul-consent/CodeSystem/mii-cs-consent-consent_category",
							"code" => "2.16.840.1.113883.3.1937.777.24.2.184"					
						]
					]
			]
		],
	"patient" => [
		"reference" => $patient_id
		],
	"dateTime" => $date,
	"policy" => [
			[
				"uri" => "urn:oid:2.16.840.1.113883.3.1937.777.24.2.1791"
			]
		],
	"provision" => [
		"type" => "deny",
		"period" => [
			"start" => $start,
			"end" => $end
			],
		"provision" => [
			[
				"code" => [
					[
						"coding" => [
							[
								"system" => "urn:oid:2.16.840.1.113883.3.1937.777.24.5.3",
								"code" => "2.16.840.1.113883.3.1937.777.24.5.3.1",
								"display" => "PATDAT_erheben_speichern_nutzen"
							]
						]
					]
				],
				"type" => ($research_use_allowed ? "permit" : "deny"),
				"period" => [
					"start" => $start,
					"end" => $end
					],
			]
		]
	]
];

//create meta data JSON
print "  creating metadata...\n";
$json = [];
$json['$schema'] = "https://raw.githubusercontent.com/BfArM-MVH/MVGenomseq/refs/tags/v1.1.7/GRZ/grz-schema.json";
$json['submission'] = [
	"submissionDate" => $db_mvh->getValue("SELECT date FROM submission_grz WHERE id='{$sub_id}'"),
	"submissionType" => $db_mvh->getValue("SELECT type FROM submission_grz WHERE id='{$sub_id}'"),
	"tanG" => $tan_g,
	"localCaseId" => $patient_id,
	"coverageType" => convert_coverage($gl_data->accounting_mode),
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
		"donorPseudonym" => $patient_id,
		"gender" => convert_gender($info["gender"]),
		"relation" => "index",
		"mvConsent" => [
			"presentationDate" => xml_str($cm_data->mvconsentpresenteddate),
			"version" => xml_str($cm_data->version_teilnahme),
			"scope" => [
					0 => [
						"type" => (xml_str($cm_data->particip_4)=="Ja" ? "permit" : "deny"),
						"date" => xml_str($cm_data->datum_teilnahme),
						"domain" => "mvSequencing"
						],
					1 => [
						"type" => (xml_str($cm_data->particip_4_1)=="Ja" ? "permit" : "deny"),
						"date" => xml_str($cm_data->datum_teilnahme),
						"domain" => "caseIdentification"
						],
					2 => [
						"type" => (xml_str($cm_data->particip_4_2)=="Ja" ? "permit" : "deny"),
						"date" => xml_str($cm_data->datum_teilnahme),
						"domain" => "reIdentification"
						],
				]
			],
		"researchConsents" => [
				0 => [
					"schemaVersion" => "2025.0.1",
					"presentationDate" => xml_str($cm_data->bc_date),
					"scope" => $research_consent
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
$config = "/mnt/storage2/MVH/config/config_".($test ? "test_phase" : "production").".txt";
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

//print submission ID for wrapper script
$submission_id = trim(implode("", file($stdout)));
print "SUBMISSION ID GRZ: {$submission_id}\n";

//clean up export folder if successfull
if (!$test)
{
	exec2("rm -rf {$folder} {$qc_folder}");
}

/*
//TODO:
- add tests when first final version is done
- remove test data when no longer needed:
	- WGS: NA12878x3_20 - 99999
	- WES T/N: NA12878x3_44 / NA12877_32 - 99998
	- lrGS: 24067LRa008_02 - 99997

#################### DOCU: Installation of tools (for running this script without GRZ QC workflow) #####################

#Installation of mosdepth

> wget https://github.com/brentp/mosdepth/releases/download/v0.3.11/mosdepth --no-check-certificate

#Installation of fastp

> wget http://opengene.org/fastp/fastp
> chmod a+x ./fastp

#Installation of fastplong

> wget http://opengene.org/fastplong/fastplong
> chmod a+x ./fastplong

#Installation notes GRZ-CLI (see https://github.com/BfArM-MVH/grz-tools/blob/main/packages/grz-cli/README.md):

	- Install miniforge at /mnt/storage2/megSAP/tools/miniforge3/
		> curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
		> bash Miniforge3-$(uname)-$(uname -m).sh
	- Install GRZ-CLI		
		> /mnt/storage2/MVH/tools/miniforge3/bin/conda create -n grz-tools -c conda-forge -c bioconda "grz-cli"
		> /mnt/storage2/MVH/tools/miniforge3/bin/conda activate grz-tools
	- Updates with:
		> /mnt/storage2/MVH/tools/miniforge3/bin/conda update -n base -c conda-forge conda
		> /mnt/storage2/MVH/tools/miniforge3/bin/conda update -n grz-tools -c conda-forge -c bioconda grz-cli
		> cd /mnt/storage2/MVH/tools/GRZ_QC_Workflow && git pull
	- Activate:
		> /mnt/storage2/MVH/tools/miniforge3/bin/conda activate grz-tools

#Installation of python3

> python3 -m venv /mnt/storage2/MVH/tools/python
> /mnt/storage2/MVH/tools/python3/bin/pip install grz-pydantic-models pandas argparse importlib


#################### DOCU: Installation of GRZ QC workflow #####################

#Installation nexflow

> cd /mnt/storage2/MVH/tools/nextflow
> curl -s https://get.nextflow.io | bash

#Installation of nf-core tools

> pipx install nf-core

#Installation of GRZ QC pipeline

> export PATH=$PATH:/mnt/storage2/MVH/tools/nextflow/
> nf-core pipelines download BfArM-MVH/GRZ_QC_Workflow --container-system singularity
> nextflow plugin install nf-schema@2.1.1

##################### DOCU: SQL database #####################

CREATE TABLE `case_data` (
 `id` int(11) NOT NULL AUTO_INCREMENT,
 `cm_id` varchar(20) NOT NULL,
 `cm_data` text NOT NULL COMMENT 'case managment data in XML format as proviced by the RedCap API',
 `se_id` text DEFAULT NULL COMMENT 'ID in SE RedCap',
 `se_data` text DEFAULT NULL COMMENT 'Entries in SE RedCap (several are possible)',
 `rc_data` text DEFAULT NULL COMMENT 'Research consent data from meDIC converted to XML',
 `rc_data_json` text DEFAULT NULL COMMENT 'Research consent data in JSON format as provided by meDIC',
 `gl_data` text DEFAULT NULL COMMENT 'GenLab data',
 `sap_id` varchar(20) NOT NULL,
 `ps` varchar(23) DEFAULT NULL COMMENT 'germline sample',
 `ps_t` varchar(23) DEFAULT NULL COMMENT 'tumor sample for tumor-normal',
 PRIMARY KEY (`id`),
 UNIQUE KEY `cm_id` (`cm_id`)
) ENGINE=InnoDB AUTO_INCREMENT=134673 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci

CREATE TABLE `submission_grz` (
 `id` int(11) NOT NULL AUTO_INCREMENT,
 `case_id` int(11) NOT NULL,
 `date` date NOT NULL,
 `type` enum('initial','followup','addition','correction') NOT NULL,
 `tang` varchar(64) NOT NULL,
 `status` enum('pending','done','failed') NOT NULL,
 `submission_id` text DEFAULT NULL,
 `submission_output` text DEFAULT NULL,
 PRIMARY KEY (`id`),
 KEY `case_id` (`case_id`),
 CONSTRAINT `submission_grz_ibfk_1` FOREIGN KEY (`case_id`) REFERENCES `case_data` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=12 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci

CREATE TABLE `submission_kdk_se` (
 `id` int(11) NOT NULL AUTO_INCREMENT,
 `case_id` int(11) NOT NULL,
 `date` date NOT NULL,
 `type` enum('initial','followup','addition','correction') NOT NULL,
 `tank` varchar(64) NOT NULL,
 `status` enum('pending','done','failed') NOT NULL,
 `submission_id` text DEFAULT NULL,
 `submission_output` text DEFAULT NULL,
 PRIMARY KEY (`id`),
 KEY `case_id` (`case_id`),
 CONSTRAINT `submission_kdk_se_ibfk_1` FOREIGN KEY (`case_id`) REFERENCES `case_data` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=2 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci
*/

?>
