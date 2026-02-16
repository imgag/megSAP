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
		$args[] = "--fasta /tmp/local_ngs_data_GRCh38/GRCh38.fa";
		$args[] = "--fast-mode -F 772"; //TODO remove on 31.03.2026
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
	
	//determine QC thresholds
	$lib_type = ($is_lrgs ? "wgs_lr" : strtolower($seq_mode));
	$seq_subtype = ($is_tumor ? "somatic" : "germline");
	$study_subtype = ($is_somatic ? "tumor+germline" : "germline-only");
	$qc_data = json_decode(file_get_contents("/mnt/storage2/MVH/tools/miniforge3/envs/grz-tools/lib/python3.12/site-packages/grz_pydantic_models/resources/thresholds.json"), true);
	$thresholds = NULL;
	foreach($qc_data as $entry)
	{
		if ($entry['libraryType']!=$lib_type) continue;
		if ($entry['sequenceSubtype']!=$seq_subtype) continue;
		if ($entry['genomicStudySubtype']!=$study_subtype) continue;
		
		$thresholds = $entry['thresholds'];
	}
	if (is_null($thresholds)) trigger_error("Could not determine thresholds for {$lib_type}/{$seq_subtype}/{$$study_subtype}!", E_USER_ERROR);
	
	//generate report
	$report = "{$qc_folder}/{$ps}_report.csv";
	$args = [];
	$args[] = "--donorPseudonym '{$patient_id}'";
	$args[] = "--sample_id '{$ps}'";
	$args[] = "--labDataName 'blood DNA'";
	$args[] = "--libraryType {$lib_type}";
	$args[] = "--sequenceSubtype {$seq_subtype}";
	$args[] = "--genomicStudySubtype {$study_subtype}";
	$args[] = "--fastp_json {$fastp_json}";
	$args[] = "--mosdepth_global_summary {$mosdepth_summary}";
	$args[] = "--mosdepth_target_regions_bed {$mosdepth_regions}";
	$args[] = "--output {$report}";
	$args[] = "--meanDepthOfCoverage 1"; //dummy value, not used on our case
	$args[] = "--targetedRegionsAboveMinCoverage 1"; //dummy value, not used on our case
	$args[] = "--percentBasesAboveQualityThreshold 1"; //dummy value, not used on our case
	$args[] = "--meanDepthOfCoverageRequired ".$thresholds['meanDepthOfCoverage']; 
	$args[] = "--qualityThreshold ".$thresholds['percentBasesAboveQualityThreshold']['qualityThreshold'];
	$args[] = "--percentBasesAboveQualityThresholdRequired ".$thresholds['percentBasesAboveQualityThreshold']['percentBasesAbove'];
	$args[] = "--minCoverage ".$thresholds['targetedRegionsAboveMinCoverage']['minCoverage'];
	$args[] = "--targetedRegionsAboveMinCoverageRequired ".$thresholds['targetedRegionsAboveMinCoverage']['fractionAbove'];
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
		if (!file_exists($checksum_file) || trim(file_get_contents($checksum_file))=="")
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

function create_lab_data_json($files, $info, $grz_qc, $is_tumor, $info_germline=null)
{
	global $megsap_ver;
	global $db_ngsd;
	global $test;
	global $is_lrgs;
	
	//determine tissue ontology ID
	$tissue = "anatomical structure";
	$tissue_id = "UBERON:0000061";
	$ontology = "UBERON";
	$ontology_version = "2025-08-15";
	if (!$is_tumor)
	{
		$tmp = trim($info["tissue"]);
		if ($tmp!="" && $tmp!="n/a")
		{
			$tissue = $tmp;
			$tissue_id = convert_tissue($tissue);
			if (starts_with($tissue_id, "CL:")) //Cell Ontology
			{
				$ontology = "CL";
				$ontology_version = "2025-10-16";
			}
		}
	}
	else //somatic: take oncotree data from NGSD
	{
		$oncotree_codes = $db_ngsd->getValues("SELECT disease_info FROM sample_disease_info WHERE sample_id='".$info['s_id']."' and type='Oncotree code'");
		if (count($oncotree_codes)>0)
		{
			$tissue_id = $oncotree_codes[0];
			$tissue = $db_ngsd->getValue("SELECT name FROM oncotree_term WHERE oncotree_code='{$tissue_id}'");
			$ontology = "OncoTree";
			$ontology_version = "2021_11_02";
		}
	}
	
	//determine sample date
	$sample_date = trim($info["s_received"]);
	if ($sample_date=="") $sample_date = trim($info["order_date"]);
	
	$output = [
				"labDataName" => "DNA ".($is_tumor ? "tumor" : "normal"),
				"tissueOntology" => [
						"name" => $ontology,
						"version" => $ontology_version
					],
				"tissueTypeId" => $tissue_id,
				"tissueTypeName" => $tissue,
				"sampleDate" => not_empty($sample_date, "sampleDate"),
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
				"sequenceData" => [
						"bioinformaticsPipelineName" => "megSAP",
						"bioinformaticsPipelineVersion" => $megsap_ver,
						"referenceGenome" => "GRCh38",
						"percentBasesAboveQualityThreshold" =>  [
								"minimumQuality" => (float)($grz_qc["qualityThreshold"]),
								"percent" => (float)($grz_qc["percentBasesAboveQualityThreshold"])
							], 
						"meanDepthOfCoverage" => (float)($grz_qc["meanDepthOfCoverage"])*1.05, //TODO remove?
						"minCoverage" => (float)($grz_qc["minCoverage"]),
						"targetedRegionsAboveMinCoverage" => (float)(number_format($grz_qc["targetedRegionsAboveMinCoverage"],2)),
						"nonCodingVariants" => "true",
						"callerUsed" => [
								0 => [
									"name" => $is_tumor ? $db_ngsd->getValue("SELECT caller FROM somatic_snv_callset WHERE processed_sample_id_tumor='".$info["ps_id"]."' AND processed_sample_id_normal='".$info_germline["ps_id"]."'") : $db_ngsd->getValue("SELECT caller FROM small_variants_callset WHERE processed_sample_id='".$info["ps_id"]."'"),
									"version" => $is_tumor ? $db_ngsd->getValue("SELECT caller_version  FROM somatic_snv_callset WHERE processed_sample_id_tumor='".$info["ps_id"]."' AND processed_sample_id_normal='".$info_germline["ps_id"]."'") : $db_ngsd->getValue("SELECT caller_version FROM small_variants_callset WHERE processed_sample_id='".$info["ps_id"]."'"),
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
		
		//tumor content - pathoglogy
		$tc = $db_ngsd->getValue("SELECT disease_info FROM sample_disease_info WHERE type='tumor fraction' AND sample_id='".$info['s_id']."' LIMIT 1", ""); 
		if (is_numeric($tc) && $tc>=0 && $tc<=100)
		{
			$tcc[] = ["count" => (float)$tc, "method" => "pathology"];
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
$parser->addInt("cm_id", "ID in case management RedCap database.", false);
$parser->addFlag("clear", "Clear export and QC folder before running this script.");
$parser->addFlag("test", "Test mode.");
extract($parser->parse($argv));

//init
$db_mvh = DB::getInstance("MVH");
$db_ngsd = DB::getInstance("NGSD");
$mvh_folder = get_path("mvh_folder");
$threads = 10;

//check that case ID is valid
$id = $db_mvh->getValue("SELECT id FROM case_data WHERE cm_id='{$cm_id}'");
if ($id=="") trigger_error("No case with ID '{$cm_id}' in MVH database!", E_USER_ERROR);

//clear
$folder = realpath($mvh_folder)."/grz_export/{$cm_id}/";
$qc_folder = realpath($mvh_folder)."/grz_qc/{$cm_id}/";
if ($clear) exec2("rm -rf {$folder} {$qc_folder}");

//get data we need from MVH database
$cm_data = get_cm_data($db_mvh, $id);

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
	$study_subtype = "germline-only";
	$disease_type = "hereditary";
}
else trigger_error("Unhandled network type '{$network}'!", E_USER_ERROR); 

//check seqencing mode
$seq_mode = xml_str($cm_data->seq_mode);
if ($seq_mode!="WGS" && $seq_mode!="lrGS" && $seq_mode!="WES") trigger_error("Unhandled seq_mode '{$seq_mode}'!", E_USER_ERROR);

//start export
print "CM ID: {$cm_id} (MVH DB id: {$id} / CM Fallnummer: ".xml_str($cm_data->case_id)." / seq_mode: {$seq_mode} / network: {$network})\n";
print "Export start: ".date('Y-m-d H:i:s')."\n";
$time_start = microtime(true);

//check germline processed sample is ok
$ps = $db_mvh->getValue("SELECT ps FROM case_data WHERE id='{$id}'");
$info = get_processed_sample_info($db_ngsd, $ps);
$sys = $info['sys_name_short'];
$sys_type = $info['sys_type'];
$is_lrgs = $sys_type=="lrGS";
print "germline sample: {$ps} (system: {$sys}, type: {$sys_type})\n";
if ($seq_mode!=$sys_type) trigger_error("Sequencing mode in CM RedCap is '{$seq_mode}', but processing system type is '{$sys_type}'!", E_USER_ERROR);

//check germline processed sample is ok
$ps_t = "";
$info_t = "";
if ($is_somatic)
{
	$ps_t = $db_mvh->getValue("SELECT ps_t FROM case_data WHERE id='{$id}'");
	if ($ps_t=="") trigger_error("Could no tumor sample for case '{$cm_id}' in network '{$network}'!", E_USER_ERROR);
	
	//check tumor and germline have same system
	$info_t = get_processed_sample_info($db_ngsd, $ps_t);
	$sys_t = $info_t['sys_name_short'];
	if ($sys!=$sys_t && !(starts_with($sys, "twistCustomExomeV2") && starts_with($sys_t, "twistCustomExomeV2"))) trigger_error("Mismatching tumor/normal processing systems for case '{$cm_id}' in network '{$network}': '{$sys}' vs '{$sys_t}'", E_USER_ERROR);
	
	print "tumor sample: {$ps_t} (system: {$sys})\n";
}

//determine ROI
$roi = "";
if ($sys=="twistCustomExomeV2" || $sys=="twistCustomExomeV2Covaris") $roi = "{$mvh_folder}/rois/twist_exome_core_plus_refseq.bed";
if ($seq_mode!="WGS" && $seq_mode!="lrGS" && $roi=="") trigger_error("Could not determine target region for sample '{$ps}' with processing system '{$sys}'!", E_USER_ERROR);

//determine tanG==VNg
$sub_ids = $db_mvh->getValues("SELECT id FROM `submission_grz` WHERE status='pending' AND case_id='{$id}'");
if (count($sub_ids)!=1)  trigger_error(count($sub_ids)." pending GRZ submissions for case {$cm_id}. Must be one!", E_USER_ERROR);
$sub_id = $sub_ids[0];
print "ID in submission_grz table: {$sub_id}\n";
$tan_g = $db_mvh->getValue("SELECT tang FROM submission_grz WHERE id='{$sub_id}'");
print "TAN: {$tan_g}\n";
$patient_id = $db_mvh->getValue("SELECT pseudog FROM submission_grz WHERE id='{$sub_id}'");
print "patient pseudonym: {$patient_id}\n";

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

print "base checks took ".time_readable(microtime(true)-$time_start)."\n";
$time_start = microtime(true);

//TODO add support for trio?
//create germline raw data (FASTQs + germline VCF)
$n_bam = $info['ps_bam'];
$n_fq1 = "{$folder}/files/normal_R1.fastq.gz";
$n_fq2 = $is_lrgs ? "" : "{$folder}/files/normal_R2.fastq.gz";
$n_vcf = "{$folder}/files/normal.vcf";
$lrgs_bam = "{$folder}/files/normal.bam";
if ($is_lrgs && !file_exists($lrgs_bam)) //for lrGS we submit BAM: convert CRAM to BAM
{
	print "  generating BAM file for germline sample {$ps}...\n";
	$parser->execTool("Tools/cram_to_bam.php", "-cram {$n_bam} -bam {$lrgs_bam} -threads {$threads}");
}
if ($is_lrgs) //fix BAM header if DS tag in @RG line is present more than once
{
	$ds_twice = false;
	list($stdout) = exec2("samtools view -H {$lrgs_bam}");  
	for ($i=0; $i<count($stdout); ++$i)
	{
		$line = $stdout[$i];
		
		if (substr_count($line, "\tDS:")>1)
		{
			$ds_twice = true;
			
			$new = [];
			$first_done = false;
			foreach(explode("\t", nl_trim($line)) as $entry)
			{
				if (starts_with($entry, "DS:"))
				{
					if (!$first_done)
					{
						$new[] = $entry;
						$first_done = true;
					}
				}
				else
				{
					$new[] = $entry;
				}
			}
			
			$stdout[$i] = implode("\t", $new)."\n";
		}
	}
	if ($ds_twice)
	{
		print "  fixing BAM header {$ps}...\n";
		
		//store header
		$header = $parser->tempFile(".txt");
		file_put_contents($header, implode("\n", $stdout));
		
		//reheader to tmp
		$tmp_bam = "{$folder}/files/tmp.bam";
		$parser->execApptainer("samtools", "samtools reheader", "{$header} {$lrgs_bam} > {$tmp_bam}", [], [dirname($lrgs_bam)]);
		
		//replace BAM and index it
		exec2("rm -rf {$lrgs_bam}*");
		$parser->moveFile($tmp_bam, $lrgs_bam);
		$parser->execApptainer("samtools", "samtools index", "-@ {$threads} {$lrgs_bam}", [], [dirname($lrgs_bam)]);
	}
}
if (!$is_lrgs && (!file_exists($n_fq1) || !file_exists($n_fq2)))
{
	print "  generating FASTQ files for germline sample {$ps}...\n";
	$parser->execApptainer("ngs-bits", "BamToFastq", "-in {$n_bam} -out1 {$n_fq1} -out2 {$n_fq2}", [$n_bam], ["{$folder}/files/"]);
}
if ($is_lrgs && !file_exists($n_fq1))
{
	print "  generating FASTQ file for germline sample {$ps} (needed for QC only)...\n";
	$parser->execApptainer("ngs-bits", "BamToFastq", "-in {$lrgs_bam} -out1 {$n_fq1}", [$lrgs_bam], ["{$folder}/files/"]);
}
if (!file_exists($n_vcf))
{
	print "  generating VCF file for germline sample {$ps}...\n";
	$vcf = $info['ps_folder']."/{$ps}_var.vcf.gz";
	$parser->execApptainer("ngs-bits", "VcfReplaceSamples", "-in {$vcf} -out {$n_vcf} -ids {$ps}=SAMPLE_GERMLINE", [$vcf], ["{$folder}/files/"]);
}
$files_to_submit = $is_lrgs ? [$lrgs_bam, $n_vcf] : [$n_fq1, $n_fq2, $n_vcf];

//create tumor raw data (FASTQs, somatic VCF)
if ($is_somatic)
{
	$t_bam = $info_t['ps_bam'];
	$t_fq1 = "{$folder}/files/tumor_R1.fastq.gz";
	$t_fq2 = "{$folder}/files/tumor_R2.fastq.gz";
	$tn_vcf = "{$folder}/files/somatic.vcf";
	if (!file_exists($t_fq1) || !file_exists($t_fq2))
	{
		print "  generating FASTQ files for tumor sample {$ps_t}...\n";
		$parser->execApptainer("ngs-bits", "BamToFastq", "-in {$t_bam} -out1 {$t_fq1} -out2 {$t_fq2}", [$t_bam], ["{$folder}/files/"]);
	}
	if (!file_exists($tn_vcf))
	{
		print "  generating VCF file for tumor/normal pair {$ps_t}/{$ps}...\n";
		
		$vcf = dirname($info_t['ps_folder'])."/Somatic_{$ps_t}-{$ps}/{$ps_t}-{$ps}_var.vcf.gz";
		$parser->execApptainer("ngs-bits", "VcfReplaceSamples", "-in {$vcf} -out {$tn_vcf} -ids {$ps}=SAMPLE_GERMLINE,{$ps_t}=SAMPLE_TUMOR", [$vcf], ["{$folder}/files/"]);
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

print "creating upload files took ".time_readable(microtime(true)-$time_start)."\n";
$time_start = microtime(true);

//run QC pipeline for germline sample
exec2("mkdir -p {$qc_folder}/checksums/");
print "QC folder: {$qc_folder}\n";
$grz_qc = run_qc_pipeline($ps, $is_lrgs ? $lrgs_bam : $n_bam, $n_fq1, $n_fq2, $roi, false);

//run QC pipeline for somatic sample
if ($is_somatic)
{
	$grz_qc_t = run_qc_pipeline($ps_t, $t_bam, $t_fq1, $t_fq2, $roi, true);
}

print "running QC pipeline took ".time_readable(microtime(true)-$time_start)."\n";
$time_start = microtime(true);

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
	$lab_data[] = create_lab_data_json($files_t, $info_t, $grz_qc_t, true, $info);
}

//if SE > check if it is a BC for children
$is_kids_bc = false;
if (xml_str($cm_data->bc_signed)=="Ja" && $network=="Netzwerk Seltene Erkrankungen")
{
	$se_data_rep = get_se_data($db_mvh, $id, true);
	list($bc_type, $bc_item) = get_bc_data_se($se_data_rep);
	$is_kids_bc = starts_with($bc_type, "Kinder");
}

//prepare research consent data - for format see https://www.medizininformatik-initiative.de/Kerndatensatz/KDS_Consent_V2025/MII-IG-Modul-Consent-TechnischeImplementierung-FHIRProfile-Consent.html
if ($is_kids_bc)
{
	$se_data = get_se_data($db_mvh, $id);
	$research_consent = convert_se_kids_bc_to_fhir($bc_item, $se_data, $parser);
}
else
{
	$active_consent_count = 0;
	$research_use_allowed = false;
	$date = "";
	$start = "";
	$end = "";
	$rc_data = get_rc_data($db_mvh, $id);
	if ($rc_data!=false)
	{
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
	}
	if ($active_consent_count>1) trigger_error("More than one active consent found in MVH data:\n{$rc_data}", E_USER_ERROR);
	if (xml_str($cm_data->bc_signed)=="Ja" && $active_consent_count==0) trigger_error("Patient has signed BC, but no active consent found in MVH data\n{$rc_data}", E_USER_ERROR);
	
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
}



//create meta data JSON
print "  creating metadata...\n";
$submission_type = $db_mvh->getValue("SELECT type FROM submission_grz WHERE id='{$sub_id}'");
$json = [];
$json['$schema'] = "https://raw.githubusercontent.com/BfArM-MVH/MVGenomseq/refs/tags/v1.2.1/GRZ/grz-schema.json";
$json['submission'] = [
	"submissionDate" => date('Y-m-d'),
	"submissionType" => $submission_type,
	"tanG" => $tan_g,
	"localCaseId" => $patient_id,
	"coverageType" => convert_coverage($cm_data->coveragetype),
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
		"researchConsents" => [],
		"labData" => $lab_data		 
		]
	];

//add consent data
if (xml_str($cm_data->bc_signed)=="Ja") //consent signed
{
	$json['donors'][0]['researchConsents'][] = [
					"schemaVersion" => "2025.0.1",
					"presentationDate" => xml_str($cm_data->bc_date),
					"scope" => $research_consent
					];
}
else //consent not signed
{
	$reason_missing = convert_bc_missing(xml_str($cm_data->bc_reason_missing));
	if ($reason_missing=="patient-inability") $reason_missing = "patient unable to consent";
	else if ($reason_missing=="patient-refusal") $reason_missing = "patient refuses to sign consent";
	else if ($reason_missing=="consent-not-returned") $reason_missing = "patient did not return consent documents";
	else if ($reason_missing=="other-patient-reason") $reason_missing = "other patient-related reason";
	else trigger_error("Count not convert reason why BC is missing: '{$reason_missing}'!", E_USER_ERROR);
	$json['donors'][0]['researchConsents'][] = [
					"schemaVersion" => "2025.0.1",
					"presentationDate" => xml_str($cm_data->bc_date),
					"noScopeJustification" => $reason_missing
					];
}

//write meta data JSON
$json_file = "{$folder}/metadata/metadata.json";
file_put_contents($json_file, json_encode($json, JSON_PRETTY_PRINT));

print "creating JSON took ".time_readable(microtime(true)-$time_start)."\n";
$time_start = microtime(true);

//print grz-cli version info
$grz_cli = "/mnt/storage2/MVH/tools/miniforge3/envs/grz-tools/bin/grz-cli";
list($stdout) = exec2("{$grz_cli} --version");
print "grz-cli version: ".trim(implode("", $stdout))."\n";
print "\n";

//validate the submission
print "running grz-cli validate...\n";
$config = "/mnt/storage2/MVH/config/config_".($test ? "test_phase" : "production").".txt";
$output = [];
$exit_code = null;
$stdout = "{$folder}/logs/grz_cli_validate.stdout";
$stderr = "{$folder}/logs/grz_cli_validate.stderr";
exec2("{$grz_cli} validate --submission-dir {$folder} --config-file {$config} > {$stdout} 2> {$stderr}", $output, $exit_code); //when using exec2 the process hangs indefinitely sometimes
if ($exit_code!=0)
{
	trigger_error("grz-cli validate failed!\nSTDOUT:\n".implode("\n", file($stdout, FILE_IGNORE_NEW_LINES|FILE_SKIP_EMPTY_LINES))."\nSTDERR:\n".implode("\n", file($stderr, FILE_IGNORE_NEW_LINES|FILE_SKIP_EMPTY_LINES)), E_USER_ERROR);
}

print "validating submission took ".time_readable(microtime(true)-$time_start)."\n";
$time_start = microtime(true);

//Encrypt the submission
print "running grz-cli encrypt...\n";
$stdout = "{$folder}/logs/grz_cli_encrypt.stdout";
$stderr = "{$folder}/logs/grz_cli_encrypt.stderr";
$output = [];
$exit_code = null;
exec("{$grz_cli} encrypt --submission-dir {$folder} --config-file {$config} > {$stdout} 2> {$stderr}", $output, $exit_code); //when using exec2 the process hangs indefinitely sometimes
if ($exit_code!=0)
{
	print_r($output);
	trigger_error("grz-cli encrypt failed!\nSTDOUT:\n".implode("\n", file($stdout, FILE_IGNORE_NEW_LINES|FILE_SKIP_EMPTY_LINES))."\nSTDERR:\n".implode("\n", file($stderr, FILE_IGNORE_NEW_LINES|FILE_SKIP_EMPTY_LINES)), E_USER_ERROR);
}

print "encrypting submission took ".time_readable(microtime(true)-$time_start)."\n";
$time_start = microtime(true);

//Upload the submission
print "running grz-cli upload...\n";
$stdout = "{$folder}/logs/grz_cli_upload.stdout";
$stderr = "{$folder}/logs/grz_cli_upload.stderr";
$output = [];
$exit_code = null;
exec("{$grz_cli} upload --submission-dir {$folder} --config-file {$config} > {$stdout} 2> {$stderr}", $output, $exit_code); //when using exec2 the process hangs indefinitely sometimes
if ($exit_code!=0)
{
	print_r($output);
	trigger_error("grz-cli upload failed!\nSTDOUT:\n".implode("\n", file($stdout, FILE_IGNORE_NEW_LINES|FILE_SKIP_EMPTY_LINES))."\nSTDERR:\n".implode("\n", file($stderr, FILE_IGNORE_NEW_LINES|FILE_SKIP_EMPTY_LINES)), E_USER_ERROR);
}

//print submission ID (used by wrapper script when updating status in MVH db)
$submission_id = trim(implode("", file($stdout)));
print "SUBMISSION ID GRZ: {$submission_id}\n";

print "uploading submission took ".time_readable(microtime(true)-$time_start)."\n";
$time_start = microtime(true);

//if upload successfull, add 'Pruefbericht' to CM RedCap
if ($submission_type=='initial' && !$test)
{
	print "Adding Pruefbericht to CM RedCap...\n";
	add_submission_to_redcap($cm_id, "G", $tan_g);
}

//archive metadata JSON
copy($json_file, $mvh_folder."/metadata_archive/GRZ/{$cm_id}.json");

//clean up export folder if successfull
if (!$test)
{
	exec2("rm -rf {$folder} {$qc_folder}");
}

print "cleanup took ".time_readable(microtime(true)-$time_start)."\n";

?>
