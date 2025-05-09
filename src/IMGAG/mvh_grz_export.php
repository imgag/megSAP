<?php
/** 
	@page mvh_grz_export
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//convert NGSD gender to JSON schema
function convert_gender($gender)
{
	if ($gender=="n/a") $gender = "unknown";
	return $gender;
}

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


//parse command line arguments
$parser = new ToolBase("mvh_grz_export", "GRZ export for Modellvorhaben.");
$parser->addInt("case_id", "'id' in 'data_data' of 'MVH' database.", false);
extract($parser->parse($argv));

//init
$db_mvh = DB::getInstance("MVH");
$db_ngsd = DB::getInstance("NGSD");
$mvh_folder = get_path("mvh_folder");

//check that case ID is valid
$id = $db_mvh->getValue("SELECT id FROM case_data WHERE id='{$case_id}'");
if ($id=="") trigger_error("No case with id '{$case}' in MVH database!", E_USER_ERROR);

//parse case management data
$cm_data = get_cm_data($case_id);

//check network > determine KDK and study_subtype
$network = $cm_data->network_title;
$kdk = "";
$study_subtype = "";
$$disease_type = "";
if ($network="Netzwerk Seltene Erkrankungen")
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
}
else if ($network=="Deutsches Konsortium Familiärer Brust- und Eierstockkrebs")
{
	$kdk = "KDKL00003"; //DK-FBREK - Leipzig
	$study_subtype = "tumor+germline";
	$disease_type = "oncological";
}
else trigger_error("Unhandled network type '{$network}'!", E_USER_ERROR); 

//check seqencing mode
$seq_mode = $cm_data->seq_mode;
if ($seq_mode!="WGS" && $seq_mode!="WES") trigger_error("Unhandled seq_mode '{$seq_mode}'!", E_USER_ERROR);

print "case: {$case_id} (seq_mode: {$seq_mode} / network: {$network})\n";

//check germline processed sample is ok
$ps = $db_mvh->getValue("SELECT ps FROM case_data WHERE id='{$case_id}'");
$info = get_processed_sample_info($db_ngsd, $ps);
$sys = $info['sys_name_short'];
$patient_id = $info['patient_identifier'];
if ($patient_id=="") trigger_error("No patient identifier set for sample '{$ps}''!", E_USER_ERROR);
print "germline sample: {$ps} (system: {$sys})\n";

//check germline processed sample is ok
$ps_t = "";
if ($study_subtype=="tumor+germline")
{
	$ps_t = $db_mvh->getValue("SELECT ps_t FROM case_data WHERE id='{$case_id}'");
	if ($ps_t=="") trigger_error("Could no tumor sample for case '{$case_id}' in network '{$network}'!", E_USER_ERROR);
	
	//check tumor and germline have same system
	$info_t = get_processed_sample_info($db_ngsd, $ps_t);
	$sys_t = $info_t['sys_name_short'];
	if ($sys!=$sys_t) trigger_error("Mismatching tumor/normal processing systems for case '{$case_id}' in network '{$network}': '{$sys}' vs '{$sys_t}'", E_USER_ERROR);
}

//determine ROI
$roi = "";
if ($sys=="twistCustomExomeV2" || $sys=="twistCustomExomeV2Covaris") $roi = "{$mvh_folder}/rois/twist_exome_core_plus_refseq.bed";
if ($sys=="hpHBOCv5") $roi = "/mnt/storage2/megSAP/data/enrichment/hpHBOCv5_2014_10_27.bed"; //TODO: only for testing > remove!
if ($seq_mode!="WGS" && $roi=="") trigger_error("Could not determine target region for sample '{$ps}' with processing system '{$sys}'!", E_USER_ERROR);

//determine tanG==VNg
$sub_ids = $db_mvh->getValues("SELECT id FROM `submission_grz` WHERE status='pending' AND case_id='$case_id}'");
if (count($sub_ids)!=1)  trigger_error(count($sub_ids)." pending GRZ submissions for case {$case_id}. Must be one!", E_USER_ERROR);
$sub_id = $sub_ids[0];
$tang = $db_mvh->getValue("SELECT tang FROM submission_grz WHERE id='{$sub_id}'");

//create export folder
$folder = "{$mvh_folder}/export_grz/{$case_id}/";
exec2("mkdir -p {$folder}");
print "export folder: {$folder}\n";
exec2("mkdir -p {$folder}/files/");
exec2("mkdir -p {$folder}/metadata/");

//TODO also handle tumor/normal
//create raw data (FASTQs + target region)
$bam = $info['ps_bam'];
$read_length = get_processed_sample_qc($db_ngsd, $ps, "QC:2000006");
while (contains($read_length, "-")) $read_length = explode("-", $read_length, 2)[1];
//TODO $parser->execApptainer("ngs-bits", "BamToFastq", "-in {$bam} -out1 {$folder}/files/{$tang}.read1.fastq.gz -out2 {$folder}/files/{$tang}.read2.fastq.gz -extend {$read_length}", [$bam], ["{$folder}/files/"]);
if ($roi!="") exec2("cp {$roi} {$folder}/files/target_regions.bed");

//create meta data
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
			]
		]
		//TODO researchConsents
	];
print_r(json_encode($json, JSON_PRETTY_PRINT));


//Validate the submission
#grz-cli validate --submission-dir EXAMPLE_SUBMISSION

//Encrypt the submission
#grz-cli encrypt --submission-dir EXAMPLE_SUBMISSION

//Upload the submission
#grz-cli upload --submission-dir EXAMPLE_SUBMISSION


/*
#TODO:
- test if export works with other user, e.g. bioinf
- store schema of MVH DB


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


