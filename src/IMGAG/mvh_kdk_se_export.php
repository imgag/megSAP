<?php
/** 
	@page mvh_grz_export
*/

require_once("mvh_functions.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

function json_metadata($cm_data, $tan_k, $rc_data_json)
{
	//add entry section to suppress PHP warning below
	if (!isset($rc_data_json["entry"])) $rc_data_json["entry"] = [];
	
	$active_rcs = [];
	foreach($rc_data_json["entry"] as $entry)
	{
		$entry = $entry['resource'];

		//skip not active
		if ($entry['status']!='active') continue;
		
		//skip before Modellvorhaben
		$date = $entry['dateTime'];
		if ($date=="") continue;
		if (new DateTime($entry['dateTime']) < new DateTime("2024-07-01")) continue;
		
		$active_rcs[] = $entry;
	}
	if (count($active_rcs)==0) trigger_error("No active consent resources found!", E_USER_ERROR);
	if (count($active_rcs)>1) trigger_error("Several active consent resources found!", E_USER_ERROR);
	
	$te_date = xml_str($cm_data->datum_teilnahme);
	
	$output = [
			"type"=>"initial", //TODO implement other types
			"transferTAN" => $tan_k,
			"modelProjectConsent" => [
				"version" => "",
				"date" => $te_date,
				"version" => xml_str($cm_data->version_teilnahme),
				"provisions" => [
					[
						"date" => $te_date,
						"purpose" => "sequencing",
						"type" => xml_str($cm_data->particip_4)=="Ja" ? "permit" : "deny"
					],
					[
						"date" => $te_date,
						"purpose" => "case-identification",
						"type" => xml_str($cm_data->particip_4_1)=="Ja" ? "permit" : "deny"
					],
					[
						"date" => $te_date,
						"purpose" => "reidentification",
						"type" => xml_str($cm_data->particip_4_2)=="Ja" ? "permit" : "deny"
					]
				],
			],
			"researchConsents" => [
				$active_rcs[0]
				]
		];
			
	return $output;
}

function json_patient($info, $gl_data, $se_data)
{
	global $patient_id;
	
	$output = [
			"id" => $patient_id,
			"gender" => [ "code" => convert_gender($info['gender'])],
			"birthDate" => xml_str($se_data->birthdate),
			"healthInsurance" => [
				"type" => [
					"code"=> convert_coverage($gl_data->accounting_mode)
					]
				],
		    "address" => [
				"municipalityCode" => get_raw_value($se_data->psn, "ags"),
				],
			];
	return $output;
}

function json_patient_ref()
{
	global $patient_id;
	return [
			"id" => $patient_id,
			"type" => "Patient"
			];
}

function json_episode_of_care($se_data)
{
	$output = [
			"id"=>"ID_EOC_1",
			"patient" => json_patient_ref(),
			"period" => [
				"start" => xml_str($se_data->datum_kontakt_zse)
				]
			];
	
	return $output;
}

function json_diagnoses($se_data, $se_data_rep)
{
	//prepare codes
	$codes = [];
	$icd10 = xml_str($se_data->diag_icd10);
	$icd10_ver = xml_str($se_data->diag_icd10_ver);
	if ($icd10!="")
	{
		$codes[] = [
			"code" => get_raw_value($se_data->psn, "diag_icd10"),
			"display" => $icd10,
			"system" => "https://www.bfarm.de/DE/Kodiersysteme/Terminologien/ICD-10-GM",
			"version" => $icd10_ver
			
			
		];
	}
	$orpha = xml_str($se_data->diag_orphacode);
	$orpha_ver = xml_str($se_data->diag_orphacode_ver);
	if ($orpha!="")
	{
		$codes[] = [
			"code" => get_raw_value($se_data->psn, "diag_orphacode"),
			"display" => "ORPHA:".$orpha,
			"system" => "https://www.orpha.net",
			"version" => $orpha_ver
		];
	}
	if (xml_bool($se_data->orphacode_undiagnostiziert___1, false))
	{
		$codes[] = [
			"code" => "ORPHA:616874",
			"display" => "Fully investigated rare disorder without a determined diagnosis",
			"system" => "https://www.orpha.net",
		];
	}
	
	$alpha_id_se = xml_str($se_data->diag_se_code);
	if ($alpha_id_se!="")
	{
		$codes[] = [
			"code" => get_raw_value($se_data->psn, "diag_se_code"),
			"display" => $alpha_id_se,
			"system" => "https://www.bfarm.de/DE/Kodiersysteme/Terminologien/Alpha-ID-SE",
		];
	}
	if (count($codes)==0) trigger_error("No disease code found in case-management data!", E_USER_ERROR);
	
	//get onset date from HPO terms
	$onset_date = [];
	foreach($se_data_rep->item as $item)
	{		
		$onset = substr(xml_str($item->beginn_symptome), 0, 7);
		if ($onset=="") continue;

		$onset_date[] = $onset;
	}
	if (count($onset_date)==0) trigger_error("No HPO onset date found in SE meta data! It is needed for the diagnosis onset!", E_USER_ERROR);
	asort($onset_date);
	$onset_date = $onset_date[0];
	
	$output = [
			"id"=>"ID_DIAG_1",
			"patient" => json_patient_ref(),
			"recordedOn" => xml_str($se_data->datum_fallkonferenz), //TODO ok so?
			"onsetDate" => $onset_date,
			"familyControlLevel" => [
				"code" => convert_diag_recommendation(xml_str($se_data->diagnostik_empfehlung)),
				],
			"verificationStatus" => [
				"code" => convert_diag_status(xml_str($se_data->bewertung_gen_diagnostik)),
				],
			"codes" => $codes,
			//missing fields: notes
			];
	
	if (count($codes)!=3)
	{
		$output["missingCodeReason"] = [
        "code" => "no-matching-code",
        "display" => "Kein geeigneter Code (ICD-10-GM, ORDO, Alpha-ID-SE) verfügbar",
        "system" => "dnpm-dip/rd/diagnosis/missing-code-reason"
		];
	}
	
	return $output;
}

function json_hpos($se_data, $se_data_rep)
{
	global $db_ngsd;
	$output = [];
	
	$num = 1;
	foreach($se_data_rep->item as $item)
	{
		$hpo = xml_str($item->hpo);
		if ($hpo=="") continue;
		
		//convert HPO name to id
		$hpo_id = $db_ngsd->getValue("SELECT hpo_id FROM hpo_term WHERE name LIKE '{$hpo}'", "");  
		if ($hpo_id=="")
		{
			trigger_error("Could not convert HPO name '{$hpo}' to ID using the NGSD hpo_term table!", E_USER_WARNING);
			continue;
		}
		
		//create entry
		$entry = [
		"id" => "ID_HPO_{$num}",
		"patient" => json_patient_ref(),
		"recordedOn" => xml_str($se_data->datum_fallkonferenz), //TODO ok so?
		"value" => [
			"code" => $hpo_id,
			"display" => $hpo,
			"system" => "https://hpo.jax.org",
			] 
		//missing fields: status/history
		];
		
		//add optional stuff to entry
		$onset = substr(xml_str($item->beginn_symptome), 0, 7);
		if ($onset!="") $entry["onsetDate"] = $onset;
		
		$hpo_ver = xml_str($item->version_hpo);
		if ($hpo_ver!="") $entry["value"]["version"] = $hpo_ver;

		$output[] = $entry;
		++$num;
	}
	
	return $output;
}

function json_gmfcs($se_data_rep)
{
	$output = [];
	$num = 1;
	
	foreach($se_data_rep->item as $item)
	{
		$gmf = xml_str($item->gmfcs_wiedervorst_aend);
		if ($gmf=="") continue;
		
		$entry = [
			"id" => "ID_GMFCS_{$num}",
			"patient" => json_patient_ref(),
			"effectiveDate" => xml_str($item->datum_wiedervorst),
			"value" => [
				"code" => $gmf,
				"display" => "Level ".$gmf,
				"system" => "Gross-Motor-Function-Classification-System",			
				]
			];
		
		$output[] = $entry;
		++$num;
	}
	
	if(count($output)==0) return null;
	return $output;
}

function json_hospitalization($se_data)
{
	$output = [
		 "numberOfStays" => [
			"code" => convert_hospitalization_stays($se_data->anzahl_stat_behandlungen_en)
			],
		 "numberOfDays" => [
			"code" => convert_hospitalization_days($se_data->dauer_stat_vortherapie_en)
			],
		];
	
	return $output;
}

function json_care_plan($se_data, $se_data_rep) //'carePlan' is misleading. This object captures the decisions of the second "Fallkonferenz".
{	
	//prepare therapy recommendations
	$therapy_recoms = [];
	$num = 1;
	foreach($se_data_rep->item as $item)
	{
		$therapy_type = xml_str($item->therapie_art);
		if ($therapy_type=="") continue;
		
		//create entry
		$entry = [
		"id" => "ID_THERAPY_{$num}",
		"patient" => json_patient_ref(),
		"issuedOn" => xml_str($se_data->klin_datum_fallkonferenz), //TODO ok so?
		"category" => [
			"code" => convert_therapy_category($therapy_type),
			],
		"type" => [
			"code" => convert_therapy_type($item->therapie_beschreibung),
			],
		//TODO supportingVariants
		//missing fields: meciation (not in SE data)
		];

		$therapy_recoms[] = $entry;
		++$num;
	}
	
	//prepare study enrollment recommendations
	$study_recoms = [];
	$num = 1;
	foreach($se_data_rep->item as $item)
	{
		$study_register = xml_str($item->studien_register);
		if ($study_register=="") continue;
		
		//create entry
		$entry = [
		"id" => "ID_STUDY_{$num}",
		"patient" => json_patient_ref(),
		"issuedOn" => xml_str($se_data->klin_datum_fallkonferenz), //TODO ok so?
		"study" => [
				[
					"id" => xml_str($item->studien_id),
					"system" => convert_study_register($study_register),
					"type" => "Study",
					"display" => xml_str($item->studienname)
				]
			],
		//TODO supportingVariants
		];

		$study_recoms[] = $entry;
		++$num;
	}
	
	//prepare clinical management recommendations
	$clin_recoms = [];
	$num = 1;
	foreach($se_data_rep->item as $item)
	{
		$clinical_mgmt = xml_str($item->klinisches_management_dr);
		if ($clinical_mgmt=="") continue;
		
		//create entry
		$entry = [
		"id" => "ID_THERAPY_{$num}",
		"patient" => json_patient_ref(),
		"issuedOn" => xml_str($se_data->klin_datum_fallkonferenz), //TODO ok so?
		"type" => [
			"code" => convert_clinincal_management($clinical_mgmt),
			]
		];

		$clin_recoms[] = $entry;
		++$num;
	}
	if (count($clin_recoms)>1) trigger_error("More than one clinical management recommendation found! Only one is allowed in KDK-SE!", E_USER_ERROR);
	
	$output = [
			"id"=>"ID_CARE_PLAN_1",
			"patient" => json_patient_ref(),
			"issuedOn" => xml_str($se_data->klin_datum_fallkonferenz), //TODO ok so?
			"therapyRecommendations" => $therapy_recoms,
			"studyEnrollmentRecommendations" => $study_recoms,
			];
	
	//optional stuff
	$recom_counceling = xml_bool($se_data->empf_hg_beratung, true);
	if (!is_null($recom_counceling))
	{
		$output["geneticCounselingRecommended"] = $recom_counceling;
	}
	
	$recom_reeval = xml_bool($se_data->empf_reeval, true);
	if (!is_null($recom_reeval))
	{
		$output["reevaluationRecommended"] = $recom_reeval;
	}
	
	if (count($clin_recoms)==1)
	{
		$output["clinicalManagementRecommendation"] = $clin_recoms[0];
	}
	
	return $output;
}

//TODO: add support for cases without sequencing (see SE RedCAP aufnahme_mvh/fallkonferenz_grund)
function json_ngs_report($cm_data, $se_data, $info)
{
	global $is_lrgs;
	global $db_ngsd;
	
	$output = [
			"id" => "ID_NGS_REPORT_1",
			"patient" => json_patient_ref(),
			"issuedOn" =>  xml_str($cm_data->gen_finding_date),
			"type" => [
				"code" => ($is_lrgs ? "genome-long-read" : "genome-short-read"),
				],
			"sequencingInfo" =>[
				"platform" => [
					"code" => convert_sequencing_platform($info["device_type"])
					],
				"kit" => $info['sys_name']
				],
			"result" => [
				
				]
		];
	//missing: autozygosity (there is no definition how to calculate it), conclusion (makes no sense at all)
	return $output;
}

//parse command line arguments
$parser = new ToolBase("mvh_kdk_se_export", "KDK-SE export for Modellvorhaben.");
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

//get patient identifer (pseudonym from case management) - this is the ID that is used to identify submissions from the same case by GRZ/KDK
$cm_data = get_cm_data($db_mvh, $case_id);
$patient_id = xml_str($cm_data->psn);
if ($patient_id=="") trigger_error("No patient identifier set for sample '{$ps}'!", E_USER_ERROR);

//create export folder
print "MVH DB id: {$case_id} (CM ID: {$cm_id} / CM pseudonym: {$patient_id})\n";
$folder = realpath($mvh_folder)."/kdk_se_export/{$case_id}/";
if ($clear) exec2("rm -rf {$folder}");
exec2("mkdir -p {$folder}/metadata/");
print "export folder: {$folder}\n";

//determine tanG==VNg
$sub_ids = $db_mvh->getValues("SELECT id FROM `submission_kdk_se` WHERE status='pending' AND case_id='{$case_id}'");
if (count($sub_ids)!=1)  trigger_error(count($sub_ids)." pending KDK-SE submissions for case {$case_id}. Must be one!", E_USER_ERROR);
$sub_id = $sub_ids[0];
print "ID in submission_kdk_se table: {$sub_id}\n";
$tan_k = $db_mvh->getValue("SELECT tank FROM submission_kdk_se WHERE id='{$sub_id}'");
print "TAN: {$tan_k}\n";

//get data
$ps = $db_mvh->getValue("SELECT ps FROM case_data WHERE id='{$case_id}'");
print "index sample: {$ps}\n";
$info = get_processed_sample_info($db_ngsd, $ps);
$is_lrgs = $info['sys_type']=="lrGS";
$gl_data = get_gl_data($db_mvh, $case_id);
$se_data = get_se_data($db_mvh, $case_id);
$se_data_rep = get_se_data($db_mvh, $case_id, true);
$rc_data_json = get_rc_data_json($db_mvh, $case_id);

//create base JSON
print "\n";
print "### creating JSON ###\n";

$json = [
	"metadata" => json_metadata($cm_data, $tan_k, $rc_data_json),
	"patient" => json_patient($info, $gl_data, $se_data),
	"episodesOfCare" => [ json_episode_of_care($se_data) ],
	"diagnoses" => [ json_diagnoses($se_data, $se_data_rep) ],
	"hpoTerms" => json_hpos($se_data, $se_data_rep),
	"hospitalization" => json_hospitalization($se_data),
	"carePlans" => [ json_care_plan($se_data, $se_data_rep) ],
	"ngsReports" => [ json_ngs_report($cm_data, $se_data, $info)],
	//TODO: followUps, therapies
	];
	
//add optional parts to JSON
$gmfcs = json_gmfcs($se_data_rep);
if (!is_null($gmfcs)) $json["gmfcsStatus"] = $gmfcs;

//write JSON
$json_file = "{$folder}/metadata/metadata.json";
file_put_contents($json_file, json_encode($json, JSON_PRETTY_PRINT));

//validate JSON
print "\n";
print "### validating JSON ###\n";

$url = "https://".($test ? "preview.dnpm-dip.net" : "dnpm-dip-p.med.uni-tuebingen.de")."/api/rd/etl/patient-record:validate";
$ch = curl_init();
curl_setopt($ch, CURLOPT_POST, true);
curl_setopt($ch, CURLOPT_URL, $url);
curl_setopt($ch, CURLOPT_POSTFIELDS, file_get_contents($json_file));
curl_setopt($ch, CURLOPT_RETURNTRANSFER, true);
$result = curl_exec($ch);
curl_close($ch);

//parse/show validation result
print "exit code: ".curl_getinfo($ch, CURLINFO_HTTP_CODE)."\n";
print "\n";
$validation_error = false;
foreach(json_decode($result)->issues as $issue)
{
	if ($issue->message=="Fehlende Angabe 'Krankenkassen-IK'") continue;
	
	$type = strtoupper($issue->severity);
	if ($type=="ERROR") $validation_error = true;
	print "{$type}: ".$issue->message."\n";
	print "  PATH: ".$issue->path."\n";
	print "\n";
}
if ($validation_error) trigger_error("Validation failed!", E_USER_ERROR);

//upload JSON
print "\n";
print "### uploading JSON ###\n";

$url = "https://".($test ? "preview.dnpm-dip.net" : "dnpm-dip-p.med.uni-tuebingen.de")."/api/rd/etl/patient-record";
$ch = curl_init();
curl_setopt($ch, CURLOPT_POST, true);
curl_setopt($ch, CURLOPT_URL, $url);
curl_setopt($ch, CURLOPT_POSTFIELDS, file_get_contents($json_file));
curl_setopt($ch, CURLOPT_RETURNTRANSFER, true);
$result = curl_exec($ch);
curl_close($ch);


//parse/show upload result
$exit_code = curl_getinfo($ch, CURLINFO_HTTP_CODE);
print "exit code: {$exit_code}\n";
print "\n";
foreach(json_decode($result)->issues as $issue)
{
	if ($issue->message=="Fehlende Angabe 'Krankenkassen-IK'") continue;
	
	$type = strtoupper($issue->severity);
	if ($type=="ERROR") $validation_error = true;
	print "{$type}: ".$issue->message."\n";
	print "  PATH: ".$issue->path."\n";
	print "\n";
}
if ($exit_code!="201") trigger_error("Upload failed!", E_USER_ERROR);



?>