<?php
/** 
	@page mvh_grz_export
*/

require_once("mvh_functions.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

function json_patient($info, $gl_data, $se_data)
{
	$output = [
			"id"=>"ID_PAT_1",
			"gender"=>[ "code" => convert_gender($info['gender'])],
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

function json_episode_of_care($se_data)
{
	$output = [
			"id"=>"ID_EOC_1",
			"patient" => [
				"id" => "ID_PAT_1",
				"type" => "Patient"
				],
			"period" => [
				"start" => xml_str($se_data->datum_kontakt_zse)
				]
			];
	
	return $output;
}

function json_diagnoses($se_data)
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
			"display" => $orpha,
			"system" => "https://www.orpha.net",
			"version" => $icd10_ver
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
	
	$output = [
			"id"=>"ID_DIAG_1",
			"patient" => [
				"id" => "ID_PAT_1",
				"type" => "Patient"
				],
			"recordedOn" => xml_str($se_data->datum_fallkonferenz), //TODO ok so?
			"onsetDate" => "TODO", //TODO min(HPO onset)!?
			"familyControlLevel" => [
				"code" => convert_diag_recommendation(xml_str($se_data->diagnostik_empfehlung)),
				],
			"verificationStatus" => [
				"code" => convert_diag_status(xml_str($se_data->bewertung_gen_diagnostik)),
				],
			"codes" => $codes,
			//TODO? missingCodeReason, notes
			];
	
	return $output;
}

//parse command line arguments
$parser = new ToolBase("mvh_kdk_se_export", "KDK-SE export for Modellvorhaben.");
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

//create export folder
$folder = realpath($mvh_folder)."/kdk_se_export/{$case_id}/";
if ($clear) exec2("rm -rf {$folder}");
exec2("mkdir -p {$folder}/metadata/");
print "export folder: {$folder}\n";


//get data
$ps = $db_mvh->getValue("SELECT ps FROM case_data WHERE id='{$case_id}'");
$info = get_processed_sample_info($db_ngsd, $ps);
$gl_data = get_gl_data($db_mvh, $case_id);
$se_data = get_se_data($db_mvh, $case_id);

//create JSON
$json = [
	"patient" => json_patient($info, $gl_data, $se_data),
	"episodesOfCare" => [ json_episode_of_care($se_data) ],
	"diagnoses" => [ json_diagnoses($se_data) ]
	];

//write JSON
file_put_contents("{$folder}/metadata/metadata.json", json_encode($json, JSON_PRETTY_PRINT));
print implode("", file("{$folder}/metadata/metadata.json"));

?>