<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

################# data from different sources #################

//returns case management data as an XML object
function get_cm_data($db_mvh, $case_id)
{
	$main_item = null;
	
	$xml_string = $db_mvh->getValue("SELECT cm_data FROM case_data WHERE id='{$case_id}'", "");
	$xml_obj = simplexml_load_string($xml_string, "SimpleXMLElement", LIBXML_NOCDATA);
	foreach($xml_obj->item as $item)
	{
		if ($item->redcap_repeat_instance!="") continue;
		
		$main_item = $item;
	}
	
	return $main_item;
}

//returns SE data as an XML object
function get_se_data($db_mvh, $case_id, $with_repeat_instances=false)
{	
	$xml_string = $db_mvh->getValue("SELECT se_data FROM case_data WHERE id='{$case_id}'", "");
	$xml_obj = simplexml_load_string($xml_string, "SimpleXMLElement", LIBXML_NOCDATA);
	
	if ($with_repeat_instances) return $xml_obj;
	
	//determine main item
	$main_item = null;
	foreach($xml_obj->item as $item)
	{
		if ($item->redcap_repeat_instance!="") continue;
		
		$main_item = $item;
	}
	return $main_item;
}
//returns GenLab data as an XML object
function get_gl_data($db_mvh, $case_id)
{
	$xml_string = $db_mvh->getValue("SELECT gl_data FROM case_data WHERE id='{$case_id}'", "");
	$xml_obj = simplexml_load_string($xml_string, "SimpleXMLElement", LIBXML_NOCDATA);
	return $xml_obj;
}


//returns research consent data as an XML object
function get_rc_data($db_mvh, $case_id)
{
	$xml_string = $db_mvh->getValue("SELECT rc_data FROM case_data WHERE id='{$case_id}'", "");
	return simplexml_load_string($xml_string, "SimpleXMLElement", LIBXML_NOCDATA);
}

//returns research consent data as an JSON object
function get_rc_data_json($db_mvh, $case_id)
{
	$json_string = $db_mvh->getValue("SELECT rc_data_json FROM case_data WHERE id='{$case_id}'", "");
	return json_decode($json_string, true);
}

################# XML functions #################

function xml_str($value)
{
	return trim((string)$value);
}

function xml_bool($value, $allow_unset)
{
	$value = strtolower(xml_str($value));
	if ($value=="yes" || $value=="checked") return true;
	if ($value=="no" || $value=="unchecked") return false;
	if ($allow_unset && $value=="") return null;
	
	trigger_error(__FUNCTION__.": Unhandled value '{$value}'!", E_USER_ERROR);
}

//returns the raw value for a drop-down file (our man query returns the label only)
function get_raw_value($record_id, $field)
{
	$data = array(
		'token' => get_path("mvh_redcap_token_se"),
		'content' => 'record',
		'action' => 'export',
		'format' => 'xml',
		'type' => 'flat',
		'csvDelimiter' => '',
		'records' => [xml_str($record_id)],
		'fields' => [$field],
		'rawOrLabel' => 'raw',
		'rawOrLabelHeaders' => 'label',
		'exportCheckboxLabel' => 'false',
		'exportSurveyFields' => 'false',
		'exportDataAccessGroups' => 'false',
		'returnFormat' => 'csv'
	);
	$ch = curl_init();
	curl_setopt($ch, CURLOPT_URL, 'https://redcap.extern.medizin.uni-tuebingen.de/api/');
	curl_setopt($ch, CURLOPT_RETURNTRANSFER, true);
	curl_setopt($ch, CURLOPT_SSL_VERIFYPEER, false);
	curl_setopt($ch, CURLOPT_VERBOSE, 0);
	curl_setopt($ch, CURLOPT_FOLLOWLOCATION, true);
	curl_setopt($ch, CURLOPT_AUTOREFERER, true);
	curl_setopt($ch, CURLOPT_MAXREDIRS, 10);
	curl_setopt($ch, CURLOPT_CUSTOMREQUEST, 'POST');
	curl_setopt($ch, CURLOPT_FRESH_CONNECT, 1);
	curl_setopt($ch, CURLOPT_POSTFIELDS, http_build_query($data, '', '&'));
	curl_setopt($ch, CURLOPT_PROXY, '');
	$result = curl_exec($ch);
	curl_close($ch);
	if ($result===false) trigger_error('CURL ERROR: '.curl_error($ch), E_USER_ERROR);
	$exit_code = curl_getinfo($ch, CURLINFO_HTTP_CODE);
	if ($exit_code!="200") trigger_error("get_raw_value: Invalid exit code '{$exit_code}'", E_USER_ERROR);
	
	$xml = simplexml_load_string($result, "SimpleXMLElement", LIBXML_NOCDATA);
	return xml_str($xml->item->$field);
}

//returns the raw value for a drop-down file (our man query returns the label only)
function add_submission_to_redcap($record_id, $data_type, $tan)
{
	//input checks
	if ($data_type!="G" && $data_type!="K") trigger_error("Invalid type '{$data_type}'", E_USER_ERROR);
	
	//TODO implement other types (report_type)
	$xml = "<records>
				<item>
					<record_id><![CDATA[{$record_id}]]></record_id>
					<redcap_repeat_instrument><![CDATA[pruefbericht]]></redcap_repeat_instrument>
					<redcap_repeat_instance><![CDATA[new]]></redcap_repeat_instance>
					<report_t_vn><![CDATA[{$tan}]]></report_t_vn>
					<report_data_type><![CDATA[{$data_type}]]></report_data_type>
					<report_type><![CDATA[0]]></report_type>
					<report_date><![CDATA[".date("Y-m-d")."]]></report_date>
				</item>
			</records>";

	$data = array(
		'token' => get_path("mvh_redcap_token_cm"),
		'content' => 'record',
		'format' => 'xml',
		'type' => 'flat',
		'data' => $xml,
	);
	$ch = curl_init();
	curl_setopt($ch, CURLOPT_URL, 'https://redcap.extern.medizin.uni-tuebingen.de/api/');
	curl_setopt($ch, CURLOPT_RETURNTRANSFER, true);
	curl_setopt($ch, CURLOPT_SSL_VERIFYPEER, false);
	curl_setopt($ch, CURLOPT_VERBOSE, 0);
	curl_setopt($ch, CURLOPT_FOLLOWLOCATION, true);
	curl_setopt($ch, CURLOPT_AUTOREFERER, true);
	curl_setopt($ch, CURLOPT_MAXREDIRS, 10);
	curl_setopt($ch, CURLOPT_CUSTOMREQUEST, 'POST');
	curl_setopt($ch, CURLOPT_FRESH_CONNECT, 1);
	curl_setopt($ch, CURLOPT_POSTFIELDS, http_build_query($data, '', '&'));
	curl_setopt($ch, CURLOPT_PROXY, '');
	$result = curl_exec($ch);
	curl_close($ch);
	if ($result===false) trigger_error('CURL ERROR: '.curl_error($ch), E_USER_ERROR);
	$exit_code = curl_getinfo($ch, CURLINFO_HTTP_CODE);
	if ($exit_code!="200")
	{
		trigger_error("get_raw_value: Invalid exit code '{$exit_code}' - result: {$result}", E_USER_ERROR);
	}
}

################# converter functions #################

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

function convert_coverage($accounting_mode)
{
	if ($accounting_mode=="gesetzliche Krankenversicherung") $converage_type = "GKV";
	else if ($accounting_mode=="private Krankenversicherung") $converage_type = "PKV";
	else if ($accounting_mode=="Berufsgenossenschaft") $converage_type = "BG";
	else if ($accounting_mode=="Selbstzahler") $converage_type = "SEL";
	else if ($accounting_mode=="Sozialamt") $converage_type = "SOZ";
	else if ($accounting_mode=="gesetzliche Pflegeversicherung") $converage_type = "GPV";
	else if ($accounting_mode=="private Pflegeversicherung") $converage_type = "PPV";
	else if ($accounting_mode=="Beihilfe") $converage_type = "Beihilfe";
	else if ($accounting_mode=="sonstiger Kostenträger") $converage_type = "SKT";
	else if ($accounting_mode=="unknown") $converage_type = "UNK";
	else trigger_error("Could not determine coverage type from GenLab accounting mode '{$accounting_mode}'!", E_USER_ERROR);
	
	return $converage_type;
}

function convert_tissue($tissue)
{
	if ($tissue=="blood") return "BTO:0000089";
	if ($tissue=="buccal mucosa") return "BTO:0003833";
	if ($tissue=="fibroblast") return "BTO:0000452";
	if ($tissue=="lymphocyte") return "BTO:0000775";
	if ($tissue=="skin") return "BTO:0001253";
	if ($tissue=="muscle") return "BTO:0000887";
		
	trigger_error(__FUNCTION__.": Unhandled tissue '{$tissue}'!", E_USER_ERROR);
}

function convert_system_type($type)
{
	if ($type=="WES") return "wes";
	if ($type=="WGS") return "wgs";
	if ($type=="lrGS") return "wgs_lr";
		
	trigger_error(__FUNCTION__.": Unhandled type '{$type}'!", E_USER_ERROR);
}

function convert_kit_manufacturer($name, $allow_ont=true)
{
	if ($name=="TruSeqPCRfree") return "Illumina";
	if ($name=="twistCustomExomeV2") return "Twist";
	if ($name=="twistCustomExomeV2Covaris") return "Twist";
	if ($name=="LR-ONT-SQK-LSK114") return $allow_ont ? "ONT" : "other";
	
	trigger_error(__FUNCTION__.": Unhandled name '{$name}'!", E_USER_ERROR);
}

function convert_sequencer_manufacturer($name)
{
	if ($name=="NovaSeq6000" || $name=="NovaSeqXPlus") return "Illumina";
	if ($name=="PromethION") return "ONT";
	
	trigger_error(__FUNCTION__.": Unhandled name '{$name}'!", E_USER_ERROR);
}

function convert_diag_recommendation($name)
{
	if ($name=="Einzelgenom") return "single-genome";
	if ($name=="Duogenom") return "duo-genome";
	if ($name=="Triogenom") return "trio-genome";
	
	trigger_error(__FUNCTION__.": Unhandled name '{$name}'!", E_USER_ERROR);
}

function convert_diag_status($name)
{
	if ($name=="keine genetische Diagnosestellung") return "unconfirmed";
	if ($name=="Genetische Verdachtsdiagnose" || $name=="3") return "provisional";
	if ($name=="Genetische Diagnose gesichert") return "confirmed";
	if ($name=="klinischer Phänotyp nur partiell gelöst") return "partial";
	
	trigger_error(__FUNCTION__.": Unhandled name '{$name}'!", E_USER_ERROR);
}

function convert_hospitalization_stays($value)
{
	$value = xml_str($value);
	
	if ($value=="" || $value=="unbekannt") return "unknown";
	if ($value=="keine") return "none";
	if ($value=="bis zu 5") return "up-to-five";
	if ($value=="bis zu 10") return "up-to-ten";
	if ($value=="bis zu 15") return "up-to-fifteen";
	if ($value==">15") return "over-fifteen";
	
	trigger_error(__FUNCTION__.": Unhandled value '{$value}'!", E_USER_ERROR);
}
	
function convert_hospitalization_days($value)
{
	$value = xml_str($value);
	
	if ($value=="" || $value=="unbekannt") return "unknown";
	if ($value=="keine") return "none";
	if ($value=="bis zu 5") return "up-to-five";
	if ($value=="bis zu 15") return "up-to-fifteen";
	if ($value=="bis zu 50") return "up-to-fifty";
	if ($value==">50") return "over-fifty";
	
	trigger_error(__FUNCTION__.": Unhandled value '{$value}'!", E_USER_ERROR);
}

function convert_therapy_category($value)
{	
	if ($value=="Symptomatische nicht-medikamentöse Therapie (Fördermaßnahmen)") return "symptomatic";
	if ($value=="Symptomatische medikamentöse Therapie (bspw. antispastische Medikation)") return "symptomatic";
	if ($value=="Symptomatische interventionelle Therapie (Operationen, Injektionen)") return "symptomatic";
	if ($value=="Kausale Therapie (medikamentös)") return "causal";
	if ($value=="Kausale Therapie (interventionell)") return "causal";
	
	trigger_error(__FUNCTION__.": Unhandled value '{$value}'!", E_USER_ERROR);
}

function convert_therapy_type($value)
{	
	if ($value=="medikamentös_systemisch") return "systemic-medication";
	if ($value=="medikamentös_zielgerichtet") return "targeted-medication";
	if ($value=="medikamentös_Prävention") return "prevention-medication";
	if ($value=="Gentherapie") return "genetic";
	if ($value=="Prophylaxe") return "prophylactic";
	if ($value=="Früherkennung") return "early-detection";
	if ($value=="Kombinationstherapie") return "combination";
	if ($value=="Ernährung") return "nutrition";
	if ($value=="andere") return "other";
	
	trigger_error(__FUNCTION__.": Unhandled value '{$value}'!", E_USER_ERROR);
}

function convert_study_register($value)
{	
	if ($value=="NCT") return "NCT";
	if ($value=="DRKS") return "DRKS";
	if ($value=="Eudra-CT") return "Eudra-CT";
	if ($value=="other") return "other";
	
	trigger_error(__FUNCTION__.": Unhandled value '{$value}'!", E_USER_ERROR);
}

function convert_clinincal_management($value)
{	
	if ($value=="Indikatorerkrankungsspezifische Ambulanz") return "disease-specific-ambulatory-care";
	if ($value=="Andere Hochschulambulanz") return "university-ambulatory-care";
	if ($value=="Eigenes ZSE") return "local-crd";
	if ($value=="Anderes ZSE") return "other-crd";
	if ($value=="Andere Ambulanz") return "other-ambulatory-care";
	if ($value=="Hausarzt") return "gp";
	if ($value=="Niedergelassener Facharzt") return "specialist";
	
	trigger_error(__FUNCTION__.": Unhandled value '{$value}'!", E_USER_ERROR);
}

function convert_sequencing_platform($name)
{
	if ($name=="NovaSeq6000" || $name=="NovaSeqXPlus") return "illu";
	if ($name=="PromethION") return "ont";
	
	trigger_error(__FUNCTION__.": Unhandled name '{$name}'!", E_USER_ERROR);
}

function convert_outcome($name)
{
	if ($name=="no significant findings" || $name=="significant findings - non-genetic" || $name=="candidate gene" || $name=="significant findings - second method") return "no-pathogenic-variant-detected";
	
	return "";
}

function convert_noseq_reason($name)
{
	if ($name=="Zieldiagnostik empfohlen") return "targeted-diagnostics-recommended";
	if ($name=="wahrscheinlich psychosomatische Erkrankung") return "psychosomatic";
	if ($name=="wahrscheinlich häufige Erkrankung") return "not-rare-disease";
	if ($name=="wahrscheinlich nicht genetische Ursache") return "non-genetic-cause";
	if ($name=="anderer Grund") return "other";
	
	trigger_error(__FUNCTION__.": Unhandled name '{$name}'!", E_USER_ERROR);
}


?>
