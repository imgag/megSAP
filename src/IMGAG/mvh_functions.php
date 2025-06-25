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
	$xml_obj = simplexml_load_string($xml_string, "SimpleXMLElement", LIBXML_NOCDATA);
	return $xml_obj;
}

################# XML functions #################

function xml_str($element)
{
	return trim((string)$element);
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
	$xml = simplexml_load_string(curl_exec($ch), "SimpleXMLElement", LIBXML_NOCDATA);;
	curl_close($ch);
	
	return xml_str($xml->item->$field);
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
	if ($accounting_mode=="Modellvorhaben") $converage_type = "GKV";
	else if ($accounting_mode=="Privat (GOÄ)") $converage_type = "PKV";
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
	if ($name=="unsolved") return "unconfirmed";
	if ($name=="unclear" || $name=="3") return "provisional";
	if ($name=="solved") return "confirmed";
	if ($name=="partially-solved") return "partial";
	
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
	//TODO fix when changed (see Email from 25.06.25: MVH Hospitalization)
	if ($value=="bis zu 15" || $value=="bis zu 10") return "up-to-fifteen";
	if ($value==">15") return "up-to-fifty";
	
	trigger_error(__FUNCTION__.": Unhandled value '{$value}'!", E_USER_ERROR);
}

?>
