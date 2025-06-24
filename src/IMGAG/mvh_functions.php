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
		if ($item->redcap_repeat_instance=="1") continue;
		
		$main_item = $item;
	}
	
	return $main_item;
}

//returns SE data as an XML object
function get_se_data($db_mvh, $case_id)
{
	$main_item = null;
	
	$xml_string = $db_mvh->getValue("SELECT se_data FROM case_data WHERE id='{$case_id}'", "");
	$xml_obj = simplexml_load_string($xml_string, "SimpleXMLElement", LIBXML_NOCDATA);
	foreach($xml_obj->item as $item)
	{
		if ($item->redcap_repeat_instance=="1") continue;
		
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
		
	trigger_error("Could not convert tissue '{$tissue}' to BTO id!", E_USER_ERROR);
}

function convert_system_type($type)
{
	if ($type=="WES") return "wes";
	if ($type=="WGS") return "wgs";
	if ($type=="lrGS") return "wgs_lr";
		
	trigger_error("Unhandled processing system type '{$type}'!", E_USER_ERROR);
}

function convert_kit_manufacturer($name, $allow_ont=true)
{
	if ($name=="TruSeqPCRfree") return "Illumina";
	if ($name=="twistCustomExomeV2") return "Twist";
	if ($name=="twistCustomExomeV2Covaris") return "Twist";
	if ($name=="LR-ONT-SQK-LSK114") return $allow_ont ? "ONT" : "other";
	
	trigger_error("Unhandled processing system name '{$name}'!", E_USER_ERROR);
}

function convert_sequencer_manufacturer($name)
{
	if ($name=="NovaSeq6000" || $name=="NovaSeqXPlus") return "Illumina";
	if ($name=="PromethION") return "ONT";
	
	trigger_error("Unhandled sequencer name '{$name}'!", E_USER_ERROR);
}

?>
