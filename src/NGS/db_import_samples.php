<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("db_import_samples", "Batch-imports samples to the NGSD.");
$parser->addInfile("in",  "Input sample list in TSV format.", false);
extract($parser->parse($argv));

function get_sender_id($sender)
{
	$db_connect = DB::getInstance("NGSD");
	$result_array = $db_connect->executeQuery("SELECT id FROM sender WHERE name ='".$sender."'");
	if (count($result_array==1))
	{
		return $result_array[0]['id'];
	}
	else
	{
		return false;
	}
}

function get_receiver_id($receiver)
{
	$db_connect = DB::getInstance("NGSD");
	$result_array = $db_connect->executeQuery("SELECT id FROM user WHERE name ='".$receiver."'");
	if (count($result_array==1))
	{
		return $result_array[0]['id'];
	}
	else
	{
		return false;
	}
}

function get_species_id($species)
{
	$db_connect = DB::getInstance("NGSD");
	$result_array = $db_connect->executeQuery("SELECT id FROM species WHERE name ='".$species."'");
	if (count($result_array==1))
	{
		return $result_array[0]['id'];
	}
	else
	{
		return false;
	}
}

function replace_empty_with_null($input)
{
	if ($input=="") return null;
	return $input;
}

function replace_empty_with_na($input)
{
	if ($input=="") return 'n/a';
	return $input;
}

function convert_yes_and_no_to_1_and_0($input)
{
	if ($input=="yes") return 1;
	else if ($input=="no") return 0;
	return $input;
	
}

function add_sample($name,$external_name,$sender,$received,$received_by,$sample_type,$tumor,$ffpe,$species,$concentration,$volume,$_260_280,$_260_230,$gender,$quality,$comment)
{	
	$db_connect = DB::getInstance("NGSD");
	$hash = $db_connect->prepare("INSERT INTO sample (name,name_external,sender_id,received,receiver_id,sample_type,tumor,ffpe,concentration,volume,species_id,od_260_280,od_260_230,gender,quality,comment) VALUES (:sam_name,:sam_external_name,:sam_sender,:sam_received,:sam_receiver_id,:sam_sample_type,:sam_tumor,:sam_ffpe,:sam_concentration,:sam_volume,:sam_species_id,:sam_260_280,:sam_260_230,:sam_gender,:sam_quality,:sam_comment);");
	$db_connect->bind($hash,"sam_name",$name);
	$db_connect->bind($hash,"sam_external_name",$external_name);
	
	$sender_id=get_sender_id($sender);
	if ($sender_id)
	{
		$db_connect->bind($hash,"sam_sender",$sender_id);
	}
	else
	{
		print "ERROR: Sender name \"".$sender."\" is not in db! NOT ADDED\n";
		return false;
	}
	$db_connect->bind($hash,"sam_received",replace_empty_with_null($received));
	
	$receiver_id=get_receiver_id($received_by);
	if ($receiver_id)
	{
		$db_connect->bind($hash,"sam_receiver_id",$receiver_id);
	}
	else
	{
		print "ERROR: Sender name \"".$received_by."\" is not in db! NOT ADDED\n";
		return false;
	}
	
	$db_connect->bind($hash,"sam_sample_type",$sample_type);
	$db_connect->bind($hash,"sam_tumor",convert_yes_and_no_to_1_and_0($tumor));
	$db_connect->bind($hash,"sam_ffpe",convert_yes_and_no_to_1_and_0($ffpe));

	$species_id=get_species_id($species);
	if ($species_id)
	{
		$db_connect->bind($hash,"sam_species_id",$species_id);
	}
	else
	{
		print "ERROR: Sender name \"".$species."\" is not in db! NOT ADDED\n";
		return false;
	}
	$db_connect->bind($hash,"sam_volume",replace_empty_with_null($volume));
	$db_connect->bind($hash,"sam_concentration",floatval($concentration));
	
	
	$db_connect->bind($hash,"sam_260_280",replace_empty_with_null($_260_280));
	$db_connect->bind($hash,"sam_260_230",replace_empty_with_null($_260_230));
	$db_connect->bind($hash,"sam_gender",replace_empty_with_na($gender));
	$db_connect->bind($hash,"sam_quality",replace_empty_with_na($quality));
	$db_connect->bind($hash,"sam_comment",$comment);
	$db_connect->execute($hash, true, false);
	return true;
}

function check_sample_name_exists($name)
{
	$db_connect = DB::getInstance("NGSD");
	$result_array = $db_connect->executeQuery("SELECT * FROM sample WHERE name ='".$name."'");
	if (count($result_array)>0)
	{
		return true;
	}
	return false;//not found
}

function check_external_exists($external)
{
	$db_connect = DB::getInstance("NGSD");
	$result_array = $db_connect->executeQuery("SELECT * FROM sample WHERE name_external= '".$external."'");
	if (count($result_array)>0)
	{
		return $result_array[0];
	}
	return false;//not found
}

$sample_lines = file($in);
foreach($sample_lines as $sample_line)
{
	list($sample_name,$external_name,$sender,$received,$received_by,$sample_type,$tumor,$ffpe,$species,$concentration,$volume,$_260_280,$_260_230,$gender,$quality,$comment)=explode("\t",trim($sample_line," \n\r\0\x0B"));//trim, but keep trailing tabs
	if (check_sample_name_exists($sample_name))
	{
		print "ERROR: Sample Name \"".$sample_name."\"already exists! NOT ADDED\n";	
		continue;
	}
	if ($external_name!="")
	{
		$external_exists=check_external_exists($external_name);
	}
	if (add_sample($sample_name,$external_name,$sender,$received,$received_by,$sample_type,$tumor,$ffpe,$species,$concentration,$volume,$_260_280,$_260_230,$gender,$quality,$comment))
	{
		print "ADDED $sample_line";
	}
}
?>