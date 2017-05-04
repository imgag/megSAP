<?php
/** 
	@page db_import_samples
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("db_import_samples", "Batch-imports samples to the NGSD.");
$parser->addInfile("in",  "Input sample list in TSV format.", false);
//optional
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
extract($parser->parse($argv));

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

//init
$db = DB::getInstance($db);
$hash = $db->prepare("INSERT INTO sample (name,name_external,sender_id,received,receiver_id,sample_type,tumor,ffpe,concentration,volume,species_id,od_260_280,od_260_230,gender,quality,comment) VALUES (:sam_name,:sam_external_name,:sam_sender,:sam_received,:sam_receiver_id,:sam_sample_type,:sam_tumor,:sam_ffpe,:sam_concentration,:sam_volume,:sam_species_id,:sam_260_280,:sam_260_230,:sam_gender,:sam_quality,:sam_comment);");

//process input
$file = file($in);
foreach($file as $line)
{
	//skip empty lines and header line
	if(trim($line)=="" || starts_with($line, "name\t")) continue;
	
	//split and trim excess whitespaces
	$parts = explode("\t",$line);
	$parts = array_map("trim", $parts);
	list($sample_name,$external_name,$sender,$received,$received_by,$sample_type,$tumor,$ffpe,$species,$concentration,$volume,$od_260_280,$od_260_230,$gender,$quality,$comment) = $parts;
	
	//check if sample already exists in DB
	$sample_id = $db->getId("sample", "name", $sample_name, false);
	if ($sample_id!=-1)
	{
		print "Notice: Sample {$sample_name} exists with id $sample_id => skipped!\n";	
		continue;
	}
	
	//import
	$db->bind($hash,"sam_name",$sample_name);
	$db->bind($hash,"sam_external_name",$external_name);
	$db->bind($hash,"sam_sender", $db->getId("sender", "name", $sender));
	$db->bind($hash,"sam_received",replace_empty_with_null($received));
	$db->bind($hash,"sam_receiver_id",$db->getId("user", "name", $received_by));
	$db->bind($hash,"sam_sample_type",$sample_type);
	$db->bind($hash,"sam_tumor",convert_yes_and_no_to_1_and_0($tumor));
	$db->bind($hash,"sam_ffpe",convert_yes_and_no_to_1_and_0($ffpe));
	$db->bind($hash,"sam_species_id",$db->getId("species", "name", $species));
	$db->bind($hash,"sam_concentration",floatval($concentration));
	$db->bind($hash,"sam_volume",replace_empty_with_null($volume));
	$db->bind($hash,"sam_260_280",replace_empty_with_null($od_260_280));
	$db->bind($hash,"sam_260_230",replace_empty_with_null($od_260_230));
	$db->bind($hash,"sam_gender",replace_empty_with_na($gender));
	$db->bind($hash,"sam_quality",replace_empty_with_na($quality));
	$db->bind($hash,"sam_comment",$comment);
	$db->execute($hash, true);
	
	print "Added $sample_name";
}
?>