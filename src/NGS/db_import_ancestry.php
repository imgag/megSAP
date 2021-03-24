<?php
/** 
	@page db_import_ancestry
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("db_import_ancestry", "Imports QC terms to NGSD.");
$parser->addString("id", "Processes sample identifier (e.g. GS000123_01).", false);
$parser->addInfile("in", "SampleAncestry output file.", false);
//optional
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
extract($parser->parse($argv));

//init
$db = DB::getInstance($db);
$psid = get_processed_sample_id($db, $id);

//parse file
$file = file($in);
if (count($file)<2) trigger_error("Ancestry file '{$in}' is invalid - it has less than two lines!", E_USER_ERROR);

$line = trim($file[1]);
$parts = explode("\t", $line);
if (count($parts)!=7) trigger_error("Ancestry file '{$in}' is invalid - line does not have 7 parts: {$line}", E_USER_ERROR);

//import
list($name, $snps, $afr, $eur, $sas, $eas, $pop) = explode("\t", $line);
if ($pop=="NOT_ENOUGH_SNPS")
{
	trigger_error("Ancestry file '{$in}' skipped (not enough SNPs).", E_USER_NOTICE);
}
else
{
	//delete if already exists
	$db->executeStmt("DELETE FROM `processed_sample_ancestry` WHERE processed_sample_id={$psid}");
	
	$hash = $db->prepare("INSERT INTO `processed_sample_ancestry`(`processed_sample_id`, `num_snps`, `score_afr`, `score_eur`, `score_sas`, `score_eas`, `population`) VALUES (:0, :1, :2, :3, :4, :5, :6)");
	$db->bind($hash, "0", $psid);
	$db->bind($hash, "1", $snps);
	$db->bind($hash, "2", $afr);
	$db->bind($hash, "3", $eur);
	$db->bind($hash, "4", $sas);
	$db->bind($hash, "5", $eas);
	$db->bind($hash, "6", $pop);
	$db->execute($hash, true);
}


?>
