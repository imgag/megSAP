<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");
error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

function gender2human($gender_code)
{	
	$gender_code = trim($gender_code);
	if($gender_code=="1")
	{
		return "female";
	}
	if($gender_code=="2")
	{
		return "male";
	}
	
	return "n/a ($gender_code)";
}

//parse command line arguments
$parser = new ToolBase("genlab_gender", "Returns gender of sample ID from Genlab.");
$parser->addOutFile("id","Sample ID (e.g DX152541/DX152541_01).",false);
//optional
extract($parser->parse($argv));

//check DB is enabled


if(!GenLabDB::isEnabled())
{
	trigger_error("GenLab database is not enabled - please add credentials to settings file!", E_USER_ERROR);
}

//get gender
$db = GenLabDB::getInstance();
$res = $db->executeQuery("SELECT geschlecht,identnr FROM v_ngs_geschlecht WHERE labornummer = '{$id}'");
if(empty($res))
{
	print "No gender in GenLab for '$id'.\n";
}
foreach($res as $row)
{
	print "gender: ".gender2human($row["geschlecht"])." (identnr: ".$row["identnr"].")\n";
}
?>
