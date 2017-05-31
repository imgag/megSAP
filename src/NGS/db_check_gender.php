<?php
/** 
	@page db_check_gender
	
	@todo Also use SRY-based method for panels that contain the SRY gene
*/
error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

$parser = new ToolBase("db_check_gender", "Checks that the gender of a sample matches the DB meta information.");
$parser->addInfile("in",  "Input file in BAM format.", false);
$parser->addString("pid",  "Processed sample ID, e.g. GS120001_01, used to determine the sample gender from the DB.", false);
//optional
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
extract($parser->parse($argv));

//get sample info from DB
$info = get_processed_sample_info($pid, false, $db);
if (is_null($info))
{
	$parser->log("Could not determine gender for processed sample '$pid' from DB: no sample entry.");
	exit(0);
}
$gender = $info['gender'];
if ($gender=='n/a')
{
	$parser->log("Could not determine gender for processed sample '$pid' from DB: gender not set in sample entry.");
	exit(0);
}

//execute SampleGender
if ($info['sys_type']=="WGS" || $info['sys_type']=="WES") 
{
	list($stdout, $stderr) = $parser->exec(get_path("ngs-bits")."SampleGender", "-method sry -in $in", true);
}
else
{
	list($stdout, $stderr) = $parser->exec(get_path("ngs-bits")."SampleGender", "-method hetx -in $in", true);
}

//determine gender
$gender2 = "";
foreach($stdout as $line)
{
	if (starts_with($line, "gender: "))
	{
		$gender2 = trim(substr($line, 7));
	}
}
if ($gender2=="")
{
	trigger_error("Could not determine gender from SampleGender output!", E_USER_ERROR);
}
else if (starts_with($gender2, "unknown"))
{
	$parser->log("Could not check gender, because we could not determine it from BAM file.");
	exit(0);
}

//make sure they match
if ($gender!=$gender2)
{
	trigger_error("Gender from DB '$gender' does not match gender from BAM file '$gender2'!", E_USER_ERROR);
}

?>