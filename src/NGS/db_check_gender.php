<?php
/** 
	@page db_check_gender
*/
error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

$parser = new ToolBase("db_check_gender", "Checks that the gender of a sample matches the DB meta information.");
$parser->addInfile("in",  "Input file in BAM format.", false);
$parser->addString("pid",  "Processed sample ID, e.g. GS120001_01, used to determine the processing system.", false);
//optional
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
extract($parser->parse($argv));

//establish database connection
$db_connect = DB::getInstance($db);

//get gender of sample from DB
$hash = $db_connect->prepare("SELECT gender FROM sample WHERE name=:sname");
$db_connect->bind($hash, "sname", substr($pid,0,-3));
$db_connect->execute($hash, TRUE);
$result = $db_connect->fetch($hash);
if (count($result)==0)
{
	$parser->log("Could not determine gender for PID $pid from DB: no sample entry.");
	exit(0);
}
$gender = $result[0]['gender'];
if ($gender=='n/a')
{
	$parser->log("Could not determine gender for PID $pid from DB: gender not set in sample entry.");
	exit(0);
}

//determine gender from BAM file
list($stdout, $stderr) = $parser->exec(get_path("ngs-bits")."SampleGender", "-method hetx -in $in", true);
$gender2 = trim(substr($stdout[count($stdout)-1], 7));
if (starts_with($gender2,"unknown"))
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