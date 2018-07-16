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
$parser->addString("gender", "Gender of the input. If unset, gender is looked up in NGSD.",true);
//optional
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
extract($parser->parse($argv));

if(isset($gender)) //gender given > no NGSD connection is used
{
	if ($gender!="male" && $gender!="female") trigger_error("Gender can only be 'male', 'female' and 'n/a'", E_USER_ERROR);
	
	$sry_in_roi = false; //we don't know > assume it is not
}
else
{
	//get sample info from DB
	$db_conn = DB::getInstance($db);
	$info = get_processed_sample_info($db_conn, $pid, true);
    $gender = $info['gender'];
    if ($gender=='n/a')
    {
        $parser->log("Could not determine gender for processed sample '$pid' from DB '$db': gender not set in sample entry.");
        exit(0);
    }
	
	if ($info['sys_type']=="WGS" || $info['sys_type']=="WES")
	{
		$sry_in_roi = true;
	}
	else //check if sry is included in target region
	{
		list($stdout, $stderr) = exec2("echo -e 'chrY\\t2655030\\t2655644' | " . get_path("ngs-bits")."BedIntersect -in2 ".$info["sys_target"], false); //works for GRCh37 only
		$sry_in_roi = ($stdout[0] == "chrY\t2655030\t2655644");
	}
}

//execute SampleGender
if ($sry_in_roi) 
{
	list($stdout, $stderr) = $parser->exec(get_path("ngs-bits")."SampleGender", "-method sry -sry_cov 10 -in $in", true);
}
else
{
	list($stdout, $stderr) = $parser->exec(get_path("ngs-bits")."SampleGender", "-method hetx -in $in", true);
}

//determine gender
$gender2 = explode("\t", $stdout[1])[1];
if ($gender2=="")
{
	trigger_error("Could not determine gender from SampleGender output!", E_USER_ERROR);
}
if (starts_with($gender2, "unknown"))
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