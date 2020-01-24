<?php
/** 
	@page db_check_gender
*/
error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

$parser = new ToolBase("db_check_gender", "Checks the gender of a sample.");
$parser->addInfile("in",  "Input file in BAM format.", false);
$parser->addString("pid",  "Processed sample ID, e.g. GS120001_01, used to determine the sample gender from the NGSD.", false);
$parser->addString("gender", "Gender of the input. If unset, gender is looked up in NGSD.", true);
//optional
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
extract($parser->parse($argv));

//determine method and method parameters
$method = "hetx"; //this is the slowest but safest method, always use this if nothing better is available
$args = "";
if(isset($gender)) //gender given > no NGSD connection is used, so nothing can be checked
{
	if ($gender!="male" && $gender!="female") trigger_error("Gender can only be 'male', 'female' and 'n/a'", E_USER_ERROR);
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
	
	if ($info['sys_type']=="WGS" || $info['sys_type']=="WGS (shallow)")
	{
		$method = "xy";
	}
	else if ($info['sys_type']=="WES")
	{
		$method = "xy"; //defaults work for WES as well, e.g. in ssHAEv7 the distribution is 15% chrY in males and 0.2% chrY in females
	}
	else if ($info['sys_type']=="RNA")
	{
		$method = "xy";
		$args = "-min_male 0.012 -max_female 0.008";
	}
	else //check if sry is included in target region
	{
		list($stdout, $stderr) = exec2("echo -e 'chrY\\t2655030\\t2655644' | ".get_path("ngs-bits")."BedIntersect -in2 ".$info["sys_target"], false); //works for GRCh37 only
		if ($stdout[0] == "chrY\t2655030\t2655644") //SRY is in roi
		{
			$method = "sry";
			$args = "-sry_cov 10";
		}
	}
}

//determine gender from BAM
list($stdout, $stderr) = $parser->exec(get_path("ngs-bits")."SampleGender", "-in {$in} -method {$method} {$args}", true);
$gender2 = explode("\t", $stdout[1])[1];
if (starts_with($gender2, "unknown"))
{
	$parser->log("Could not check gender, because we could not determine it from BAM file.");
	exit(0);
}

//make sure they match
if ($gender!=$gender2)
{
	trigger_error("Expected gender '{$gender}', but determined gender '{$gender2}' using method {$method} on BAM file {$in}", E_USER_ERROR);
}

?>