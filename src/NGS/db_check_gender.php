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

if(isset($gender)) //gender given > no NGSD connection is used
{
	if ($gender!="male" && $gender!="female") trigger_error("Gender can only be 'male', 'female' and 'n/a'", E_USER_ERROR);
	
	//we don't know and cannot check
	$is_wgs = false;
	$sry_in_roi = false;
	$is_rna = false;
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
		$is_wgs = true;
		$sry_in_roi = true;
		$is_rna = false;
	}
	else if ($info['sys_type']=="WES")
	{
		$is_wgs = false;
		$sry_in_roi = true;
		$is_rna = false;
	}
	else if ($info['sys_type']=="RNA")
	{
		$is_wgs = false;
		$sry_in_roi = false;
		$is_rna = true;
	}
	else //check if sry is included in target region
	{
		$is_wgs = false;
		list($stdout, $stderr) = exec2("echo -e 'chrY\\t2655030\\t2655644' | ".get_path("ngs-bits")."BedIntersect -in2 ".$info["sys_target"], false); //works for GRCh37 only
		$sry_in_roi = ($stdout[0] == "chrY\t2655030\t2655644");
		$is_rna = false;
	}
}

//determine gender from BAM
$extra = "";
if ($is_wgs)
{
	$method = "xy";
}
else if ($is_rna)
{
	$method = "xy";
	$extra = "-min_male 0.012 -max_female 0.008";
}
else if ($sry_in_roi)
{
	$method = "sry";
	$extra = "-sry_cov 10";
}
else
{
	$method = "hetx";
}
list($stdout, $stderr) = $parser->exec(get_path("ngs-bits")."SampleGender", "-in {$in} -method {$method} {$extra}", true);
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