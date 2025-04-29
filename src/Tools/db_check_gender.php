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
$parser->addInt("sry_cov", "Minimum SRY coverage to consider a sample as male.", true, "10");
$parser->addFlag("check_sry_cov", "Check coverage of SRY additionally to the prefered method (hetx/xy)");
$parser->addString("build", "Genome build of sample.", true, "GRCh38");
//optional
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
extract($parser->parse($argv));


function sry_in_target($parser, $sys_roi)
{
	if ($sys_roi != "")
	{
		$ngsbits_command = $parser->execApptainer("ngs-bits", "BedIntersect", "-in2 ".$sys_roi, [$sys_roi], [], true, true, false);
		list($stdout, $stderr) = exec2("echo -e 'chrY\\t2786989\\t2787603' | $ngsbits_command", false); //works for GRCh38 only

		return $stdout[0] == "chrY\t2786989\t2787603"; //SRY is in roi
	}
	return false;
}

function check_gender($gender, $gender2, $method, $in)
{
	if (starts_with($gender2, "unknown"))
	{
		global $parser;
		$parser->log("Could not check gender, because we could not determine it from BAM file with method {$method}.");
		return;
	}

	//make sure they match
	if ($gender!=$gender2)
	{
		trigger_error("Expected gender '{$gender}', but determined gender '{$gender2}' using method {$method} on BAM file {$in}", E_USER_ERROR);
	}

}

 //check parameters
if(isset($gender) && ($gender!="male" && $gender!="female"))
{
	trigger_error("Gender can only be 'male' or 'female'!", E_USER_ERROR);
}

//determine method and parameters (skip if gender is set)
$method = "hetx";
$args = "";
$sys_roi = "";

if (!isset($gender))
{
	if (! db_is_enabled("NGSD"))
	{
		trigger_error("Could not determine gender for processed sample '{$pid}': NGSD is not enabled.", E_USER_WARNING);
		exit(0);
	}
	
	//get sample info from DB
	$db = DB::getInstance($db);
	$ps_id = get_processed_sample_id($db, $pid, false);
	if ($ps_id==-1)
	{
		$parser->log("Could not determine gender for processed sample '{$pid}': processed sample is not found in NGSD.");
		exit(0);
	}
	$info = get_processed_sample_info($db, $pid);
	$gender = $info['gender'];
	$sys_type = $info['sys_type'];
	$sys_roi = trim(realpath($info["sys_target"]));

	//check gender from NGSD is set
	if ($gender=='n/a')
	{
		$parser->log("Could not determine gender for processed sample '{$pid}': gender not specified in NGSD.");
		exit(0);
	}

	//skip for cf-DNA
	if ($sys_type=="cfDNA (patient-specific)")
	{
		$parser->log("Could not determine gender for processed sample '{$pid}' : gender not specified in NGSD.");
		exit(0);
	}

	//determine method and arguments
	if ($sys_type=="WGS" || $sys_type=="WES")
	{
		$method = "hetx";
	}
	else if ($sys_type=="WGS (shallow)")
	{
		$method = "xy";
	}
	else if ($sys_type=="RNA")
	{
		$method = "xy";
		$args = "-min_male 0.012 -max_female 0.008";
	}
	else if ($sys_type=="lrGS")
	{
		$method = "hetx";
		$args = "-include_single_end_reads";
	}
	else
	{
		if (sry_in_target($parser, $sys_roi)) //SRY is in roi
		{
			$method = "sry";
			$args = "-sry_cov $sry_cov";
		}
	}
}

//determine gender from BAM
$genome = genome_fasta($build);
list($stdout, $stderr) = $parser->execApptainer("ngs-bits", "SampleGender", "-in {$in} -method {$method} {$args} -build ".ngsbits_build($build)." -ref {$genome}", [$in, $genome]);
$gender2 = explode("\t", $stdout[1])[1];

check_gender($gender, $gender2, $method, $in);

if ($check_sry_cov && $method != "sry")
{
	if (sry_in_target($parser, $sys_roi)) //SRY is in roi
	{
		list($stdout, $stderr) = $parser->execApptainer("ngs-bits", "SampleGender", "-in {$in} -method sry -sry_cov {$sry_cov} -build ".ngsbits_build($build)." -ref {$genome}", [$in, $genome]);
		$gender_sry = explode("\t", $stdout[1])[1];
	
		check_gender($gender, $gender_sry, $method, $in);
	}
	
	trigger_error("Flag 'check_sry_cov' given but target region of the processed sample {$pid} doesn't include SRY.", E_USER_WARNING);
}

?>