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
$parser->addString("gender", "Gender of the specimen to be checked. If not set data is taken from NGSD.",true);
//optional
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
extract($parser->parse($argv));

if(isset($gender) && $gender != "male" && $gender != "female")
{
    trigger_error("Gender can only be 'male', 'female' and 'n/a'",E_USER_ERROR);
}

//get sample info from DB
$db_conn = DB::getInstance($db);
$info = get_processed_sample_info($db_conn, $pid, false);

if(!isset($gender))
{
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
}

$sry_gene_in_target_region = false;
//chromosomal coordinates of SRY gene on chrY (GRch37)
$start_sry = 2655030;
$end_sry = 2655644;
//check wheter sry is included in target region
$target_bed = $info["sys_target"];
list($sry_stdout,$sry_stderr) = exec2("echo -e 'chrY\\t{$start_sry}\\t{$end_sry}' | " . get_path("ngs-bits")."BedIntersect -in2 {$target_bed}",false);
if($sry_stdout[0] == "chrY\t2655030\t2655644") $sry_gene_in_target_region = true;

//execute SampleGender
if ($info['sys_type']=="WGS" || $info['sys_type']=="WES") 
{
	list($stdout, $stderr) = $parser->exec(get_path("ngs-bits")."SampleGender", "-method sry -sry_cov 10 -in $in", true);
}
elseif($sry_gene_in_target_region)
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