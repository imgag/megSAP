<?php
/** 
	@page db_import_proSamples
	
	@todo implement normal_sample support!
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("db_import_proSamples", "Batch-imports processed samples to the NGSD.");
$parser->addInfile("in",  "Input sample list in TSV format.", false);
//optional
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
extract($parser->parse($argv));

//checks if the sequence given for an MID is correct
function check_mid_seq(&$db, $name, $sequence)
{
	if ($sequence=="") return;
	
	$db_seq = trim($db->getValue("SELECT sequence FROM mid WHERE name='{$name}'"));
	if ($sequence!=$db_seq)
	{
		trigger_error("MID '$name' has sequence '$db_seq' in NGSD, but '$sequence' was given. Aborting!'", E_USER_ERROR);
	}
}

//init
$db = DB::getInstance($db);
$hash = $db->prepare("INSERT INTO processed_sample (sample_id,process_id,sequencing_run_id,lane,mid1_i7,mid2_i5,operator_id,processing_system_id,project_id,molarity,comment) VALUES (:ps_sample_id,:ps_process_id,:ps_sequencing_run_id,:ps_lane,:ps_mid1_i7,:ps_mid2_i5,:ps_operator_id,:ps_processing_system_id,:ps_project_id,:ps_molarity,:ps_comment);");

//process input
$file = file($in);
foreach($file as $line)
{
	//skip empty lines and header line
	if(trim($line)=="" || starts_with($line, "sample\t")) continue;
	
	//split and trim excess whitespaces
	$parts = explode("\t",$line);
	$parts = array_map("trim", $parts);
	list($sample_name,$project,$run,$lane,$mid1_i7_name,$mid1_i7_seq,$mid2_i5_name,$mid2_i5_seq,$normal_sample,$operator,$pro_sys,$molarity,$status,$comment) = $parts;
	
	//check if processed sample already exists for this run/lane
	$sample_id = $db->getId("sample", "name", $sample_name);
	$run_id = $db->getId("sequencing_run", "name", $run);
	$ps_names = $db->getValues("SELECT CONCAT(s.name,'_',LPAD(ps.process_id,2,'0')) as ps_id FROM processed_sample ps, sequencing_run r, sample s WHERE ps.sequencing_run_id=r.id AND ps.sample_id=s.id AND s.id='$sample_id' AND r.id='$run_id' AND ps.lane='$lane'");
	if (count($ps_names)>0)
	{
		print "Notice: Sample {$sample_name} exists on lane $lane of run $run as ".implode(", ", $ps_names)." => skipped!\n";	
		continue;
	}
		
	//warn about ignored data
	if ($status!="") print "Warning: Status '$status' given, but it is ignored!\n";	
	if ($normal_sample!="") print "Warning: Normal sample '$normal_sample' given, but it is ignored!\n";	
		
	//curate input data
	$molarity = ($molarity=="" ? null : floatval($molarity));
	
	$mid1_i7_id = null;
	if($mid1_i7_name!="")
	{
		$mid1_i7_id = $db->getId("mid", "name", $mid1_i7_name);
		check_mid_seq($db, $mid1_i7_name, $mid1_i7_seq);
	}
	
	$mid2_i5_id = null;
	if($mid2_i5_name!="")
	{
		$mid2_i5_id = $db->getId("mid", "name", $mid2_i5_name);
		check_mid_seq($db, $mid2_i5_name, $mid2_i5_seq);
	}
	
	$next_process_id = $db->getValue("SELECT MAX(process_id)+1 FROM `processed_sample` WHERE sample_id=$sample_id", 1);

	//import
	$db->bind($hash, "ps_sample_id", $sample_id);
	$db->bind($hash, "ps_process_id", $next_process_id);
	$db->bind($hash, "ps_sequencing_run_id", $run_id);
	$db->bind($hash, "ps_lane", $lane);
	$db->bind($hash, "ps_mid1_i7", $mid1_i7_id);
	$db->bind($hash, "ps_mid2_i5", $mid2_i5_id);
	$db->bind($hash, "ps_operator_id", $db->getId("user", "name", $operator));
	$db->bind($hash, "ps_processing_system_id", $db->getId("processing_system", "name_manufacturer", $pro_sys));
	$db->bind($hash, "ps_project_id", $db->getId("project", "name", $project));
	$db->bind($hash, "ps_molarity", $molarity);
	$db->bind($hash, "ps_comment", $comment);
	$db->execute($hash, true);

	print "Added $sample_name\n";
}
?>