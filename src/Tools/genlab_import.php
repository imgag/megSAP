<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");
error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("genlab_import", "Imports HPO/ICD10 data from GenLab into the NGSD.");
extract($parser->parse($argv));

if(!db_is_enabled("GL8"))
{
	trigger_error("Could not connect to GenLab database.", E_USER_ERROR);
}
$db = DB::getInstance("GL8");
$db2 = DB::getInstance("NGSD");
$user_id = $db2->getValue("SELECT id FROM user WHERE user_id='genlab_import'");

//get data from GenLab
$genlab_file = "genlab_import.tsv";
if (!file_exists($genlab_file))
{
	$output = array();
	$res = $db->executeQuery("SELECT labornummer, HPOTERM1, HPOTERM2, HPOTERM3, HPOTERM4, ICD10DIAGNOSE, TUMORANTEIL FROM v_ngs_sap");
	foreach($res as $row)
	{
		$output[] = $row['labornummer']."\t".$row['HPOTERM1']."\t".$row['HPOTERM2']."\t".$row['HPOTERM3']."\t".$row['HPOTERM4']."\t".$row['ICD10DIAGNOSE']."\t".$row['TUMORANTEIL']."\n";
	}
	file_put_contents($genlab_file, $output);
}

//load GenLab data
$i = 0;
$file = file($genlab_file);
foreach($file as $line)
{
	++$i;
	$line = nl_trim($line);
	$hpo = array();
	list($lab, $hpo[], $hpo[], $hpo[], $hpo[], $icd10, $fraction) = explode("\t", $line);
	print "LAB $i: $lab - ";
	
	//determine what to import
	$import = array();
	
	$icd10 = trim($icd10);
	if ($icd10!="" && $icd10!="-") $import[] = array("ICD10 code", $icd10);
	
	foreach($hpo as $term)
	{
		$term = trim($term);
		if (starts_with($term, "HP:")) $import[] = array("HPO term id", $term);
	}
	
	$fraction = trim($fraction);
	if ($fraction!="" && $fraction!="-" && $fraction!=0.0) $import[] = array("tumor fraction", $fraction);
	
	if (count($import)==0)
	{
		print " SKIPPED: nothing to import\n";
		continue;
	}
	
	//fix processed sample instead of sample
	$sample_name = $lab;
	if (ends_with($sample_name, "_01") || ends_with($sample_name, "_02"))
	{
		$sample_name = substr($sample_name, 0, -3);
	}
	
	//find sample ID
	$sample_id = $db2->getValue("SELECT id FROM sample WHERE name='$sample_name'", -1);
	if ($sample_id=="-1")
	{
		print " SKIPPED: sample '$sample_name' not found in NGSD\n";
		continue;	
	}
	
	//determine which information is missing
	$import_missing = array();
	foreach($import as list($type, $entry))
	{
		$entry_id = $db2->getValue("SELECT id FROM sample_disease_info WHERE sample_id='$sample_id' AND disease_info='$entry'", -1);
		if ($entry_id==-1)
		{
			$import_missing[] = array($type, $entry);
		}
	}
	if (count($import_missing)==0)
	{
		print " SKIPPED: all information already present\n";
		continue;
	}
	
	// import data into NGSD
	foreach($import_missing as list($type, $entry))
	{
		$db2->executeStmt("INSERT INTO sample_disease_info (`sample_id`, `disease_info`, `type`, `user_id`, `date`) VALUES ($sample_id, '$entry', '$type', $user_id, CURRENT_TIMESTAMP())");
	}
	print " IMPORTED (".count($import_missing).")\n";
}

?>