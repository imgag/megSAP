<?php
/** 
	@page db_delete_variants
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("db_delete_variants", "Deletes variants from NGSD.");
$parser->addInfile("in", "Input VCF file.", false);
$parser->addFlag("del_class", "Also delete entries in variant_classification if present.", false);
//optional
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
extract($parser->parse($argv));

//init
$db = DB::getInstance($db);
$genome = genome_fasta("GRCh38");

//check input VCF
list($stderr, $stdout, $exit_code) = $parser->execApptainer("ngs-bits", "VcfCheck", "-in {$in} -ref $genome", [$in, $genome], [], false, true, false, false);
if ($exit_code!=0)
{
	trigger_error("Invalid input VCF:\n".implode("\n", $stderr), E_USER_ERROR);
}

//import data
$file = file($in);
foreach($file as $line)
{
	$line = nl_trim($line);
	if($line=="" || $line[0]=="#") continue;
	
	//split and check column count
	$parts = explode("\t", $line);
	if (count($parts)<7) trigger_error("Input VCF line with less than 7 columns found: ".$line, E_USER_ERROR);
	list($chr, $pos, , $ref, $obs, , , $info) = $parts;
	$tag = "$chr:$pos $ref>$obs";
	
	//convert VCF to GSvar format
	$start = trim($pos);
	$end = $start;
	$ref = strtoupper(trim($ref));
	$obs = strtoupper(trim($obs));
	if(strlen($ref)>1 || strlen($obs)>1) //correct indels
	{
		list($start, $end, $ref, $obs) = correct_indel($start, $ref, $obs);
	}
	
	//get variant ID	
	$v_id = get_variant_id($db, $chr, $start, $end, $ref, $obs);
	if ($v_id==-1)
	{
		print "$tag skipped: not in NGSD!\n";
		continue;
	}
	
	//check if variant is used in other tables
	$ref_tables = ["variant_publication","somatic_vicc_interpretation","detected_variant","detected_somatic_variant","somatic_report_configuration_variant","somatic_report_configuration_germl_var","report_configuration_variant","variant_validation","variant_literature"];
	if (!$del_class) $ref_tables[] = "variant_classification";
	$hits = [];
	foreach($ref_tables as $table)
	{
		$count = $db->getValue("SELECT COUNT(*) FROM $table WHERE variant_id=$v_id");
		if ($count==0) continue;
		
		if (!in_array($table, $hits)) $hits[] = $table;
	}
	if (count($hits)>0)
	{
		print "$tag skipped: used in the tables ".implode(", ", $hits)."\n";
		continue;
	}

	print "$tag ($v_id) deleting...\n";
	if ($del_class)
	{
		//print "DELETE FROM variant_classification WHERE variant_id=$v_id\n";
		$db->executeStmt("DELETE FROM variant_classification WHERE variant_id=$v_id");
	}
	
	//print "DELETE FROM variant WHERE id=$v_id\n";
	$db->executeStmt("DELETE FROM variant WHERE id=$v_id");
}

?>