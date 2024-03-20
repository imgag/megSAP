<?php
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("dna2ps", "Add a column with processed samples to a TSV file with DNA numbers.");
$parser->addInfile("in", "Input TSV file with DNA number if first column.", false);
$parser->addString("system", "Processing system short name.", true);
$parser->addString("project_type", "Project type.", true);
extract($parser->parse($argv));

$db = DB::getInstance("NGSD");

foreach(file($in) as $line)
{
	$line = nl_trim($line);
	if ($line=="") continue;
	
	if ($line[0]=="#")
	{
		if (!starts_with($line, "##"))
		{
			print $line."\tprocessed_samples\n";
		}
	}
	
	list($dna_nr) = explode("\t", $line, 2);
	$dna_nr = strtr($dna_nr, ["DNA-"=>"", "A1"=>"", "B1"=>"", "A2"=>"", "B2"=>"", "A3"=>"", "B3"=>""]);
	$dna_nr = trim($dna_nr);
	
	$tables = [];
	$tables[] = "sample s";
	$tables[] = "processed_sample ps";
	
	$conditions = [];
	$conditions[] = "ps.sample_id=s.id";
	$conditions[] = "ps.quality!='bad'";
	$conditions[] = "ps.id NOT IN (SELECT processed_sample_id FROM merged_processed_samples)";
	$conditions[] = "(s.name_external LIKE 'DNA{$dna_nr}%' OR s.name_external LIKE 'DX{$dna_nr}%' OR s.name_external LIKE 'DNA-{$dna_nr}%' OR s.name_external LIKE 'DX-{$dna_nr}%')";
	if ($system!="")
	{
		$sys_id = $db->getValue("SELECT id FROM `processing_system` WHERE name_short='{$system}'");
		$conditions[] = "ps.processing_system_id='$sys_id'";
	}
	if ($project_type!="")
	{
		$tables[] = "project p";
		$conditions[] = "ps.project_id=p.id";
		$conditions[] = "p.type='{$project_type}'";
	}
	
	$ps_ids = $db->getValues("SELECT CONCAT(s.name, '_0', ps.process_id) FROM ".implode(", ", $tables)." WHERE ".implode(" AND ", $conditions)." ORDER BY ps.id ASC");
	print $line."\t".implode(", ", $ps_ids)."\n";

}

?>