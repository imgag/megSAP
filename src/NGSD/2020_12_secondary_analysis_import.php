<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("2020_12_secondary_analysis_import", "Imports secondary anayses.");
$parser->addInfile("folders", "File containing the paths of secondary analyses to import.", false);
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
extract($parser->parse($argv));

//init
$folders = file($folders);
$db = DB::getInstance($db);

foreach($folders as $folder)
{
	$folder = trim($folder);
	if ($folder=="") continue;
	
	//check that folder exists
	if (!file_exists($folder))
	{
		print "$folder: missing > skipped!\n";
		continue;
	}
	
	//check that folder contains
	$basename = basename($folder);
	if (starts_with($basename, "Trio_"))
	{
		$type = "trio";
		$gsvar = "{$folder}/trio.GSvar";
	}
	else if (starts_with($basename, "Multi_"))
	{
		$type = "multi sample";
		$gsvar = "{$folder}/multi.GSvar";
	}
	else if (starts_with($basename, "Somatic_"))
	{
		$type = "somatic";
		$gsvar = "{$folder}/".substr($basename, 8).".GSvar";
	}
	else
	{
		trigger_error("Unhandled folder type: $folder", E_USER_ERROR);
	}
	
	//check GSvar exists
	if (!file_exists($gsvar))
	{
		print "$gsvar: missing > skipped!\n";
		
		$files = glob("{$folder}/*.*");
		if(count($files)==1 && ends_with($files[0], ".log"))
		{
			print "  only log file > deleting it!\n";
			exec2("rm ".$files[0]);
			exec2("rmdir ".$folder);
		}
		else
		{
			//print_r($files);
		}
		continue;
	}
	
	//check processed sample names are in NGSD
	$ps_list = [];
	$tmp = explode("_", $basename, 2)[1]."-";
	$start = 0;
	for ($i=0; $i<strlen($tmp); ++$i)
	{
		if ($tmp[$i]=="-")
		{
			$ps_list[] = substr($tmp, $start, $i-$start);
			$start = $i+1;
		}
		if ($tmp[$i]=="_" && $tmp[$i-3]=="_")
		{
			$ps_list[] = substr($tmp, $start, $i-$start);
			$start = $i+1;
		}
	}
	$not_found = [];
	foreach($ps_list as $ps)
	{
		list($sample, $process_id) = explode("_", $ps."_");
		$ps_id = $db->getValue("SELECT ps.id FROM processed_sample as ps, sample as s WHERE ps.sample_id = s.id AND s.name='{$sample}' AND ps.process_id='{$process_id}'", -1);
		if ($ps_id==-1)
		{
			$not_found[] = $ps;
		}
	}
	if (count($not_found)>0)
	{
		print "$gsvar: processed samples not found: ".implode(", ", $not_found)." > skipped!\n";
		continue;
	}
	
	//add to NGSD
	$db->executeStmt("INSERT INTO `secondary_analysis`(`type`, `gsvar_file`) VALUES ('{$type}','{$gsvar}')");
}

?>