<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("batch_reanalysis", "Batch reanalysis of single processed samples.");
$parser->addStringArray("samples", "Processed sample names.", false);
$parser->addString("steps",  "Analysis steps to perform.", false);
extract($parser->parse($argv));

$count = count($samples);
for ($i=1; $i<=$count; ++$i)
{
	$ps = trim($samples[$i-1]);
	
	//determine folder name
	$info = get_processed_sample_info($ps, false);
	if (is_null($info))
	{
		print "$i/$count: Skipped '$ps' because it was not found in NGSD!\n";
		continue;
	}
	$folder = $info['ps_folder'];
	
	//check folder exists
	if(!file_exists($folder))
	{
		print "$i/$count: Skipped '$ps' because folder does not exist: $folder\n";
		continue;
	}
		
	//variant calling with multi-sample pipeline
	print "$i/$count: Analyzing '$ps'...\n";
	list($stdout, $stderr, $return) = $parser->execTool("Pipelines/analyze.php", "-name $ps -folder $folder -steps $steps", false);
	if ($return!=0)
	{
		print "  Error occurred:\n";
		print "  ".implode("  ", array_map("trim", array_merge($stdout, $stderr)));
	}
}

?>