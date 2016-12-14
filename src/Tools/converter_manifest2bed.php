<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("converter_manifest2bed", "Convert manifest to bed.");
$parser->addInfile("m",  "Illumina manifest file.", false);
$parser->addOutfile("o",  "Output bed file.", false);
//optional
$parser->addInt("t",  "Target number (column 3 within manifest file).", true, 1);
extract($parser->parse($argv));

$manifest_file = file($m, FILE_IGNORE_NEW_LINES | FILE_SKIP_EMPTY_LINES);

//move to 
$targets = false;
$output = array();
for($i=0;$i<count($manifest_file);++$i)
{
	
	//find targets
	if(strpos($manifest_file[$i], "[Targets]")===0) 
	{
		$targets = true;
		++$i; //-> skip headers
		continue;
	}
	
	$cols = explode("\t", $manifest_file[$i]);
	if($targets && $cols[2]==$t)
	{
		$output[] = $cols[3]."\t".($cols[4]-1)."\t".$cols[5]."\t".$cols[0]; //-> for targets
	}
}

file_put_contents($o, implode("\n",$output));
