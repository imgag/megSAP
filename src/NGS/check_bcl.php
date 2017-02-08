<?php

/**
	@page check_bcl	
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("check_bcl", "Checks which tiles of an Illumina GAIIx lane are missing or corrupt.");
$parser->addInfile("in",  "Run folder.", false);
$parser->addFloat("lane",  "Lane number.", false);
extract($parser->parse($argv));

// find missing tiles in all cycles
$folder = $in."/Data/Intensities/BaseCalls/L00".$lane."/";

$count_cycles = 0;
$missing = array();
$empty = array();
$files = scandir($folder);
foreach($files as $cycle)
{
	$dir = $folder.$cycle."/";
	
	if (is_dir($dir) && $cycle[0]=="C")
	{
		++$count_cycles;
				
		$files2 = scandir($dir);
		
		for ($i=1; $i<121; ++$i)
		{
			$bcl_file = "s_".$lane."_".$i.".bcl";
			$sta_file = "s_".$lane."_".$i.".stats";
			if (!in_array($bcl_file, $files2) || !in_array($sta_file, $files2))
			{
				$missing[] = "s_".$lane."_".$i;
			}
			else if (filesize($dir.$bcl_file)<20 || filesize($dir.$sta_file)<20)
			{
				$empty[] = "s_".$lane."_".$i;
			}
		}
	}
}

//make unique and sort
sort($missing);
$missing = array_unique($missing);
sort($empty);
$empty = array_unique($empty);

//output
print "Folder: $in\n";
print "Lane: $lane\n";
print "Cycles: $count_cycles\n";
print "Missing: ".implode(", ", $missing)."\n";
print "Empty: ".implode(", ", $empty)."\n";

//print cassava option
$present = array();
for ($i=1; $i<121; ++$i)
{
	if (!in_array("s_".$lane."_".$i, $missing) && !in_array("s_".$lane."_".$i, $empty))
	{
		$present[] =  "s_".$lane."_".str_pad($i, 4, "0", STR_PAD_LEFT);
	}
}
print "\n";
print "--tiles ".implode(",", $present);
print "\n";


?>
