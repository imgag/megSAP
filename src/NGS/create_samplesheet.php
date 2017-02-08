<?php

/**
	@page create_samplesheet
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("create_samplesheet", "Creates a CASAVA-compatible sample sheet from a TSV file.");
$parser->addInfile("in",  "Input sample file with four columns: lane, sample name, p7 index, p5 index.", false);
$parser->addOutfile("out", "Output sample sheet file name.", false);
$parser->addString("run",  "Run number, e.g. '304'.", false);
$parser->addString("fcid",  "Flow cell identifier, e.g. 'C3NPWACXX'.", false);
extract($parser->parse($argv));

$run = str_pad($run, 5, "0", STR_PAD_LEFT);

//determine index length
$p7l = 100;
$p5l = 100;
$data = array();
$file = file($in);
foreach($file as $line)
{
	if ($line=="" || $line[0]=="#") continue;
	$parts = explode("\t", $line);
	if (count($parts)!=4)
	{
		trigger_error("Invalid input line: $line has ".count($parts)." columns", E_USER_ERROR);
	}
	list($lanes, $sample, $p7, $p5) = $parts;
	$p5=trim($p5);
	$sample=str_replace(array("Ä","Ü","Ö","ä","ü","ö"," ","ß"),array("AE","UE","OE","ae","ue","oe","_","ss"),$sample); //replace german special characters and spaces
	$sample=preg_replace("/[^a-zA-Z\d-_]/", '_', $sample); //replace all remaining chars that are neither alphanumeric nor underscore nor minus
	$data[$lanes][$sample] = array($p7, $p5);
	$p7l = min($p7l, strlen($p7));
	$p5l = min($p5l, strlen($p5));
}
print "p7 length: $p7l\n";
print "p5 length: $p5l\n";

//create joined index sequence
foreach($data as $lanes => $tmp)
{
	foreach($tmp as $sample => $tmp2)
	{
		list($p7, $p5) = $tmp2;
		$p5_rev = strtr(strtoupper(strrev($p5)), array("A"=>"T","C"=>"G","G"=>"C","T"=>"A"));	
		$data[$lanes][$sample] = substr($p7, 0, $p7l).substr($p5_rev, 0, $p5l);
	}
}

//check that the indices don't overlap
foreach($data as $lanes => $tmp)
{
	$mids = array_values($tmp);
	for ($i=0; $i<count($mids); ++$i)
	{
		for ($j=$i+1; $j<count($mids); ++$j)
		{
			$dist = levenshtein($mids[$i], $mids[$j]);
			if ($dist==0)
			{
				trigger_error("MID clash in lanes $lanes ".$mids[$i]." <=> ".$mids[$j], E_USER_ERROR);
			}
			else if ($dist<2)
			{
				print "Warning: MID similarity in lanes $lanes ".$mids[$i]." <=> ".$mids[$j]." (distance $dist)\n";
			}
		}
	}
}

//create samplesheet
$output = array();
$output[] = "FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,Project\n";
foreach($data as $lanes => $tmp)
{
	$lanes = explode("+", $lanes);
	foreach($tmp as $sample => $mid)
	{
		foreach($lanes as $lane)
		{
			$output[] = "$fcid,".trim($lane).",$sample,hg19,$mid,-,N,-,-,$run\n";
		}
	}
}
file_put_contents($out, $output);

?>
