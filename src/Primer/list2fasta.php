<?php

/** 
	@page list2fasta
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("list2fasta", "Converts text input format to FASTA.");
$parser->addInfile("in",  "Input primer TXT file (name_for, name_rev, seq_for, seq_ref, chr, prod_length).", false);
$parser->addOutfile("out",  "Output FASTA file.", false);
extract($parser->parse($argv));

$names = array();
$input = file($in);
$output = array();
foreach($input as $line)
{
	if (trim($line)=="" || $line[0]=="#" || starts_with($line, "Name FOR")) continue;
	
	//skip invalid lines
	$parts = explode("\t", trim($line));
	if (count($parts)<6)
	{
		trigger_error("Skipping invalid line: $line", E_USER_WARNING);
		continue;
	}
	
	list($n1, $n2, $s1, $s2, $chr, $len) = $parts;
	$n1 = trim($n1);
	$n2 = trim($n2);
	$s1 = trim($s1);
	$s2 = trim($s2);
	
	//warn when strange sequences are inserted
	if(!preg_match("/^[ACTGacgt]+$/", $s1))
	{
		trigger_error("Primer $n1 contains invalid bases: $s1", E_USER_WARNING);
	}
	if(!preg_match("/^[ACTGacgt]+$/", $s2))
	{
		trigger_error("Primer $n2 contains invalid bases: $s2", E_USER_WARNING);
	}
	
	//warn when names are used several times
	if (in_array($n1, $names))
	{
		trigger_error("Sequence name $n1 used several times!", E_USER_WARNING);
	}
	else
	{
		$names[] = $n1;
	}
	if (in_array($n2, $names))
	{
		trigger_error("Sequence name $n2 used several times!", E_USER_WARNING);
	}
	else
	{
		$names[] = $n2;
	}
	
	//output
	$output[] = ">$n1";
	$output[] = "$s1";
	$output[] = ">$n2";
	$output[] = "$s2";
}
file_put_contents($out, implode("\n", $output));

?>