<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");
	
function get_af($info)
{
	$parts = explode(";", $info);
	foreach($parts as $entry)
	{
		if (starts_with($entry, "AF="))
		{
			return substr($entry, 3);
		}
	}
	
	trigger_error("No AF entry found in '$entry'!", E_USER_ERROR);
}

//merge same variants to one line (sum of AF) and remove all unneeded information (ID, AC, AN, ...)
$chr = null;
$pos = null;
$afs = array();
$handle = fopen("php://stdin", "r");
while($line = fgets($handle))
{
	//header
	$line = trim($line);
	if ($line=="") continue;
	if($line[0]=="#")
	{
		print $line."\n";
		continue;
	}
	
	//content
	$parts = explode("\t", $line);
	list($c, $p, , $r, $a, , , $i) = $parts;
	
	//skip special variants (containing "N", copy-number variants, etc.)
	if (contains($r, "N") || contains($a, "N") || contains($a, "<")) continue;
	
	$f = get_af($i);
	if ($p!=$pos || $c!=$chr) //position changed 
	{
		//write out variants of last position
		if (!is_null($chr))
		{
			foreach($afs as $tag => $af)
			{
				list($ref, $alt) = explode(":", $tag);
				print "$chr	$pos	.	$ref	$alt	.	.	AF=".number_format($af, 7)."\n";
			}
		}
		
		//init position/AFs
		$chr = $c;
		$pos = $p;
		$afs = array("$r:$a" => $f);
	}
	else //same position as last variant
	{
		//init/sum AFs
		$tag = "$r:$a";
		$afs[$tag] = isset($afs[$tag]) ? $afs[$tag] + $f : $f;
	}
}
fclose($handle);

//write out variants of last position
if (!is_null($chr))
{
	foreach($afs as $tag => $af)
	{
		list($ref, $alt) = explode(":", $tag);
		print "$chr	$pos	.	$ref	$alt	.	.	AF=".number_format($af, 7)."\n";
	}
}

?>