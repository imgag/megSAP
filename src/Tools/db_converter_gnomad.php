<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

$handle = fopen("php://stdin", "r");
while($line = fgets($handle))
{
	$line = trim($line);
	if ($line=="") continue;
	if($line[0]=="#")
	{
		if (in_array("-header", $argv))
		{
			print $line."\n";
		}
		continue;
	}
	
	$parts = explode("\t", $line);
	list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info) = $parts;
	$info = explode(";", $info);
	
	//get AF
	$info_new = array();
	$af = null;
	$an = null;
	$hom = 0;
	$hemi = 0;
	foreach ($info as $entry)
	{
		//overall AF
		if (starts_with($entry, "AF="))
		{
			$af = substr($entry, 3);
		}
		else if (starts_with($entry, "AN="))
		{
			$an = substr($entry, 3);
			$info_new[] = $entry;
		}
		else if (starts_with($entry, "Hom="))
		{
			$info_new[] = $entry;
		}
		else if (starts_with($entry, "Hemi="))
		{
			$info_new[] = $entry;
		}
	}
	$info_new[] = 'AF='.($an<200 ? '0.0' : number_format($af, 5));

	print "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t".implode(";", $info_new)."\n";
}
fclose($handle);

?>