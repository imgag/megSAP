<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

$handle = fopen("php://stdin", "r");
while($line = fgets($handle))
{
	$line = trim($line);
	if ($line=="") continue;
	if($line[0]=="#")
	{
		print $line."\n";
		continue;
	}
	
	$parts = explode("\t", $line);
	list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info) = $parts;
	$info = explode(";", $info);
	$info_new = array();
	$ac = null;
	$an = null;
	foreach ($info as $entry)
	{
		//keep
		if (starts_with($entry, "AC_Hom=") || starts_with($entry, "Hom_NFE=") || starts_with($entry, "Hom_AFR="))
		{
			$info_new[] = $entry;
		}
		//rename original AF
		if (starts_with($entry, "AF="))
		{
			$info_new[] = "OLD_".$entry;
		}
		//calcualte new AF
		if (starts_with($entry, "AC_Adj="))
		{
			$info_new[] = $entry;
			$ac = substr($entry, 7);
		}
		if (starts_with($entry, "AN_Adj="))
		{
			$info_new[] = $entry;
			$an = substr($entry, 7);
		}		
	}
	//print $line;
	$info_new[] = "AF=".(is_null($ac) || is_null($an) || $an<200 ? "0.0" : $ac/$an);
	print "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t".implode(";", $info_new)."\n";
}
fclose($handle);

?>