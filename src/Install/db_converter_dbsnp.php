<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

$handle = fopen2("php://stdin", "r");
while($line = fgets($handle))
{
	$line = nl_trim($line);
	if ($line=='') continue;
	
	//headers
	if($line[0]=='#')
	{
		if (starts_with($line, "#CHROM\t") || starts_with($line, '##fileformat=') || starts_with($line, '##fileDate=') || starts_with($line, '##dbSNP_BUILD_ID='))
		{
			print $line."\n";
		}
		if (starts_with($line, '##INFO=<ID=RS,'))
		{
			print "##INFO=<ID=RS,Number=1,Type=String,Description=\"dbSNP ID\">\n";
		}
		
		continue;
	}
	
	//variants
	$parts = explode("\t", $line);
	
	//get chromosome
	$chr = chr2NC($parts[0], true, false);
	if ($chr=='') continue; //non-standard chromosome
	$parts[0] = $chr;
	
	//set RS number (from ID)
	$parts[7] = 'RS='.$parts[2];
	
	//we need the RS number on in INFO
	$parts[2] = '.';
	
	print implode("\t", $parts)."\n";
}


?>