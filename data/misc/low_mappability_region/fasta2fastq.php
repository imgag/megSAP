<?php

$h = fopen("php://stdin", "r");
while(!feof($h))
{
	$line = trim(fgets($h));
	if ($line=="") continue;
	
	if ($line[0]==">")
	{
		print "@".substr(strtr($line, ":-", "__"), 1)."\n";
	}
	else
	{
		print $line."\n";
		print "+\n";
		print str_repeat("I", strlen($line))."\n";
	}
}


?>