<?php

$h = fopen("php://stdin", "r");
while(!feof($h))
{
	$line = trim(fgets($h));
	if ($line=="") continue;
	list($name, $chr, $start, $mapq) = explode("\t", $line);
	
	list($chr, $start, $end) = explode("_", $name);
	print $chr."\t".($start-1)."\t".$end."\t".$mapq."\n";
}


?>