<?php

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
	if ($parts[0]=="MT")$parts[0]="M";

	if (isset($parts_last))
	{
		list($chr, $start, $dbsnp, $ref, $alt, , , $info) = $parts_last;
		list($chr2, $start2, $dbsnp2, $ref2, $alt2, , , $info2) = $parts;
		if ($chr==$chr2 && $start==$start2 && $ref==$ref2 && $alt==$alt2)
		{
			$info = explode(";dbSNPBuildID=", $info);
			$v1 = explode(";", $info[1], 2);
			$v1 = $v1[0];
			$info2 = explode(";dbSNPBuildID=", $info2);
			$v2 = explode(";", $info2[1], 2);
			$v2 = $v2[0];
			//print "1) $v1 ".implode(" ",$parts_last)."\n";
			//print "2) $v2 ".implode(" ",$parts)."\n";
			if ($v1>$v2)
			{
				$parts = $parts_last;
			}
			//print ">      ".implode(" ",$parts)."\n\n";
		}
		else
		{
			print implode("\t", $parts_last)."\n";
		}
	}
	$parts_last = $parts;
}
fclose($handle);
print implode("\t", $parts_last)."\n";

?>