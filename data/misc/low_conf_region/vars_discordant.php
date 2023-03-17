<?php

include("/mnt/storage2/megSAP/pipeline/src/Common/all.php");

$file = file($argv[1]);
foreach($file as $line)
{
	$line = trim($line);
	if ($line=="" || $line[0]=="#") continue;
	
	list($ps1, $ps2, $system, $path) = explode("\t", $line);
	print "{$path} \n";
	$out_disc = [];
	$c_vars = 0;
	$c_disc = 0;
	$i_ps1 = -1;
	$i_ps2 = -1;
	$h = gzopen2($path."/all.vcf.gz", "r");
	while(!gzeof($h))
	{
		$line =  trim(gzgets($h));
		if ($line=="") continue;
		
		$parts = explode("\t", $line);
		if ($line[0]=="#")
		{
			$out_disc[] = $line."\n";
			
			if ($line[1]!="#")
			{
				$i_ps1 = array_search($ps1, $parts);
				if ($i_ps1===false) $i_ps1 = array_search(substr($ps1, 0, -3), $parts);
				$i_ps2 = array_search($ps2, $parts);
				if ($i_ps2===false) $i_ps2 = array_search(substr($ps2, 0, -3), $parts);
				//print "  indices: $i_ps1 / $i_ps2\n";
			}
			
			continue;
		}
		
		$chr = $parts[0];
		if ($chr=="chrY" || $chr=="chrMT") continue;
		
		++$c_vars;
		$gt1 = strtr(explode(":", $parts[$i_ps1])[0], ["|"=>"/", "."=>"0"]);
		$gt2 = strtr(explode(":", $parts[$i_ps2])[0], ["|"=>"/", "."=>"0"]);
		if ($gt1!=$gt2)
		{
			++$c_disc;
			$out_disc[] = $line."\n";
		}
	}
	gzclose($h);
	print "  discordant: $c_disc/$c_vars (".number_format(100.0*$c_disc/$c_vars, 2)."%)\n";
	
	file_put_contents("twins/{$ps1}_{$ps2}_disc.vcf", $out_disc);
}

?>