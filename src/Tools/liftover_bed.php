<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("liftover_bed", "Lift-over from.");
$parser->addInfile("in", "Input BED file.", false);
$parser->addOutfile("out", "Input BED file.", false);
$parser->addInt("max_size_dev", "Maximum size deviation.", true, 5);
$parser->addInt("max_size_dev_perc", "Maximum size deviation in percent.", true, 2);
extract($parser->parse($argv));

$tmp_in = temp_file(".bed");
$tmp_unmapped = temp_file(".bed");
$tmp_out = temp_file(".bed");

$ho = fopen($out, "w");
foreach(file($in) as $line)
{
	$line = nl_trim($line);
	if ($line=="") continue;
	
	//header
	if ($line[0]=="#")
	{
		fputs($ho, $line."\n");
		continue;
	}
	
	$parts = explode("\t", $line);
	list($chr, $start, $end) = $parts;

	//execute liftover tool
	file_put_contents($tmp_in, "$chr\t$start\t$end");
	list($stdout, $stderr, $exit) = exec2("/mnt/share/opt/liftOver/liftOver {$tmp_in} /mnt/share/opt/liftOver/GRCh37_to_GRCh38.chain.gz {$tmp_out} {$tmp_unmapped}");
	
	//check unmapped
	$unmapped = trim(file_get_contents($tmp_unmapped));
	if ($unmapped!="")
	{
		print "Error - not mappable: $line\n";
		continue;
	}
	
	//check size
	list($chr2, $start2, $end2) = explode("\t", trim(file_get_contents($tmp_out)));
	$size_before = $end-$start;
	$size_after = $end2-$start2;
	$size_diff = abs($size_before-$size_after);
	$size_diff_perc = 100.0 * $size_diff / $size_before;
	if ($size_diff>$max_size_dev && $size_diff_perc>$max_size_dev_perc)
	{
		print "Error - different size after mapping (from ".($end-$start)." to ".($end2-$start2)."): $line\n";
		continue;
	}
	
	//add prefix to chromosome
	if (!starts_with($chr2, "chr")) $chr2 = "chr".$chr2;
	
	//write output
	$parts[0] = $chr2;
	$parts[1] = $start2;
	$parts[2] = $end2;
	
	fputs($ho, implode("\t", $parts)."\n");
}

?>
