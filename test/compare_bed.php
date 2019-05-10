<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/all.php");

function load_bed_coords($filename)
{
	$output = array();
	$file = file($filename);
	foreach($file as $line)
	{
		$line = trim($line);
		if ($line=="" || $line[0]=="#") continue;
		list($chr, $start, $end) = explode("\t", $line);	
		$output[] = "$chr\t$start\t$end\n";
	}
	
	return $output;
}

//store data to tmp files
$file1 = temp_file(".bed", "compare_bed");
file_put_contents($file1, load_bed_coords($argv[1]));
$file2 = temp_file(".bed", "compare_bed");
file_put_contents($file2, load_bed_coords($argv[2]));

//diff
list($stdout, $stderr, $code) = exec2("diff -b $file1 $file2", false);
foreach($stdout as $line)
{
	print $line."\n";
}
foreach($stderr as $line)
{
	print $line."\n";
}

exit($code);
?>