<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/all.php");

function seg2tsv($filename)
{
	$output = [];
	$file = file($filename);
	foreach($file as $line)
	{
		$line = trim($line);
		if ($line=="" || $line[0]=="#") continue;
				
		if (starts_with($line, "ID\tchr\tstart\tend"))
		{
			$output[] = "#".$line."\n";
		}
		else
		{
			$output[] = $line."\n";
		}
	}
	
	return $output;
}

//store data to tmp files
$file1 = temp_file(".tsv", "compare_seg");
file_put_contents($file1, seg2tsv($argv[1]));
$file2 = temp_file(".tsv", "compare_seg");
file_put_contents($file2, seg2tsv($argv[2]));


//exec
$args = [];
$args[] = "-in1 {$file1}";
$args[] = "-in2 {$file2}";
if (isset($argv[3]))
{
	$args[] = "-diff_abs ".$argv[3];
}

list($stdout, $stderr, $exit_code) = execApptainer("ngs-bits", "TsvDiff", implode(" ", $args), [$file1, $file2], [], false, false, false);
foreach($stdout as $line)
{
	if ($line=="") continue;
	print $line."\n";
}
foreach($stderr as $line)
{
	if ($line=="") continue;
	print $line."\n";
}

exit($exit_code);
?>