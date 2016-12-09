<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/all.php");

list($stdout, $stderr) = exec2(get_path("ngs-bits")."/SampleDiff -window 0 -ei -in1 ".$argv[1]." -in2 ".$argv[2]);

foreach($stdout as $line)
{
	print $line."\n";
}
foreach($stderr as $line)
{
	print $line."\n";
}

?>