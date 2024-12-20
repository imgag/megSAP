<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("compare_variants", "Wrapper for SampleDiff.");
$parser->addInfile("in1", "First GSvar file.", false, false);
$parser->addInfile("in2", "Second GSvar file.", false, false);
extract($parser->parse($argv));

//exec
list($stdout, $stderr) = execApptainer("ngs-bits", "SampleDiff", "-window 0 -in1 {$in1} -in2 {$in2} 2>&1", [$in1, $in2]);
foreach($stdout as $line)
{
	print $line."\n";
}
foreach($stderr as $line)
{
	print $line."\n";
}

?>