<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("compare_variants", "Comparison of two GSvar files.");
$parser->addInfile("in1", "First GSvar file.", false, false);
$parser->addInfile("in2", "Second GSvar file.", false, false);
$parser->addString("add", "Comma-separated list of columns to compare in addition to chr,start,end,ref,obs.", true);
extract($parser->parse($argv));

//exec
$args = [];
$args[] = "-in1 {$in1}";
$args[] = "-in2 {$in2}";
$args[] = "-skip_comments";
$args[] = "-comp chr,start,end,ref,obs".($add=="" ? "" : ",".$add);
list($stdout, $stderr, $code) = execApptainer("ngs-bits", "TsvDiff", implode(" ", $args), [$in1, $in2], [], false, false, false);
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

exit($code);
?>