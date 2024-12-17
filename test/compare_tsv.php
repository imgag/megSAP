<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("compare_tsv", "Wrapper for TsvDiff.");
$parser->addInfile("in1", "First TSV file.", false, false);
$parser->addInfile("in2", "Second TSV file.", false, false);
$parser->addString("skip_cols", "Passed to TsvDiff.", true);
extract($parser->parse($argv));

//exec
$args = [];
$args[] = "-in1 {$in1}";
$args[] = "-in2 {$in2}";
if ($skip_cols!="") $args[] = "-skip_cols {$skip_cols}";

list($stdout, $stderr, $exit_code) = execApptainer("ngs-bits", "TsvDiff", implode(" ", $args), [$in1, $in2], [], false, false, false);
foreach($stdout as $line)
{
	print $line."\n";
}
foreach($stderr as $line)
{
	print $line."\n";
}

exit($exit_code);

?>