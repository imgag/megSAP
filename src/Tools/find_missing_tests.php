<?php
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

//init
$src = repository_basedir()."/src/";
$test = repository_basedir()."/test/";

//find tests
list($tests) = exec2("ls $test | grep tool_test_");
$tests = array_map(function($value) { return substr($value, 10); }, $tests);

//find tools under SVN control
$tools = array();
$tools2path = array();;
list($tmp) = exec2("find {$src}NGS/ {$src}Tools/ {$src}Primer/ -name \"*.php\"");
sort($tmp);
foreach($tmp as $t)
{
	$tool = basename($t);
	$tools[] = $tool;
	$tools2path[$tool] = substr($t, strlen($src));
}

//find missing tools
$missing_tools = array_diff($tests, $tools);
foreach($missing_tools as $tool)
{
	print "Missing tool: $tool\n";
}


//find missing tests
$excluded_patterns = array(
				  "queue_sample.php", //needs SGE
				  "queue_trio.php", //needs SGE
				  "queue_multi.php", //needs SGE
				  "qbic_copy.php", //needs datamover
				  
				  "db_converter_.*", //converters to set up annotation databases
				  "find_unused_tools.php", //utility
				  "backup_run.php", //utility
				  "check_bcl.php", //utility
				  "find_missing_tests.php", //utility
				  "find_wrong_adapters.php", //utility
				  );
$missing_tests = array_diff($tools, $tests);
foreach($missing_tests as $test)
{
	//skip excluded tools
	$skip = false;
	foreach($excluded_patterns as $pattern)
	{
		if(preg_match('/^'.$pattern.'$/', $test)) $skip = true;
	}
	if ($skip) continue;
	
	//count usage in php
	$usage = array();
	$hits = array();
	exec("find {$src}NGS/ {$src}Tools/ {$src}Primer/ {$src}Pipelines/ -name \"*.php\" | xargs grep $test", $hits);
	if (count($hits)>0)
	{
		$usage[] = count($hits)."x in php";
	}
	
	//count usge in DB
	$hits = array();
	exec("find /mnt/users/ahsturm1/SVN/DB/ -name \"*.php\" | xargs grep $test", $hits);
	if (count($hits)>0)
	{
		$usage[] = count($hits)."x in DB";
	}
	
	//count usge in webservices
	$hits = array();
	exec("find /mnt/users/ahsturm1/SVN/webservices/ -name \"*.php\" | xargs grep $test", $hits);
	if (count($hits)>0)
	{
		$usage[] = count($hits)."x in webservices";
	}
	
	print "Missing test: ".$tools2path[$test];
	if (count($usage)!=0) print " (used: ".implode(", ", $usage).")";
	print "\n";
}


?>
