<?php
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

//init
$src = repository_basedir()."/src/";
$test = repository_basedir()."/test/";

//find tests
list($tests) = exec2("ls $test | grep tool_test_");
$tests = array_map(function($value) { return substr($value, 10, -4); }, $tests);

//find tools under version control
$tools2path = array();
list($tmp) = exec2("find {$src}NGS/ {$src}Tools/ {$src}Primer/ -name \"*.php\" -or -name \"*.py\"");
sort($tmp);
foreach($tmp as $t)
{
	$tool = basename($t, ".php");
	if (ends_with($tool, ".py")) $tool = substr($tool, 0, -3);
	$tools2path[$tool] = substr($t, strlen($src));
}

//(1) find missing tools
$missing_tools = array_diff($tests, array_keys($tools2path));
foreach($missing_tools as $tool)
{
	print "Test but no tool: $tool\n";
}

//(2) find unused test data
list($test_files) = exec2("ls test/data/");
foreach($test_files as $test_file)
{
	$test_match = false;
	foreach($tests as $test)
	{
		if (starts_with($test_file, $test))
		{
			$test_match = true;
			break;
		}
	}
	if (!$test_match)
	{
		list($matches) = exec2("grep {$test_file} test/*.php || true");
		if (trim(implode("", $matches))=="")
		{
			print "Unused test file: test/data/{$test_file}\n";
		}
	}
}

//(3) find missing tests
$excluded_patterns = array(
				  "qbic_copy", //needs datamover
				  
				  "db_converter_.*", //converters to set up annotation databases
				  "find_unused_tools", //utility
				  "backup_.*", //backup scripts
				  "low_cov_regions", //utility for low-coverage statistics
				  "check_bcl", //utility
				  "find_missing_tests", //utility
				  "find_wrong_adapters", //utility
				  );
$missing_tests = array_diff(array_keys($tools2path), $tests);
foreach($missing_tests as $test)
{
	//skip excluded tools
	$skip = false;
	foreach($excluded_patterns as $pattern)
	{
		if(preg_match('/^'.$pattern.'$/', $test)) $skip = true;
	}
	if ($skip) continue;
	
	//count usage in megSAP
	$usage = array();
	$hits = array();
	exec("find {$src}NGS/ {$src}Tools/ {$src}Primer/ {$src}Pipelines/ -name '*.php' -or -name '*.py' | xargs grep $test", $hits);
	if (count($hits)>0)
	{
		$usage[] = count($hits)."x in megSAP";
	}
	
	//count usge in webservices
	$hits = array();
	exec("find /mnt/storage1/users/bioinf/http/ -name '*.php' -or -name '*.py' | grep -v '/+old/' | grep -v '/tmp/' | xargs grep $test", $hits);
	if (count($hits)>0)
	{
		$usage[] = count($hits)."x in webservices";
	}
	
	print "Missing test: ".$tools2path[$test];
	if (count($usage)!=0) print " (used: ".implode(", ", $usage).")";
	print "\n";
}


?>
