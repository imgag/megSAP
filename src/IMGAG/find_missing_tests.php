<?php
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

//init
$src = repository_basedir()."/src/";
$test = repository_basedir()."/test/";

//find functions
list($functions) = exec2("egrep '^function ' {$src}Common/functions.php {$src}Common/genomics.php | cut -f2 -d' ' | cut -f1 -d'('");

//find tests
list($tests) = exec2("ls $test | grep tool_test_");
$tests = array_map(function($value) { return substr($value, 10, -4); }, $tests);

//find main tools
$tools2path = array();
list($tmp) = exec2("find {$src}Tools/ -name \"*.php\" -or -name \"*.py\"");
sort($tmp);
foreach($tmp as $t)
{
	$tool = basename2($t);
	$tools2path[$tool] = substr($t, strlen($src));
}
list($tools_imgag) = exec2("find {$src}IMGAG/ {$src}Auxilary/ {$src}Deprecated/ -name \"*.php\" -or -name \"*.py\"");
$tools_imgag = array_map('basename2', $tools_imgag);

//(1) find left-over tool tests (tool already deleted)
$missing_tools = array_diff($tests, array_keys($tools2path));
$missing_tools = array_diff($missing_tools, $tools_imgag);
foreach($missing_tools as $tool)
{
	if ($tool=="barcode_correction") continue; //barcode_correction.by from umiVar has no wrapper in megSAP
	print "Test but no tool: $tool\n";
}

//(2) find unused test data
list($test_files) = exec2("ls test/data/");
foreach($test_files as $test_file)
{
	$test_match = false;
	foreach($tests as $test)
	{
		if (starts_with($test_file, $test)) $test_match = true;
	}
	foreach($functions as $function)
	{
		if (starts_with($test_file, $function)) $test_match = true;
	}
	foreach(['manta_sv_normal', 'manta_sv_tumor', 'somatic_pipeline_rna_counts', 'somatic_pipeline_tumor_rna'] as $other) //used in pipeline tests
	{
		if (starts_with($test_file, $other)) $test_match = true;
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
$missing_tests = array_diff(array_keys($tools2path), $tests);
foreach($missing_tests as $test)
{
	if($test=="annotate" || $test=="mapping" || $test=="data_setup") continue; //base functionality tested in pipeline tests
	if($test=="vc_dragen_somatic" || $test=="mapping_dragen") continue; //DRAGEN cannot be easily tested
	
	//count usage in megSAP
	$usage = array();
	$hits = array();
	exec("find {$src}Install/ {$src}Tools/ {$src}Pipelines/ {$src}Auxilary/ {$src}IMGAG/ -name '*.php' -or -name '*.py' | xargs egrep '{$test}.php|{$test}.py'", $hits);
	if (count($hits)>0)
	{
		$usage[] = count($hits)."x in megSAP";
	}
	
	//count usage in webservices
	$hits = array();
	exec("find /var/www/bioinf/http/ -name '*.php' -or -name '*.py' | grep -v '/+old/' | grep -v '/tmp/' | xargs egrep '{$test}.php|{$test}.py'", $hits);
	if (count($hits)>0)
	{
		$usage[] = count($hits)."x in webservices";
	}
	
	print "Missing test: ".$tools2path[$test];
	if (count($usage)!=0) print " (used: ".implode(", ", $usage).")";
	print "\n";
}

?>
