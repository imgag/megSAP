<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/all.php");

$from = $argv[1];
if ($argc>2)
{
	$to = $argv[2];
}
else
{
	$to = ".";
}

$copy_file = false; //copy file instead of linking
if ($argc>3)
{
	if ($argv[3] == "-copy") $copy_file = true;
	else trigger_error("Invalid parameter '".$argv[3]."' provided. Only valid parameter is '-copy'!", E_USER_ERROR);
}

if ($copy_file)
{
	list($stdout, $stderr) = exec2("cp -f ".get_path("test_data_folder")."{$from} {$to}");

}
else
{
	list($stdout, $stderr) = exec2("ln -fs ".get_path("test_data_folder")."{$from} {$to}");

}
 
?>