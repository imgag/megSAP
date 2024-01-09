<?php
/** 
	@page correct_file_rights
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("correct_file_rights", "Corrects file owner and permissions.");
$parser->addString("types", "Comma-separated list of project types.", false);
$parser->addString("owner", "File owner.", true, "bioinf");
$parser->addString("group", "File group.", true, "domänen-benutzer");
$parser->addString("permissions", "File permissions.", true, "775");
extract($parser->parse($argv));

//init
$db = DB::getInstance("NGSD");
$project_folders  = get_path("project_folder");

//check types
$project_types = explode(",", $types);
$project_types = array_map('trim', $project_types);
$valid_types = $db->getEnum("project", "type");
foreach($project_types as $type)
{
	if (!in_array($type, $valid_types)) trigger_error("Invalid type '$type'!", E_USER_ERROR);
}

foreach($project_types as $type)
{
	$folder = $project_folders[$type];
	
	$time_start = microtime(true);
	print "Changing owner of '{$folder}' to '{$owner}'...\n";
	exec2("chown -R {$owner} {$folder}");
	print "  took ".time_readable(microtime(true)-$time_start)."\n";
	
	$time_start = microtime(true);
	print "Changing group of '{$folder}' to '{$group}'...\n";
	exec2("chgrp -R {$group} {$folder}");
	print "  took ".time_readable(microtime(true)-$time_start)."\n";
	
	$time_start = microtime(true);
	print "Changing permissions of '{$folder}' to '{$permissions}'...\n";
	exec2("chmod -R {$permissions} {$folder}");
	print "  took ".time_readable(microtime(true)-$time_start)."\n";
}

?>