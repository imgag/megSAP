<?php
/** 
	@page correct_file_rights
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("correct_file_rights", "Corrects file owner and permissions.");
$parser->addStringArray("folders", "Folders to change.", false);
$parser->addString("owner", "File owner.", true, "bioinf");
$parser->addString("group", "File group.", true, "f_ad_bi_l_medgen_access_storages");
$parser->addString("permissions", "File permissions.", true, "775");
extract($parser->parse($argv));

foreach($folders as $folder)
{
	$folder = trim($folder);
	
	//change permissions (before user and group to make sure the files that are currently being written are still writable)
	$time_start = microtime(true);
	print "Changing permissions of '{$folder}' to '{$permissions}'...\n";
	exec2("chmod -R {$permissions} {$folder}");
	print "  took ".time_readable(microtime(true)-$time_start)."\n";
		
	//change group (everyone should have this group)
	$time_start = microtime(true);
	print "Changing group of '{$folder}' to '{$group}'...\n";
	exec2("chgrp -R {$group} {$folder}");
	print "  took ".time_readable(microtime(true)-$time_start)."\n";
	
	//change user
	$time_start = microtime(true);
	print "Changing owner of '{$folder}' to '{$owner}'...\n";
	exec2("chown -R {$owner} {$folder}");
	print "  took ".time_readable(microtime(true)-$time_start)."\n";
	
}

?>