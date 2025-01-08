<?php
/** 
	@page find_unused_tools 
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("find_unused_tools", "Checks for unused ngs-bits and megSAP tools (execute on SRV005).");
$parser->addInfile("ngsbits", "ngs-bits repository path.", false);
$parser->addInfile("megsap", "ngs-bits repository path.", false);
//optional
$parser->addInt("n", "Number of usages below which a tool is reported.", true, 1);
extract($parser->parse($argv));

//get list of ngs-bits tools
list($tmp) = exec2("grep SUBDIRS {$ngsbits}/src/tools.pro | grep -v cpp | grep -v TEST | cut -f3 -d' '");
$tools = [];
foreach($tmp as $t) $tools[$t] = "ngs-bits";

//get list of megSAP tools
list($tmp) = exec2("find {$megsap}/src/Tools/ -name '*.php' | xargs -l1 basename");
foreach($tmp as $t) $tools[$t] = "megSAP";

###################### find unused tools ######################

//get scripts
$scripts_file = $parser->tempFile(".txt", "find_unused_tools");
exec2("find {$megsap} -name 'Makefile' -or -name '*.sh' -or -name '*.php' -or -name '*.py' -or -name '*.md' | grep -v tool_test_ > {$scripts_file}");
exec2("find /mnt/storage2/users/ahsturm1/scripts/ -maxdepth 2 -name '*.php' >> {$scripts_file}");
exec2("find /mnt/storage1/share/kasp/ -name '*.php' >> {$scripts_file}");
exec2("find /var/www/bioinf/http/ -name 'Makefile' -or -name '*.sh' -or -name '*.php' -or -name '*.py' | grep -v 'NightlyTests/tmp/' >> {$scripts_file}");

//get matches
foreach($tools as $tool => $type)
{
	list($hits) = exec2("cat {$scripts_file} | xargs grep {$tool}", false);

	//filter hits
	$good_hits = array();
	foreach($hits as $hit)
	{
		if (trim($hit)=="") continue;
		
		//skip php error messages
		if (strpos($hit, "E_USER_ERROR")!==false) continue;
		
		$good_hits[] = $hit;
	}
	
	//output
	if (count($good_hits)<$n)
	{
		print "$type\t$tool\t".count($good_hits)."\n";
		foreach($good_hits as $hit)
		{
			print "  HIT: ".$hit."\n";
		}
	}
}

###################### find unused function ######################

//find functions
list($functions) = exec2("egrep '^function ' {$megsap}/src/Common/functions.php {$megsap}/src/Common/genomics.php | cut -f2 -d' ' | cut -f1 -d'('");

//get scripts
$scripts_file = $parser->tempFile(".txt", "find_unused_tools");
exec2("find {$megsap} /mnt/storage1/share/kasp/ /var/www/bioinf/http/ -name '*.php' | grep -v 'NightlyTests/tmp/' > {$scripts_file}");
exec2("find /mnt/storage2/users/ahsturm1/scripts/ -maxdepth 2 -name '*.php' >> {$scripts_file}");

//get matches
foreach($functions as $function)
{
	list($hits) = exec2("cat {$scripts_file} | xargs grep '{$function}(' | grep -v 'function {$function}('", false);
	
	//filter hits
	$good_hits = array();
	foreach($hits as $hit)
	{
		if (trim($hit)=="") continue;
		
		//skip php error messages
		if (strpos($hit, "E_USER_ERROR")!==false) continue;
		
		$good_hits[] = $hit;
	}
	
	//output
	if (count($good_hits)<$n)
	{
		print "megSAP function\t$function\t".count($good_hits)."\n";
		foreach($good_hits as $hit)
		{
			print "  HIT: ".$hit."\n";
		}
	}
}
?>