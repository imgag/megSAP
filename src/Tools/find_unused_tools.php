<?php
/** 
	@page find_unused_tools 
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("find_unused_tools", "Checks for unused ngs-bits tools.");
$parser->addInfile("ngsbits", "ngs-bits repository path.", false);
$parser->addInfile("megsap", "ngs-bits repository path.", false);
//optional
$parser->addInt("n", "Number of usages below which a tool is reported.", true, 1);
extract($parser->parse($argv));

//get list of ngs-bits tools
list($tools) = exec2("find {$ngsbits} -name \"tools.pro\" | xargs grep SUBDIRS | grep -v cpp | grep -v TEST | cut -f3 -d' '");

//get list of megSAP tools
list($tools2) = exec2("find {$megsap}/src/ -name \"*.php\" | grep -v /Common/ | grep -v /Pipelines/ | xargs -l1 basename");
$tools = array_merge($tools, $tools2);

//find hits
$scripts_file = $parser->tempFile(".txt", "find_unused_tools");
exec2("find {$megsap} /mnt/storage1/share/data/ /mnt/storage1/share/chips/ /mnt/storage1/share/doc/ /mnt/storage1/share/kasp/ /mnt/storage1/share/primer/ /mnt/storage3/users/ahsturm1/Sandbox/ /mnt/storage1/users/bioinf/http/ -name \"Makefile\" -or -name \"*.sh\" -or -name \"*.php\" -or -name \"*.py\"> {$scripts_file}");

foreach($tools as $tool)
{
	list($hits) = exec2("cat {$scripts_file} | xargs grep {$tool}", false);
	
	//filter hits
	$good_hits = array();
	foreach($hits as $hit)
	{
		if (trim($hit)=="") continue;
		
		//skip c++ build folder hits
		if (strpos($hit, "build-tools-Linux")!==false) continue;
		
		//skip php error messages
		if (strpos($hit, "E_USER_ERROR")!==false) continue;
		
		$good_hits[] = $hit;
	}
	
	//output
	if (count($good_hits)<$n)
	{
		print "$tool\t".count($good_hits)."\n";
		foreach($good_hits as $hit)
		{
			print "  HIT: ".$hit."\n";
		}
	}
}

?>