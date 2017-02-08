<?php
/** 
	@page find_unused_tools 
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("find_unused_tools", "Checks for unused ngs-bits tools.");
$parser->addInfile("path", "ngs-bits repository path.", false);
//optional
$parser->addInt("n", "Number of usages below which a tool is reported.", true, 3);
extract($parser->parse($argv));

//get list of tools
list($tools) = exec2("find $path -name \"tools.pro\" | xargs grep SUBDIRS | grep -v cpp | grep -v TEST | cut -f3 -d' '");

foreach($tools as $tool)
{
	//find hits
	list($hits) = exec2("find /mnt/share/data/ /mnt/share/chips/ /mnt/share/doc/ /mnt/share/kasp/ /mnt/share/primer/ /mnt/users/ahsturm1/ -name \"Makefile\" -or -name \"*.sh\" -or -name \"*.php\" 2> /dev/null | xargs grep $tool");
	
	//filter hits
	$good_hits = array();
	foreach($hits as $hit)
	{
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