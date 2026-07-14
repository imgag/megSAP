<?php
/** 
	@page check_urls
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("check_urls", "Checks if URLS are broken in files.");
$parser->addInfileArray("files", "Input files.", false);
extract($parser->parse($argv));


//check for broken links
foreach($files as $file)
{
	$is_sh = ends_with($file, ".sh");
	$is_php = ends_with($file, ".php");
	if (!$is_sh && !$is_php) trigger_error("Unsupported file type for $file", E_USER_ERROR);
	
	foreach(file($file) as $line)
	{
		$line = trim($line);
		
		//skip empty and comment lines
		if ($line=="") continue;		
		if ($is_sh && starts_with($line, "#")) continue;
		if ($is_php && starts_with($line, "//")) continue;
		
		preg_match_all('#\b(https|http|ftp)://[^\s<>"\']+#', $line, $matches);
		if (count($matches[0])==0) continue;

		$url = $matches[0][0];
		
		print basename($file)." - $url\n";

		if(!url_exists($url))
		{
			print "\tMISSING!\n";
		}
	}
}

?>