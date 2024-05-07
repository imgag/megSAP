<?php
/** 
	@page export_file 
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("export_file", "Exports a file (single file of tar.");
$parser->addInfile("file", "File to export.", false);
$parser->addFlag("internal", "Use internal webserver.");
extract($parser->parse($argv));

//init
if (!file_exists($file))
{
	trigger_error("File does not exist: {$file}", E_USER_ERROR);
}

//determine password and folder
$share_url = $internal ? "https://datashare.img.med.uni-tuebingen.de/" : "https://download.imgag.de/DataShare/";
list($stdout) = exec2("curl --noproxy '*' ".($internal ? " -k" : "")." '{$share_url}/index.php?action=request&filename={$file}'");
if (contains(implode(" ", $stdout), "ERROR:")) trigger_error(implode(" ", $stdout), E_USER_ERROR);
$folder = trim($stdout[0]);
$password = trim($stdout[1]);
 
//move
print "You can move the file to the webserver using:\n";
if($internal)
{
	print "  > mv {$file} /mnt/storage1/share/http_shareukt/DataShare/data/{$folder}/\n";
}
else
{
	print "  > scp {$file} imgag.de:/var/www/html/download/DataShare/data/{$folder}/\n";
}

print "The URL of the file is:\n";
print "  {$share_url}/index.php?filename={$file}\n";

print "Password:\n";
print "  {$password}\n";

?>