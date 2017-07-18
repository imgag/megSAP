<?php 
/** 
	@page copy_sample
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("test", "test tool");
extract($parser->parse($argv));

?>