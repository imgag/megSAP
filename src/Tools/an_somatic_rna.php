<?php 
/** 
	@page an_somatic_rna
	
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("an_somatic_rna", "RNA annotation for somatic pipeline.");


extract($parser->parse($argv));


?>
