<?php
/** 
	@page export_system 
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("export_system", "Exports a processing system INI file.");
$parser->addString("name", "Processing system short name.", false);
$parser->addOutfile("out", "Output file in TSV format.", false);
extract($parser->parse($argv));

//check output file does not exist
if (file_exists($out))
{
	trigger_error("Output file '$out' already exists!", E_USER_ERROR);
}

//get processing system id
$ps_id = get_processed_sample_name_by_processing_system($name);
if ($ps_id===false)
{
	trigger_error("Processing system short name '$name' unknown!", E_USER_ERROR);
}

//store system
$tmp = "";
load_system($tmp, $ps_id);

//move system from temp to output file
$parser->moveFile($tmp, $out);

?>