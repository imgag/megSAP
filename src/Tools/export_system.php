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

//get processed sample that was processed with this processing system
$db_conn = DB::getInstance('NGSD');
$results = $db_conn->executeQuery("SELECT CONCAT(s.name,'_',LPAD(ps.process_id,2,'0')) as ps_id FROM processed_sample as ps, processing_system as sys, sample as s WHERE ps.sample_id = s.id AND ps.processing_system_id = sys.id AND sys.name_short = :system LIMIT 1", array('system' => $name));
if(count($results)==0)
{
	trigger_error("Could not find processing system with short name '$name' in NGSD!", E_USER_ERROR);
}

//store system
$tmp = "";
load_system($tmp, $results[0]['ps_id']);

//move system from temp to output file
$parser->moveFile($tmp, $out);

?>