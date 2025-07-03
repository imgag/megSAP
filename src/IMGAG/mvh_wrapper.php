<?php
/** 
	@page mvh_wrapper
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("mvh_wrapper", "Wrapper for GRZ/KDK export.");
$parser->addInt("case_id", "'id' in 'data_data' of 'MVH' database.", false);
$parser->addEnum("type", "Export type", false, ["GRZ", "KDK_SE"]);
$parser->addFlag("clear", "Clear date (passed to export script).");
$parser->addFlag("test", "Test mode (passed to export script).");
extract($parser->parse($argv));

//execute export script
$args = [];
$args[] = "-case_id {$case_id}";
if ($clear) $args[] = "-clear";
if ($test) $args[] = "-test";
$script = dirname(realpath($_SERVER['SCRIPT_FILENAME']))."/mvh_".strtolower($type)."_export.php";
$command = "php {$script} ".implode(" ", $args)." 2>&1";
list($stdout, , $exit_code) = exec2($command);
	
//get submission IDs from output
$id_mvh = "";
$id_sub = "";
foreach($stdout as $line)
{
	$line = trim($line);
	
	if (starts_with($line, "SUBMISSION ID MVH DB:"))
	{
		$id_mvh = trim(explode(":", $line, 2)[1]);
	}
	
	if (starts_with($line, "SUBMISSION ID {$type}:"))
	{
		$id_sub = trim(explode(":", $line, 2)[1]);
	}
}

//update submission status in MVH database
$status = $exit_code==0 ? "done" : "failed";
$db_mvh = DB::getInstance("MVH");
$hash = $db_mvh->prepare("UPDATE submission_".strtolower($type)." SET status=:status, submission_id=:id_sub, output=:output WHERE id=:id");
$db_mvh->bind($hash, "status", $status);
$db_mvh->bind($hash, "submission_id", $id_sub);
$db_mvh->bind($hash, "output", $output);
$db_mvh->bind($hash, "id", $id);
$db_mvh->execute($hash, true);

?>
