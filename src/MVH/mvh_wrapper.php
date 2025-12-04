<?php
/** 
	@page mvh_wrapper
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("mvh_wrapper", "Wrapper for GRZ/KDK export.");
$parser->addInt("cm_id", "ID in case management RedCap database.", false);
$parser->addEnum("type", "Export type", false, ["GRZ", "KDK_SE"]);
$parser->addFlag("clear", "Clear date (passed to export script).");
$parser->addFlag("test", "Test mode (passed to export script).");
extract($parser->parse($argv));

//init
$mvh_folder = realpath(get_path("mvh_folder"));

//execute export script
print "Performing export\n";
$args = [];
$args[] = "-cm_id {$cm_id}";
if ($clear) $args[] = "-clear";
if ($test) $args[] = "-test";
$args[] = "--log {$mvh_folder}/".strtolower($type)."_export/logs/{$cm_id}.log";
$script = dirname(realpath($_SERVER['SCRIPT_FILENAME']))."/mvh_".strtolower($type)."_export.php";
$command = "php {$script} ".implode(" ", $args)." 2>&1";
list($stdout, , $exit_code) = exec2($command, false);

//print output to command line
print "\n";
print "Export output:\n";
foreach($stdout as $line)
{
	$line = nl_trim($line);
	if ($line=="") continue;
	print "  {$line}\n";
}
print "Export exit code: $exit_code\n";

//get submission IDs from output (available for GRZ only)
$id_sub = "";
if ($type=="GRZ" && $exit_code==0)
{
	foreach($stdout as $line)
	{
		$line = trim($line);
		if (starts_with($line, "SUBMISSION ID GRZ:"))
		{
			$id_sub = trim(explode(":", $line, 2)[1]);
		}
	}
	print "\n";
	print "Submission ID of GRZ: {$id_sub}\n";
}

//get ID for submission status in MVH database
$id_status = "";
foreach($stdout as $line)
{
	$line = trim($line);
	if (starts_with($line, "ID in submission_".strtolower($type)." table:"))
	{
		$id_status = trim(explode(":", $line, 2)[1]);
	}
}
print "\n";
print "Submission ID in MVH database: {$id_status}\n";

//update submission status in MVH database
print "\n";
print "Updating MVH database\n";
if (!$test)
{
	$status = $exit_code==0 ? "done" : "failed";
	$db_mvh = DB::getInstance("MVH");
	$hash = $db_mvh->prepare("UPDATE submission_".strtolower($type)." SET status=:status, submission_id=:submission_id, submission_output=:submission_output, metadata=:metadata, date=CURDATE() WHERE id=:id");
	$db_mvh->bind($hash, "status", $status);
	$db_mvh->bind($hash, "submission_id", $id_sub);
	$db_mvh->bind($hash, "submission_output", implode("\n", $stdout));
	$db_mvh->bind($hash, "metadata", file_get_contents("{$mvh_folder}/metadata_archive/{$type}/{$cm_id}.json"));
	$db_mvh->bind($hash, "id", $id_status);
	$db_mvh->execute($hash, true);
}
else
{
	print "  Skipped because this is a test run\n";
}

?>
