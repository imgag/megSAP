<?php
/** 
	@page db_add_mids
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("db_add_mids", "Adds MIDs to the NGSD database.");
$parser->addInfile("in",  "Input file in TSV format (name, sequence).", false);
$parser->addEnum("db",  "Database to connect to.", false, db_names());
extract($parser->parse($argv));

//init DB
$db_connect = DB::getInstance($db);
$db_connect->beginTransaction();
$hash = $db_connect->prepare("INSERT INTO mid (name, sequence) VALUES (:name, :sequence);");

$file = file($in);
foreach($file as $line)
{
	$line = trim($line);
	if ($line=="" || $line[0]=="#") continue;
	$parts = explode("\t", $line);
	if (count($parts)!=2)
	{
		$db->rollBack();
		trigger_error("Invalid input line: $line", E_USER_ERROR);
	}
	$db_connect->bind($hash, "name", $parts[0]);
	$db_connect->bind($hash, "sequence", $parts[1]);
	$db_connect->execute($hash, true);
}
$db_connect->endTransaction();
