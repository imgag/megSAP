<?php
/** 
	@page db_export
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("db_export", "Exports NGSD query as TSV file.");
$parser->addString("query", "Query.", false);
extract($parser->parse($argv));

if (!starts_with($query, "SELECT")) trigger_error("Query has to start with SELECT!", E_USER_ERROR);
if (contains($query, ";")) trigger_error("Query must not contain ';'!", E_USER_ERROR);

//execute query
$db = DB::getInstance("NGSD");
$res = $db->executeQuery($query);

if (count($res)==0)
{
	print "No rows selected!\n";
	exit(0);
}

//print header
$cols = array_keys($res[0]);
print "#".implode("\t", $cols)."\n";


//print contents
foreach($res as $row)
{
	for($c=0; $c<count($cols); ++$c)
	{
		$col = $cols[$c];
		if ($c>0) print "\t";
		print $row[$col];
	}
	print "\n";
}

?>