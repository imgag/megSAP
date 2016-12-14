<?php
/** 
	@page db_diff 
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("db_diff", "Checks for schema differences between test and productive NGSD.");
//optional
$parser->addFlag("init", "Freshly initialize test database.");
$parser->addFlag("local", "Create database TSV in current directory instead of in /tmp/ folder.");
extract($parser->parse($argv));

//Returns a description of the database in TSV format
function get_desc($db_name)
{
	//init
	$db = DB::getInstance($db_name);
	$tables = $db->getValues("SHOW TABLES");
	sort($tables);
	$fields = array("Field","Type","Null","Key","Default","Extra");
	
	$output = array();
	$output[] = "#Table\t".implode("\t", $fields)."\n";
	
	foreach($tables as $table)
	{
		$res = $db->executeQuery("DESCRIBE $table");
		foreach($res as $row)
		{
			$line = "$table";
			foreach($fields as $field)
			{
				$line .= "\t".$row[$field];
			}
			$line .= "\n";
			$output[] = $line;
		}
	}
	return $output;
}

//init
if ($init)
{
	exec2(get_path("ngs-bits")."NGSDInit -test");
}

//diff columns
if ($local)
{
	$file_p = "ngsd_prod.tsv";
	$file_t = "ngsd_test.tsv";
}
else
{
	$file_p = $parser->tempFile("_prod.tsv");
	$file_t = $parser->tempFile("_test.tsv");
}
file_put_contents($file_t, get_desc("NGSD"));
file_put_contents($file_p, get_desc("NGSD_TEST"));
exec("diff $file_p $file_t", $diff);
foreach($diff as $line)
{
	print $line."\n";
}


?>