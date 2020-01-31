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
extract($parser->parse($argv));

//Returns a description of the database in TSV format
function get_desc($db_name)
{
	//init
	$db = DB::getInstance($db_name);
	$tables = $db->getValues("SHOW TABLES");
	sort($tables);
	$description_fields = array("Field","Type","Null","Key","Default","Extra");
	
	$output = array();	
	foreach($tables as $table)
	{
		$res = $db->executeQuery("DESCRIBE $table");
		foreach($res as $row)
		{
			$field_name = $row["Field"];
			$tag = $table."/".$field_name;
			
			//DESCRIBE
			foreach($description_fields as $description_field)
			{
				$value = $row[$description_field];
				//ignore order of enums
				if ($description_field=="Type" && starts_with($value, "enum("))
				{
					$enum = substr($row["Type"], 6, -2);
					$enum = explode("','", $enum);
					sort($enum);	
					$value = "enum('".implode("','", $enum)."')";
				}
				$output[$tag][$description_field] = $value;
			}
			
			//comment
			$output[$tag]["Comment"] = $db->getValue("SELECT COLUMN_COMMENT FROM INFORMATION_SCHEMA.COLUMNS WHERE TABLE_SCHEMA=database() AND TABLE_NAME='{$table}' AND COLUMN_NAME='{$field_name}'");			
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
$db_prod = get_desc("NGSD");
$db_test = get_desc("NGSD_TEST");

//determine all tags
$tags = array_merge(array_keys($db_prod), array_keys($db_test));
sort($tags);
$tags = array_unique($tags);

//search for additional/missing fields
foreach($tags as $tag)
{
	if (!isset($db_prod[$tag]))
	{
		print "Missing table field in production database: $tag\n";
	}
	else if (!isset($db_test[$tag]))
	{
		print "Extra table field in production database: $tag\n";
	}
	else
	{
		$keys = array_keys($db_prod[$tag]);
		foreach($keys as $key)
		{
			if ($db_prod[$tag][$key] != $db_test[$tag][$key])
			{
				print "Difference in {$tag} in property {$key}: is '".$db_prod[$tag][$key]."', but should be '".$db_test[$tag][$key]."'\n";
			}
		}
	}
}

?>