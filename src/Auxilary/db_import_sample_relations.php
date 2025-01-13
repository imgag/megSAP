<?php
/** 
	@page db_import_sample_relations
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("db_import_sample_relations", "Imports sample relations (batch).");
$parser->addInfile("in", "TSV input file (sample1, relation, sample2).", false);
//optional
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
extract($parser->parse($argv));

function get_sample_id($sample_name)
{
	global $db;
	
	$res = $db->executeQuery("SELECT id FROM sample WHERE name=:name", array('name' => $sample_name));
	if (count($res)!=1)
	{
		trigger_error("Sample name not found: ".$sample_name, E_USER_WARNING);
		$sampleid = "";
	}
	else
	{
		$sampleid = $res[0]['id'];
	}
	
	# return: sample found, sample id
	return [count($res)===1, $sampleid];
}

function is_contained($id1, $rel, $id2)
{
	global $db;
	
	$res = $db->executeQuery("SELECT id FROM sample_relations WHERE sample1_id=:id1 AND relation=:rel AND sample2_id=:id2", array('id1'=>$id1, 'rel'=>$rel, 'id2'=>$id2));
	
	return count($res)>=1;
}

//init
$db = DB::getInstance($db);
$valid_relations = 	$db->getEnum("sample_relations", "relation");

$file = file($in);
foreach($file as $line)
{
	$line = trim($line);
	if ($line=="" || $line[0]=="#") continue;
	
	//split columns
	$parts = explode("\t", $line);
	if (count($parts)<3) trigger_error("Invalid line (less than 3 tab-separated columns): ".$line, E_USER_ERROR);
	list($s1, $rel, $s2) = $parts;
		
	//check input data
	list($sample1_found, $id1) = get_sample_id($s1);
	list($sample2_found, $id2) = get_sample_id($s2);
	if (!$sample1_found || !$sample2_found)
	{
		continue;
	}
	if (!in_array($rel, $valid_relations)) trigger_error("Invalid sample relation '$rel' in line : ".$line, E_USER_ERROR);
	
	//check if already contained
	$is_commutative = ($rel=="sample sample") || ($rel=="siblings");
	if (is_contained($id1, $rel, $id2) || ($is_commutative && is_contained($id2, $rel, $id1)))
	{
		print "Skipping '$s1 $rel $s2' - already contained!\n";
		continue;
	}
	
	print "Adding '$s1 $rel $s2'...\n";
	$db->executeStmt("INSERT INTO sample_relations (`sample1_id`, `relation`, `sample2_id`) VALUES (:id1, :rel, :id2)", array('id1'=>$id1, 'rel'=>$rel, 'id2'=>$id2));
}

?>