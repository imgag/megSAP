<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("db_import_variant_classifications", "\$Rev: 804 $", "Import variant classifications into the NGSD.");
$parser->addInfile("in",  "Input sample list in TSV format (chr, start, end, ref, obs, class, comment).", false);
$parser->addString("user",  "NGSD user that created the classification (for comment header).", false);
$parser->addEnum("db",  "Database to connect to.", true, array("NGSD", "NGSD_TEST"), "NGSD");
extract($parser->parse($argv));

//check user id
$db = DB::getInstance($db);
$users = $db->getValues("SELECT user_id FROM user WHERE active='1' ORDER BY user_id");
if (!in_array($user, $users))
{
	trigger_error("Unknown NGSD user '$user'. Valid users are: '".implode("','", $users)."'.", E_USER_ERROR);
}

//get valid classes
$valid_classes = $db->getEnum("variant_classification", "class");

//import data
$file = file($in);
foreach($file as $line)
{
	$line = nl_trim($line);
	if($line=="" || $line[0]=="#") continue;
	
	//split and check column count
	$parts = explode("\t", $line);
	if (count($parts)<7)
	{
		trigger_error("Input file line with less than 7 columns found: ".$line, E_USER_ERROR);
	}
	list($chr, $start, $end, $ref, $obs, $class, $comment) = $parts;
	
	//check variant exists
	$v_id = $db->getValue("SELECT id FROM variant WHERE chr='".$chr."' AND start='".$start."' AND end='".$end."' AND ref='".$ref."' AND obs='".$obs."'", -1);
	if($v_id==-1)
	{
		trigger_error("Variant $chr:$start $ref>$obs not in NGSD!", E_USER_ERROR);
	}
	
	//check enum is valid
	if (!in_array($class, $valid_classes))
	{
		trigger_error("Classification '$class' is not valid. Valid classes are '".implode("','", $valid_classes)."'", E_USER_ERROR);
	}
	
	//set classification
	$c_data = $db->executeQuery("SELECT id, comment, class FROM variant_classification WHERE variant_id='".$v_id."'");
	$comment = "[$class] $user ".date("Y-m-d")."\n$comment\n";
	if (count($c_data)==0) //new classification
	{
		$db->executeStmt("INSERT INTO variant_classification (variant_id, class, comment) VALUES ('$v_id', '$class', '$comment')");
		print "$chr:$start $ref>$obs\tAdded new classification $class\n";
	}
	else //classification exists => update
	{
		$c_id = $c_data[0]['id'];
		$c_comment = $c_data[0]['comment'];
		$c_class = $c_data[0]['class'];
		$db->executeStmt("UPDATE variant_classification SET class='$class', comment='$comment\n$c_comment' WHERE id='$c_id'");
		print "$chr:$start $ref>$obs\tUpdated classification to $class (was $c_class) \n";
	}
}

?>