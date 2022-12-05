<?php
/** 
	@page db_import_variant_classes
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("db_import_variant_classes", "Imports variants with classification into the NGSD.");
$parser->addString("user",  "NGSD user used in the classification comment.", false);
$parser->addInfile("in", "Input VCF file.", false);
$parser->addString("text", "Classification text added to NGSD entry.", true, "");
//optional
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
extract($parser->parse($argv));

//init
$db = DB::getInstance($db);
$valid_chrs = $db->getEnum("variant", "chr");
$valid_classes = $db->getEnum("variant_classification", "class");
$hash_var = $db->prepare("INSERT INTO `variant`(`chr`, `start`, `end`, `ref`, `obs`) VALUES (:0, :1, :2, :3, :4, :5)");

//check input VCF
list($stderr, $stdout, $exit_code) = $parser->exec(get_path("ngs-bits")."VcfCheck", "-in {$in}", true, false, false);
if ($exit_code!=0)
{
	trigger_error("Invalid input VCF:\n".implode("\n", $stderr), E_USER_ERROR);
}

//check user id
$users = $db->getValues("SELECT user_id FROM user WHERE active='1' ORDER BY user_id");
if (!in_array($user, $users))
{
	trigger_error("Unknown NGSD user '$user'. Valid users are: '".implode("','", $users)."'.", E_USER_ERROR);
}

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
		trigger_error("Input VCF line with less than 7 columns found: ".$line, E_USER_ERROR);
	}
	list($chr, $pos, , $ref, $obs, , , $info) = $parts;
	
	//convert variant to GSvar format
	if (!in_array($chr, $valid_chrs))
	{
		trigger_error("Chromosome '$chr' is not valid. Valid chromosomes are '".implode("','", $valid_chrs)."'", E_USER_ERROR);
	}
	$start = trim($pos);
	$end = $start;
	$ref = strtoupper(trim($ref));
	$obs = strtoupper(trim($obs));
	if(strlen($ref)>1 || strlen($obs)>1) //correct indels
	{
		list($start, $end, $ref, $obs) = correct_indel($start, $ref, $obs);
	}
	
	//get classification
	$class = null;
	$info = explode(";", trim($info));
	foreach($info as $entry)
	{
		if (contains($entry, "="))
		{
			list($key, $value) = explode("=", $entry);
		}
		else
		{
			$key = $entry;
		}
		
		if ($key=="CLASS")
		{
			$class = $value;
		}
	}
	
	//check class is valid
	if (!in_array($class, $valid_classes))
	{
		trigger_error("Classification '$class' is not valid. Valid classes are '".implode("','", $valid_classes)."'", E_USER_ERROR);
	}
	
	//check if variant exists (insert otherwise)
	$v_id = $db->getValue("SELECT id FROM variant WHERE chr='{$chr}' AND start='{$start}' AND end='{$end}' AND ref='{$ref}' AND obs='{$obs}'", -1);
	if($v_id==-1)
	{
		$db->executeStmt("INSERT INTO `variant` SET chr='{$chr}', start='{$start}', end='{$end}', ref='{$ref}', obs='{$obs}'");
		$v_id = $db->lastInsertId();
	}
	
	//set classification
	$c_data = $db->executeQuery("SELECT id, comment, class FROM variant_classification WHERE variant_id='{$v_id}'");
	$comment = "[{$class}] {$user} ".date("Y-m-d")."\nBatch import from command line.\n{$text}\n";
	if (count($c_data)==0) //new classification
	{
		$db->executeStmt("INSERT INTO variant_classification SET variant_id='{$v_id}', class='{$class}', comment='{$comment}'");
		print "$chr:$start $ref>$obs\tAdded new classification {$class}\n";
	}
	else //classification exists => update
	{
		$c_id = $c_data[0]['id'];
		$c_comment = $db->quote($comment."\n".$c_data[0]['comment']);
		$c_class = $c_data[0]['class'];
		$db->executeStmt("UPDATE variant_classification SET class='{$class}', comment={$c_comment} WHERE id='{$c_id}'");
		print "$chr:$start $ref>$obs\tUpdated classification to {$class} (was {$c_class}) \n";
	}
}

?>