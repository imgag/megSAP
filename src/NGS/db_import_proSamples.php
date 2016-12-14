<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("db_import_proSamples", "Batch-imports processed samples to the NGSD.");
$parser->addInfile("in",  "Input sample list in TSV format.", false);
extract($parser->parse($argv));

function get_sample_id($sample_name)
{
	$db_connect = DB::getInstance("NGSD");
	$result_array = $db_connect->executeQuery("SELECT id FROM sample WHERE name ='".$sample_name."'");
	if (count($result_array)==1)
	{
		return $result_array[0]['id'];
	}
	else
	{
		print "ERROR: sample '$sample_name' does not exist";
		return false;
	}
}

function get_run_id($run_number)
{
	$db_connect = DB::getInstance("NGSD");
	$result_array = $db_connect->executeQuery("SELECT id FROM sequencing_run WHERE name ='".$run_number."'");
	if (count($result_array)==1)
	{
		return $result_array[0]['id'];
	}
	else
	{
		print "ERROR: run '$run_number' does not exist";
		return false;
	}
}


function get_mid_id($name)
{
	$db_connect = DB::getInstance("NGSD");
	$result_array = $db_connect->executeQuery("SELECT id FROM mid WHERE name ='".$name."'");
	if (count($result_array)==1)
	{
		return $result_array[0]['id'];
	}
	else
	{
		print "ERROR: mid '$name' does not exist";
		return false;
	}
}

function get_operator_id($name)
{
	$db_connect = DB::getInstance("NGSD");
	$result_array = $db_connect->executeQuery("SELECT id FROM user WHERE name ='".$name."'");
	if (count($result_array)==1)
	{
		return $result_array[0]['id'];
	}
	else
	{
		print "WARNING: operator '$name' does not exist. Set to NULL";
		return null;
	}
}

function get_pro_sys_id($name)
{
	$db_connect = DB::getInstance("NGSD");
	$result_array = $db_connect->executeQuery("SELECT id FROM processing_system WHERE name_manufacturer ='".$name."'");
	if (count($result_array)==1)
	{
		return $result_array[0]["id"];
	}
	else
	{
		print "ERROR: processing system '$name' does not exist";
		return false;
	}
}

function get_project_id($name)
{
	$db_connect = DB::getInstance("NGSD");
	$result_array = $db_connect->executeQuery("SELECT id FROM project WHERE name ='".$name."'");
	if (count($result_array)==1)
	{
		return $result_array[0]['id'];
	}
	else
	{
		print "ERROR: project '$name' does not exist";
		return false;
	}
}

function check_mid($name,$sequence)
{
	$db_connect = DB::getInstance("NGSD");
	$result_array = $db_connect->executeQuery("SELECT sequence FROM mid WHERE name ='".$name."'");
	if (count($result_array)!=1)
	{
		print "ERROR: Mid $name not found in db!\n";
		return false;
	}
	if (trim($result_array[0]['sequence'])==$sequence)
	{
		return true;
	}
	else
	{
		print "ERROR: Given mid sequence '$sequence' for MID $name does not match db entry '".trim($result_array[0]['sequence'])."\n";
		return false;
	}
}

function generate_pro_id($sample_id)
{
	$db_connect = DB::getInstance("NGSD");
	$result_array = $db_connect->executeQuery("SELECT process_id FROM processed_sample WHERE sample_id ='".$sample_id."'");
	$process_ids=array(0);
	foreach($result_array as $result)
	{
		$process_ids[]=$result["process_id"];
	}
	return (max($process_ids)+1);
}

function get_normal_sample_id()
{
	//Todo: find normal sample
	print "WARNING: Functionality to add a normal sample to a tumor sample not added yet. Information NOT added.";
}

function add_pro_sample($sample_name,$project,$run,$lane,$mid1_i7_name,$mid1_i7_seq,$mid2_i5_name,$mid2_i5_seq,$normal_sample,$operator,$pro_sys,$molarity,$status,$comment)
{
	$mid1_i7_id=null;
	if($mid1_i7_name!="")
	{
		if ($mid1_i7_seq!="")
		{
			if (!(check_mid($mid1_i7_name,$mid1_i7_seq)))
			{
				return false;
			}
		}
		$mid1_i7_id=get_mid_id($mid1_i7_name);
	}
	$mid2_i5_id=null;
	if($mid2_i5_name!="")
	{
		if ($mid2_i5_seq!="")
		{
			if (!(check_mid($mid2_i5_name,$mid2_i5_seq)))
			{
				return false;
			}
		}
		$mid2_i5_id=get_mid_id($mid2_i5_name);
	}
	$sample_id=get_sample_id($sample_name);
	$run_id=get_run_id($run);
	$operator_id=get_operator_id($operator);
	$pro_sys_id=get_pro_sys_id($pro_sys);
	$project_id=get_project_id($project);
	$pro_id=generate_pro_id($sample_id);
	if ((!($sample_id))||(!($run_id))||(!($pro_sys_id))||(!($project_id)))
	{
		return false;
	}

	$db_connect = DB::getInstance("NGSD");
	$hash = $db_connect->prepare("INSERT INTO processed_sample (sample_id,process_id,sequencing_run_id,lane,mid1_i7,mid2_i5,operator_id,processing_system_id,project_id,molarity,comment) VALUES (:ps_sample_id,:ps_process_id,:ps_sequencing_run_id,:ps_lane,:ps_mid1_i7,:ps_mid2_i5,:ps_operator_id,:ps_processing_system_id,:ps_project_id,:ps_molarity,:ps_comment);");

	$db_connect->bind($hash,"ps_sample_id",$sample_id);
	$db_connect->bind($hash,"ps_process_id",$pro_id);
	$db_connect->bind($hash,"ps_sequencing_run_id",$run_id);
	$db_connect->bind($hash,"ps_lane",$lane);
	$db_connect->bind($hash,"ps_mid1_i7",$mid1_i7_id);
	$db_connect->bind($hash,"ps_mid2_i5",$mid2_i5_id);
	$db_connect->bind($hash,"ps_operator_id",$operator_id);
	$db_connect->bind($hash,"ps_processing_system_id",$pro_sys_id);
	$db_connect->bind($hash,"ps_project_id",$project_id);
	$db_connect->bind($hash,"ps_molarity",floatval($molarity));
	$db_connect->bind($hash,"ps_comment",$comment);
	
	$db_connect->execute($hash, true, false);

	return true;
}

$proc_sample_lines = file($in);
foreach($proc_sample_lines as $proc_sample_line)
{
	list($sample_name,$project,$run,$lane,$mid1_i7_name,$mid1_i7_seq,$mid2_i5_name,$mid2_i5_seq,$normal_sample,$operator,$pro_sys,$molarity,$status,$comment)=explode("\t",trim($proc_sample_line," \n\r\0\x0B"));//trim, but keep trailing tabs
	if ($molarity=="") $molarity=null;
	else $molarity=floatval($molarity);
	
	if(add_pro_sample($sample_name,$project,$run,$lane,$mid1_i7_name,$mid1_i7_seq,$mid2_i5_name,$mid2_i5_seq,$normal_sample,$operator,$pro_sys,$molarity,$status,$comment))
	{
		print "ADDED Processed sample to $sample_name\n";
	}
}
?>