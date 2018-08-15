<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");
error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("db_get_gender_genlab", "Returns gender of sample ID from Genlab.");
$parser->addOutFile("id","Sample ID (e.g DX152541_01).",false);
//optional
extract($parser->parse($argv));
if(!db_is_enabled("GL8"))
{
	trigger_error("Could not connect to genlab db.\n",E_USER_ERROR);
}

$db = DB::getInstance("GL8");

//Resolve patient ID from GenLab
$res = $db->executeQuery("SELECT PATIENT_ID FROM t_untersuchung WHERE labornummer = '{$id}'");
$res = array_unique($res);
if(empty($res))
{
	trigger_error("Sample ID {$id} not found in t_untersuchung.",E_USER_ERROR);
}
if(count($res) > 1)
{
	trigger_error("Found multiple patient IDs for sample ID {$id}",E_USER_ERROR);
}

$patient_id = $res[0]["PATIENT_ID"];

//Resolve Gender
$res2 = $db->executeQuery("SELECT GESCHLECHT FROM t_patient WHERE id = '{$patient_id}'");
if(empty($res2))
{
	trigger_error("No data found for patient id {$patient_id}.",E_USER_ERROR);
}

$gender_code = $res2[0]["GESCHLECHT"];
if($gender_code == 1) //1:female
{
	print("female");
}
elseif($gender_code == 2) //2:male
{
	print("male");
}
else
{
	print("n/a");
}
?>