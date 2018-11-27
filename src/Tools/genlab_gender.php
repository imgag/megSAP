<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");
error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("genlab_gender", "Returns gender of sample ID from Genlab.");
$parser->addOutFile("id","Sample ID (e.g DX152541_01).",false);
//optional
extract($parser->parse($argv));
if(!db_is_enabled("GL8"))
{
	trigger_error("Could not connect to GenLab database.", E_USER_ERROR);
}

$db = DB::getInstance("GL8");

//Resolve patient ID from GenLab
$patient_ids = $db->getValues("SELECT PATIENT_ID FROM t_untersuchung WHERE labornummer = '{$id}'");
$patient_ids = array_unique($patient_ids);
if(count($patient_ids)!=1)
{
	trigger_error("Found no/multiple patient IDs for sample ID {$id}: ".implode(", ", $patient_ids), E_USER_ERROR);
}
$patient_id = $patient_ids[0];

//Resolve Gender
$res2 = $db->executeQuery("SELECT GESCHLECHT FROM t_patient WHERE id = '{$patient_id}'");
if(empty($res2))
{
	trigger_error("No gender in GenLab for patient ID {$patient_id}.", E_USER_ERROR);
}

$gender_code = $res2[0]["GESCHLECHT"];
if($gender_code == 1) //1:female
{
	print("female\n");
}
elseif($gender_code == 2) //2:male
{
	print("male\n");
}
else
{
	print("n/a\n");
}
?>