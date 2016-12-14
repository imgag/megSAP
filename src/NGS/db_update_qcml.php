<?php

/** 
	@page db_update_qcml
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("db_update_qcml", "Update qcML terms in NGSD.");
//optional
$parser->addEnum("db",  "Database to connect to.", true, array("NGSD", "NGSD_TEST"), "NGSD");
extract($parser->parse($argv));

// database connection
$db = DB::getInstance($db);
$db->executeStmt("SET SESSION sql_mode = 'STRICT_TRANS_TABLES';");
$db->beginTransaction();

// setup ontology
load_qc_terms();
$terms = $GLOBALS["qcml"];
foreach($terms as $id => $data)
{
	//skip terms for plots
	if ($data['type']=="base64Binary") continue;
	
	$name = addslashes($data['name']);
	$desc = addslashes($data['def']);
	$db->executeStmt("INSERT INTO qc_terms (qcml_id, name, description) VALUES ('$id','{$name}','{$desc}') ON DUPLICATE KEY UPDATE name='{$name}',description='{$desc}'");
}

//disconnect
$db->endTransaction();
unset($db);

?>