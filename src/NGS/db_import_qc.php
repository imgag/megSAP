<?php
/** 
	@page db_import_qc
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("db_import_qc", "Imports QC terms to NGSD.");
$parser->addString("id", "Processes sample identifier (e.g. GS000123_01).", false);
$parser->addInfileArray("files", "qcML files to import.", false, true);
//optional
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
$parser->addFlag("force", "Overwrites already existing DB entries instead of throwing an error.");
extract($parser->parse($argv));

//establish database connection
$db = DB::getInstance($db);
$psid = get_processed_sample_id($db, $id);

// check if QC parameters were already imported for this pid
$count_old =  $db->getValue("SELECT count(id) FROM processed_sample_qc WHERE processed_sample_id='$psid'");
$parser->log("Found $count_old QC terms for processed sample '$id' already in DB.");

//throw error if not forced to overwrite and terms already exist
if($count_old!=0 && !$force)
{
	trigger_error("QC terms were already imported for processed sample '$id'. Use '-force' to overwrite them.", E_USER_ERROR);
}

//remove old terms
if ($count_old!=0 && $force)
{
	$db->executeStmt("DELETE FROM processed_sample_qc WHERE processed_sample_id='$psid'");
	$parser->log("Deleted old QC terms because the flag '-force' was used.");
}

// read QC terms from qcML files
$qc_par = array();
foreach($files as $file)
{
	$xdoc = new DomDocument;
	$xdoc->LoadXML(file_get_contents($file));
	$enties = $xdoc->GetElementsByTagName('qualityParameter');
	foreach($enties as $entry)
	{
		$key = $entry->getAttribute('accession');
		$value = $entry->getAttribute('value');

		//check term exists in the NGSD
		$type = $db->getValue("SELECT type FROM qc_terms WHERE qcml_id='$key'", "");
		if ($type=="")
		{
			trigger_error("QC term '$key' not found in NGSD > skipped!", E_USER_NOTICE);
			continue;
		}
		
		//skip numeric terms that are not available
		if ($type!="string" && starts_with($value, "n/a"))
		{
			continue;
		}
				
		//error if term is found several times
		if(isset($qc_par[$key])) trigger_error("Found QC term '$key' more than once!", E_USER_ERROR);
		
		//check that value is of the correct type
		switch($type)
		{
			case "int":
				if(!is_int((int)$value))
				{
					trigger_error("Value of '$key' is not an integer '$value'!", E_USER_ERROR);
				}
				$qc_par[$key] = (int)$value;
				break;
			case "float":
				if(!is_float((float)$value))
				{
					trigger_error("Value of '$key' is not a float '$value'!", E_USER_ERROR);
				}
				$qc_par[$key] = (float)$value;
				break;
			case "string":
				if(!is_string($value))
				{
					trigger_error("Value of '$key' is not a string '$value'!", E_USER_ERROR);
				}
				break;
			case "base64Binary": //skip plots
				continue(2);
				break;
			default:
				trigger_error("Internal error: Unknown QC term type '$type'!", E_USER_ERROR);
		}
		
		$qc_par[$key] = $value;
	}
}
$parser->log("Found ".count($qc_par)." QC terms for processed sample '$id'.");

// insert QC terms into the DB
$hash = $db->prepare("INSERT INTO processed_sample_qc (processed_sample_id, qc_terms_id, value) VALUES (:0, :1, :2)");
$db->beginTransaction();
foreach($qc_par as $key => $value)
{
	$db->bind($hash, "0", $psid);
	$qc_id = $db->getValue("SELECT id FROM qc_terms WHERE qcml_id='$key'");
	$db->bind($hash, "1", $qc_id);
	$db->bind($hash, "2", $value);
	$db->execute($hash, true);
}
$db->endTransaction();
$db->unsetStmt($hash);

?>
