<?php
/** 
	@page db_import_qc
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("db_import_qc", "Imports QC terms to NGSD.");
$parser->addString("id", "Processing ID (e.g. GS000123_01).", false);
$parser->addInfileArray("files", "qcML files to import.", false, true);
//optional
$parser->addEnum("db",  "Database to connect to.", true, array("NGSD", "NGSD_TEST"), "NGSD");
$parser->addFlag("force", "Overwrites already existing DB entries instead of throwing an error.");
$parser->addInt("min_depth", "Minimum average read depth required.", true, "40");
$skip_qual_parameters = array();
$parser->addString("skip_parameters", "Comma-separated list of quality parameters that are not inserted into the Database.", true, implode(",", $skip_qual_parameters));
extract($parser->parse($argv));

//load valid term list from OBO file
load_qc_terms();
$terms = $GLOBALS["qcml"];

// Terms that are not uploaded
$skip_parameters = explode(",", $skip_parameters);

//check pid format
if(!preg_match("/^([A-Za-z0-9]{4,})_(\d{2})$/", $id, $matches))
{
	trigger_error("'$id' is not a valid processing ID!", E_USER_ERROR);
}
$sid = $matches[1];
$pid = (int)$matches[2];

//establish database connection
$db = DB::getInstance($db);

// check if this PID exists in the DB
$hash = $db->prepare("SELECT s.id as sid, ps.id as psid FROM sample as s, processed_sample as ps WHERE s.id = ps.sample_id and s.name = :sid and ps.process_id = :pid");
$db->bind($hash, 'sid', $sid);
$db->bind($hash, 'pid', $pid);
$db->execute($hash, true); 
$result = $db->fetch($hash); 
$db->unsetStmt($hash);
if(count($result) != 1)
{
	trigger_error("PID '$id' not found in DB. Sample and Processed Sample entered to DB?", E_USER_ERROR);
}
$psid = $result[0]['psid'];

// check if QC parameters were already imported for this pid
$result = $db->executeQuery("SELECT count(id) FROM processed_sample_qc WHERE processed_sample_id='$psid'");
$count_old = $result[0]['count(id)'];
$parser->log("Found $count_old QC terms for PID '$id' already in DB.");

//throw error if not forced to overwrite and terms already exist
if($count_old!=0 && !$force)
{
	trigger_error("QC terms were already imported for PID '$id'. Use '-force' to overwrite them.", E_USER_ERROR);
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

		//check term exists in the OBO file
		if (!isset($terms[$key]))
		{
			trigger_error("qcML term '$key' not found in OBO file!", E_USER_ERROR);
		}
		$term = $terms[$key];
		
		//skip numeric terms that are not available
		if ($term["type"]!="string" && starts_with($value, "n/a"))
		{
			continue;
		}
		
		//skip terms that should be skipped
		if(in_array($key, $skip_parameters)) continue;
		
		//error if term is found several times
		if(isset($qc_par[$key])) trigger_error("Found QC term '$key' more than once!", E_USER_ERROR);
		
		//check that value is of the correct type
		switch($term["type"])
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
				continue;
				break;
			default:
				trigger_error("Internal error: Unknown QC term type '$type'!", E_USER_ERROR);
		}
		
		$qc_par[$key] = $value;
	}
}
$parser->log("Found ".count($qc_par)." QC terms for PID '$id' in log files.");

// insert QC terms into the DB
$hash = $db->prepare("INSERT INTO processed_sample_qc (processed_sample_id, qc_terms_id, value) VALUES (:0, (SELECT id FROM qc_terms WHERE qcml_id = :1), :2);");
$db->beginTransaction();
foreach($qc_par as $key => $value)
{
	$db->bind($hash, "0", $psid);
	$db->bind($hash, "1", $key);
	$db->bind($hash, "2", $value);
	$db->execute($hash, true, false);
}
$db->endTransaction();
$db->unsetStmt($hash);

?>
