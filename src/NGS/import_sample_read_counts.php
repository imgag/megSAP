<?php
/** 
	@page import_sample_read_counts
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("import_sample_read_counts", "Parses the read counts of all samples from a given JSON stats file and imports it into the NGSD.");
$parser->addInfile("stats",  "JSON file containing the demultiplexing results.", false);
$parser->addFlag("csv_mode", "Use CSV instead of JSON stat file (NovaSeq X)");
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
extract($parser->parse($argv));

function read_novaseq_json_file($filename)
{
	$read_counts = array();
	$json_file_content = json_decode(file_get_contents($filename),true);

	// extract read counts
	// iterate over all lanes
	foreach($json_file_content["ConversionResults"] as $lane_info)
	{
		print "parsing lane ".$lane_info["LaneNumber"]." ... \n";
		// iterate over all samples
		foreach($lane_info["DemuxResults"] as $sample_info)
		{
			$sample_name = $sample_info["SampleName"];
			// * 2 since a read pair contains two reads
			$sample_read_count = $sample_info["NumberReads"] * 2;
			if (array_key_exists($sample_name, $read_counts))
			{
				// add read counts to read count from previous lanes
				$read_counts[$sample_name] = $read_counts[$sample_name] + $sample_read_count;
			}
			else
			{
				// create new entry for the sample
				$read_counts[$sample_name] = $sample_read_count;
			}
		}
	}
	
	return $read_counts;
}

function read_novaseqx_csv_file($filename)
{
	$read_counts = array();
	$rows = file($filename);
	foreach ($rows as $row) 
	{
		//skip header
		if (starts_with($row, "Lane,")) continue;
		//skip empty lines
		if (trim($row) == "") continue;

		//extract read counts
		$parts = explode(",", $row);
		$sample_name = trim($parts[1]);
		$n_reads = intval($parts[3]);
		//skip Undetermined
		if($sample_name == "Undetermined") continue;
		if (isset($read_counts[$sample_name]))
		{
			// * 2 since a read pair contains two reads
			$read_counts[$sample_name] += $n_reads * 2;
		}
		else
		{
			// * 2 since a read pair contains two reads
			$read_counts[$sample_name] = $n_reads * 2;
		}

	}
	return $read_counts;
}

// Abort if no connection to the NGSD is available
if (!db_is_enabled($db))
{
	trigger_error("No Connection to the NGSD available! Couldn't store read counts.", E_USER_ERROR);
}

// get id for read count in qc_terms
$db_conn = DB::getInstance($db);
$query = "SELECT id FROM qc_terms WHERE qcml_id = 'QC:2000005'";
$res = $db_conn->executeQuery($query);
$qc_term_id = $res[0]["id"];

// read JSON/CSV file
if (!file_exists($stats))
{
	trigger_error("Stats file not found in $stats! No read count import possible.", E_USER_ERROR);
}
if ($csv_mode)
{
	trigger_error("NOTICE: CSV mode (NovaSeq X)", E_USER_NOTICE);
	$read_counts = read_novaseqx_csv_file($stats);
	print "Reads of ".count($read_counts)." samples extracted.\n";
}
else
{
	trigger_error("NOTICE: JSON mode (e.g. NovaSeq 6000)", E_USER_NOTICE);
	$read_counts = read_novaseq_json_file($stats);
	print "Reads of ".count($read_counts)." samples extracted.\n";
}


// import read counts into the database
foreach($read_counts as $sample_name => $read_count)
{
	print "$sample_name: \t $read_count \n";
	// get psid
	$ps_id = get_processed_sample_id($db_conn, $sample_name);
	// generate query:
	$query = "REPLACE INTO processed_sample_qc (processed_sample_id, qc_terms_id, value) "
			."VALUES ($ps_id, $qc_term_id, $read_count)";
	$db_conn->executeStmt($query); 
}

?>