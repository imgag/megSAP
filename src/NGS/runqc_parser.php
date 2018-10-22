<?php
/**
 * @page runqc_parser
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("runqc_parser", "Parses and imports run QC metrics into the NGSD.");
$parser->addString("name", "Name of the run.", false);
$parser->addString("run_dir", "Run directory.", false);
$parser->addFlag("force", "Overwrites already existing DB entries instead of throwing an error.");
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
extract($parser->parse($argv));

//init
$db_connect = DB::getInstance($db);
$run_dir = realpath($run_dir);

function get_col_value($cols, $index)
{
	$value = $cols[$index];
	if ($value=="nan") return NULL;
	
	$pos = strpos($value, " +/-");
	if ($pos!==FALSE)
	{
		$value = substr($value, 0, $pos);
	}
	
	return $value;
}

//check that flowcell ID is set in the NGSD
$result = $db_connect->executeQuery("SELECT id, fcid FROM sequencing_run WHERE name = '$name'");
if(count($result) != 1)	trigger_error("Sequencing run with the name $name not found in DB.", E_USER_ERROR);
$seq_run_id = $result[0]['id'];
$flowcell_id = $result[0]['fcid'];
if(!contains($run_dir, $flowcell_id))
{
	trigger_error("Flowcell ID from NGSD ($flowcell_id) was not found in directory '$run_dir'.", E_USER_ERROR);
}

if (!file_exists("Fq"))
{
	//parse number of cycles/index for each read
	$run_info_parsed = new SimpleXMLElement(file_get_contents($run_dir . "/RunInfo.xml"));
	$read_metrics = array();
	foreach ($run_info_parsed->{'Run'}->{'Reads'}->{'Read'} as $read)
	{
		$read_num = intval($read["Number"]);
		$read_metrics[$read_num]["cycles"] = intval($read["NumCycles"]);
		$read_metrics[$read_num]["index"] = $read["IsIndexedRead"] == "Y" ? "1" : "0";
	}

	//parse QC data from InterOp summary output
	list($stdout) = exec2(get_path("interop") . " $run_dir --level=3");
	foreach ($stdout as $line)
	{
		$line = trim($line);

		//split line
		$cols = str_getcsv($line);
		$cols = array_map('trim', $cols);

		//parse read number (for lane-wise QC)
		if (count($cols) == 1 && starts_with($cols[0], "Read "))
		{
			$read_num = $cols[0][5];
			continue;
		}
		if (count($cols) < 7) continue;

		//parse header
		if ($cols[0] == "Level" || $cols[0] == "Lane")
		{
			$indices = array_flip($cols);
			continue;
		}

		//parse content lines (read-wise)
		if (starts_with($cols[0], "Read "))
		{
			$read_num = $cols[0][5];
			$read_metrics[$read_num]["error_rate"] = get_col_value($cols, $indices['Error Rate']);
			$read_metrics[$read_num]["q30"] = get_col_value($cols, $indices['%>=Q30']);
		}

		//parse content lines (lane-wise)
		if (is_numeric($cols[0]) && $cols[1] == "-")
		{
			$lane = intval($cols[0]);
			$dens = get_col_value($cols, $indices['Density']) * 1000;
			$read_metrics[$read_num]['lane_metrics'][$lane]['cluster_dens'] = $dens;
			$read_metrics[$read_num]['lane_metrics'][$lane]['passed_filter'] = get_col_value($cols, $indices['Cluster PF']) / 100.0 * $dens;
			$read_metrics[$read_num]['lane_metrics'][$lane]['yield'] = get_col_value($cols, $indices['Yield']) * pow(1000, 3);
			$read_metrics[$read_num]['lane_metrics'][$lane]['error_rate'] = get_col_value($cols, $indices['Error']);
			$read_metrics[$read_num]['lane_metrics'][$lane]['q30'] = get_col_value($cols, $indices['%>=Q30']);
		}
	}
}
else
{
	//get number of cycles from BioInfo.csv
	$mgi_general = array_column(load_tsv("Fq/L01/BioInfo.csv", ","), 1, 0);
	$read_metrics[1]["cycles"] = intval($mgi_general["Read1 Cycles"]);
	$read_metrics[1]["index"] = "0";
	$read_metrics[2]["cycles"] = intval($mgi_general["Read2 Cycles"]);
	$read_metrics[2]["index"] = "0";
	$read_metrics[3]["cycles"] = intval($mgi_general["Barcode"]);
	$read_metrics[3]["index"] = "1";

	//parse summary files
	foreach (glob("Fq/L??/summaryTable.csv") as $f)
	{
		$lane = intval(substr(basename(dirname($f)), 1));
		$summary[$lane] = array_column(load_tsv($f, ","), 1, 0);
	}

	//parse Q30 and number of clusters
	$q30 = [];
	foreach (glob("Fq/L??/BasecallQC.txt") as $f)
	{
		$handle = fopen($f, "r");
		$lane = intval(substr(basename(dirname($f)), 1));
		$q30[1][$lane] = [];
		$q30[2][$lane] = [];
		$q30[3][$lane] = [];

		$cluster[$lane] = 0;

		$parse = "";
		while ( ($line = fgets($handle)) !== FALSE )
		{
			$line = trim($line);

			if ($parse === "Q30")
			{
				$values = explode("\t", $line);
				array_push($q30[1][$lane], array_sum(array_slice($values, 0, $read_metrics[1]["cycles"])));
				array_push($q30[2][$lane], array_sum(array_slice($values, $read_metrics[1]["cycles"], $read_metrics[2]["cycles"])));
				array_push($q30[3][$lane], array_sum(array_slice($values, $read_metrics[2]["cycles"], $read_metrics[3]["cycles"])));
				$parse = "";
			}
			elseif ($parse === "cluster")
			{
				$cluster[$lane] += intval($line);
				$parse = "";
			}

			if (starts_with($line, "#CycQ30"))
			{
				//per-cycle Q30 value
				$parse = "Q30";
			}
			elseif (starts_with($line, "#NUMDNBS"))
			{
				//number of DNBs
				$parse = "cluster";
			}
		}
	}

	foreach ($q30 as $read => $lane_q30_values)
	{
		$read_metrics[$read]["q30"] = 0;

		foreach ($lane_q30_values as $lane => $values)
		{
			//per lane and read
			$read_metrics[$read]["lane_metrics"][$lane]["q30"] = array_sum($values) / ($summary[$lane]["ImageArea"] * $read_metrics[$read]["cycles"]) * 100;
			$read_metrics[$read]["q30"] += $read_metrics[$read]["lane_metrics"][$lane]["q30"];

			$read_metrics[$read]["lane_metrics"][$lane]["cluster_dens"] = $cluster[$lane];
			$read_metrics[$read]["lane_metrics"][$lane]["passed_filter"] = floatval($summary[$lane]["TotalReads(M)"]) * 1e6;
			$read_metrics[$read]["lane_metrics"][$lane]["yield"] = floatval($summary[$lane]["TotalReads(M)"]) * 1e6 * $read_metrics[$read]["cycles"];
			$read_metrics[$read]["lane_metrics"][$lane]["error_rate"] = .0;
		}

		//per read
		$read_metrics[$read]["q30"] /= count($lane_q30_values);
		$read_metrics[$read]["error_rate"] = .0;
	}
}

//check if runqc data was already imported for this run
$result = $db_connect->executeQuery("SELECT id FROM runqc_read WHERE sequencing_run_id = {$seq_run_id}");
if(count($result)!=0)
{
	if ($force) //remove data
	{
		foreach($result as $row)
		{
			$db_connect->executeStmt("DELETE FROM runqc_lane WHERE runqc_read_id = ".$row['id']);
			$db_connect->executeStmt("DELETE FROM runqc_read WHERE id = ".$row['id']);
		}
	}
	else //abort
	{
		trigger_error("QC data for this run and read was already imported. Use the flag '-force' to overwrite it.", E_USER_ERROR);
	}
}

//import QC data into NGSD
foreach($read_metrics as $read_num => $read_metric)
{
	//insert read
	$hash = $db_connect->prepare("INSERT INTO runqc_read (read_num, cycles, is_index, q30_perc, error_rate, sequencing_run_id) VALUES (:read_num, :cycles, :is_index, :q30_perc, :error_rate, :seq_run_id);");
	$db_connect->bind($hash, "read_num", $read_num);
	$db_connect->bind($hash, "cycles", $read_metric['cycles']);
	$db_connect->bind($hash, "is_index", $read_metric['index']);
	$db_connect->bind($hash, "q30_perc", $read_metric['q30']);
	$db_connect->bind($hash, "error_rate", $read_metric['error_rate']);
	$db_connect->bind($hash, "seq_run_id", $seq_run_id);
	$db_connect->execute($hash, true);
	
	//ID of the inserted read
	$read_id = $db_connect->lastInsertId();
	
	//insert lane metrics for that read
	foreach($read_metric['lane_metrics'] as $lane_num => $lane_metric)
	{
		$hash = $db_connect->prepare("INSERT INTO runqc_lane (lane_num, cluster_density, cluster_density_pf, yield, error_rate, q30_perc, runqc_read_id) VALUES (:lane_num, :cl_dens, :cl_dens_pf, :yield, :error_rate, :q30_perc, :runqc_read_id);");
		$db_connect->bind($hash, "lane_num", $lane_num);
		$db_connect->bind($hash, "cl_dens", $lane_metric['cluster_dens']);
		$db_connect->bind($hash, "cl_dens_pf", $lane_metric['passed_filter']);
		$db_connect->bind($hash, "yield", $lane_metric['yield']);
		$db_connect->bind($hash, "error_rate", $lane_metric['error_rate']);
		$db_connect->bind($hash, "q30_perc", $lane_metric['q30']);
		$db_connect->bind($hash, "runqc_read_id", $read_id);
		$db_connect->execute($hash, true);
	}
}
?>