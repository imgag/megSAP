<?php
/**
 * @page runqc_parser
 */

$basedir = dirname($_SERVER['SCRIPT_FILENAME'])."/../";

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("runqc_parser", "\$Rev: 2$", "Parsers Illumina InterOp files and imports RunQC metrics into the NGSD database.");
$parser->addString("name", "Name of the run.", false);
$parser->addString("run_dir", "Absolute path to run directory.", false);
$parser->addFlag("force", "Overwrites already existing DB entries instead of throwing an error.");
$parser->addEnum("db",  "Database to connect to.", true, array("NGSD", "NGSD_TEST"), "NGSD_TEST");
extract($parser->parse($argv));



print "=================\n";
print basename($run_dir);
print "=================\n";
$db_connect = DB::getInstance($db);
$hash = $db_connect->prepare("SELECT id, fcid FROM sequencing_run WHERE name = :name");
$db_connect->bind($hash, 'name', $name);
$db_connect->execute($hash, true); 
$result = $db_connect->fetch($hash);
if(count($result) != 1)	trigger_error("Sequencing run with the name $name not found in DB.", E_USER_ERROR);
$seq_run_id = $result[0]['id'];
$flowcell_id = $result[0]['fcid'];
if(!(contains($run_dir,$flowcell_id))) trigger_error("Flowcell ID from DB ($flowcell_id) was not found in directory name.", E_USER_ERROR);

// Check if the following two xml-files are present
// - RunInfo.xml
// - RunParameters.xml
$run_info = $run_dir."/RunInfo.xml";
$run_parameters = $run_dir."/runParameters.xml";
if(!is_readable($run_info))
{
	trigger_error("Input directory '$run_dir' does not contain the file RunInfo.xml", E_USER_ERROR);
}
if(!is_readable($run_parameters)) {
	// For some sequencer the file has a capital RTA
	$run_parameters = $run_dir."RunParameters.xml";
	if(!is_readable($run_parameters)) {
		trigger_error("Input directory '$run_dir' does not contain the file RunParameters.xml or runParameters.xml", E_USER_ERROR);
	}
}

// Read and parse the two xml-files
$run_info_string = file_get_contents($run_info);
$run_parameters_string = file_get_contents($run_parameters);
$run_info_parsed = new SimpleXMLElement($run_info_string);
$run_parameters_parsed = new SimpleXMLElement($run_parameters_string);
# The contents of the InterOp file is described in this pdf-document:
# http://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/sav/sequencing-analysis-viewer-v1_8_46-guide-15066069-a.pdf
// Get the sequencer and RTA version from RunParameters.xml
// Type		Name		RTAVersion		MetricsFileName		file version
// NextSeq	JAMESBOND	2.4.11			QMetricsOut.bin		5
// HiSeq	MORPHEUS	1.18.64			QMetricsOut.bin		5
// MiSeq	NEO			1.18.54			QMetricsOut.bin		4
// MiSeq	TRINITY		1.18.54			QMetricsOut.bin		4
$sequencer = explode(' ', $run_parameters_parsed->Setup->ApplicationName)[0];
// Get the RTA version
if($sequencer == "HiSeq") {
	$RTA_version = $run_parameters_parsed->Setup->RTAVersion;
} else {
	$RTA_version = $run_parameters_parsed->RTAVersion;
}
print "Sequencer:\t".$sequencer."\n";
print "RTAVersion:\t".$RTA_version."\n";

// Check if all InterOp files are present
// Format of the quality metrics file name (based on the RTAversion)
$quality_metrics = $run_dir."/InterOp/QMetricsOut.bin";
if (version_compare($RTA_version, '1.18.64', ">=")) {
    $quality_metrics_version = 5;
} else {
	$quality_metrics_version = 4;
}
$error_metrics = $run_dir."/InterOp/ErrorMetricsOut.bin";
$tile_metrics = $run_dir."/InterOp/TileMetricsOut.bin";
if(!is_readable($quality_metrics)) trigger_error("Input directory '$run_dir' does not contain the InterOp/QMetricsOut.bin or InterOp/QualityMetricsOut.bin", E_USER_ERROR);
if(!is_readable($tile_metrics)) trigger_error("Input directory '$run_dir' does not contain the file InterOp/TileMetricsOut.bin", E_USER_ERROR);

# Get information about the reads (number of cycles for each read, indexed read?)
$read_metrics = array();
$num_reads = 0;
foreach ($run_info_parsed->{'Run'}->{'Reads'}->{'Read'} as $read) {
	$num_reads += 1;
	$read_num = intval($read["Number"]);
	$read_metrics[$read_num]["num_cycles"] = intval($read["NumCycles"]);
	$index_parsed = (string) $read["IsIndexedRead"];
	if($index_parsed == "Y") {
		$indexed = "1";
	}
	else {
		$indexed = "0";
	}
	$read_metrics[$read_num]["indexed"] = $indexed;
}

// Parse TileMetricsOut.bin that contains cluster density, cluster density (PF), number of clusters, and number of clusters (PF)
$tile_metrics_handle = fopen($tile_metrics, "rb");
// Parse the file version number
$data = fread($tile_metrics_handle, 1);
$tile_metrics_file_version = unpack("c", $data)[1];
// Parse the length of each record
$data = fread($tile_metrics_handle, 1);
$tile_metrics_file_record_length = unpack("c", $data)[1];

// Parse all records
// TileMetrics data will be stored as two dimensional array of lanes and tiles. Each entry is an associative array with the following keys:
$tile_metrics_keys = array();
$tile_metrics_keys[100] = 0; // cluster_density     K/mm^2
$tile_metrics_keys[101] = 1; // cluster_density_pf  K/mm^2
$tile_metrics_keys[102] = 2; // number_of_clusters
$tile_metrics_keys[103] = 3; // number_of_clusters_pf
$tile_metrics_data = array();

while(!feof($tile_metrics_handle)) {
	// 2 bytes lane number (uint16)
	$data = fread($tile_metrics_handle, 2);
	if(feof($tile_metrics_handle)) break;	// feof returns true if one has tried to read from an empty file
	$cur_lane = unpack("S", $data)[1];
	// 2 bytes tile number (uint16)
	$data = fread($tile_metrics_handle, 2);
	$cur_tile = unpack("S", $data)[1];
	// 2 bytes metric code (uint16)
	$data = fread($tile_metrics_handle, 2);
	$cur_metric_code = unpack("S", $data)[1];
	// 4 bytes metric value (float)
	$data = fread($tile_metrics_handle, 4);
	$cur_metric_value = unpack("f", $data)[1];
	
	// Save the read record (do not consider phasing or percent aligned reads)
	if($cur_metric_code < 200) {
		$tile_metrics_data[$cur_lane][$cur_tile][$tile_metrics_keys[$cur_metric_code]] = $cur_metric_value;
	}
}
fclose($tile_metrics_handle);

// Compute the cluster density, cluster density (PF), and yield per lane
// Each entry in the tile_metrics_data contains the data of a whole lane
$lane_metrics = array();
foreach($tile_metrics_data as $lane => $lane_data) {
	// The lane quality metrics
	$cluster_dens = 0;
	$cluster_dens_pf = 0;
	$num_cluster = 0;
	$num_cluster_pf = 0;
	$num_tiles = 0;
	// Iterate over the tiles of a lane
	foreach($lane_data as $tile_data) {
		$num_tiles += 1;
		$cluster_dens += $tile_data[0];
		$cluster_dens_pf += $tile_data[1];
		$num_cluster += $tile_data[2];
		$num_cluster_pf += $tile_data[3];
	}
	// Normalize cluster density by number of tiles
	$cluster_dens /= $num_tiles;
	$cluster_dens_pf /= $num_tiles;
	// Set the computed metrics for each read for the current lane
	for($i=1; $i<=$num_reads; $i++) {
		$lane_metrics[$i][$lane]["cluster_dens"] = $cluster_dens;
		$lane_metrics[$i][$lane]["cluster_dens_pf"] = $cluster_dens_pf;
		$lane_metrics[$i][$lane]["yield"] = ($num_cluster_pf * $read_metrics[$i]["num_cycles"]);
		$lane_metrics[$i][$lane]["num_cluster"] = $num_cluster;
		$lane_metrics[$i][$lane]["num_cluster_pf"] = $num_cluster_pf;
	}
}

// Parse QMetricsOut.bin/QualityMetricsOut.bin that contains Q1 up to Q50 for each lane/tile/cycle
$quality_metrics_handle = fopen($quality_metrics, "rb");
// Skip the file version number and the length of each record
$data = fread($quality_metrics_handle, 2);

// File version 5 needs a check if score binning was on
if($quality_metrics_version == 5) {
	
	// Was score binning on?
	$data = fread($quality_metrics_handle, 1);
	$score_binning = unpack("c", $data)[1];
	
	// If score_binning was on, skip the following  bytes
	if($score_binning) {
		// Parse the number of quality score bins
		$data = fread($quality_metrics_handle, 1);
		$num_bins = unpack("c", $data)[1];
		// From the number of bins one can compute how many bytes should be skippde
		$bytes_to_skip = 3 * $num_bins;
		$data = fread($quality_metrics_handle, $bytes_to_skip);
	}
}

// Parsing the quality scores is the same for both file versions
$quality_metrics_data = array();
while(!feof($quality_metrics_handle)) {
	// 2 bytes lane number (uint16)
	$data = fread($quality_metrics_handle, 2);
	if(feof($quality_metrics_handle)) break;	// feof returns true if one has tried to read from an empty file
	$cur_lane = unpack("S", $data)[1];
	// 2 bytes tile number (uint16)
	$data = fread($quality_metrics_handle, 2);
	$cur_tile = unpack("S", $data)[1];
	// 2 bytes cycle number (uint16)
	$data = fread($quality_metrics_handle, 2);
	$cur_cycle = unpack("S", $data)[1];
	// 4 * 50 bytes number of clusters assigned score (uint32) Q1 through Q50
	$data = fread($quality_metrics_handle, 200);
	$cur_quality_scores = unpack("L*", $data);
	
	// Save the quality scores
	$quality_metrics_data[$cur_lane][$cur_tile][$cur_cycle] = $cur_quality_scores;
}
fclose($quality_metrics_handle);

// Compute the %bases>Q30 for each read and each lane
$last_cycle = 1;
for($cur_read=1; $cur_read<=$num_reads; $cur_read++) {
	$first_cycle = $last_cycle;												// first cycle of the current read
	$last_cycle = $first_cycle + $read_metrics[$cur_read]["num_cycles"];	// last cycle of the current read		
	// Iterate over the lanes
	foreach($quality_metrics_data as $lane => $lane_data) {
		// Iterate over the tiles of the lane and the cycles that belong the the actual read
		$over_q30 = 0;															// number of clusters >Q30
		$less_or_equal_q30 = 0;     											// number of clusters <=Q30
		foreach($lane_data as $tile_data) {
			for($cycle=$first_cycle; $cycle<$last_cycle; $cycle++) {
				for($i=1; $i<30; $i++) {
					if(isset($tile_data[$cycle][$i]))
						$less_or_equal_q30 += $tile_data[$cycle][$i];
				}
				for($i=30; $i<=50; $i++) {
					if(isset($tile_data[$cycle][$i]))
						$over_q30 += $tile_data[$cycle][$i];
				}
			}
		}
	
		// Compute the percentage of bases >Q30
		$total = $less_or_equal_q30 + $over_q30;
		if($total != 0)
			$q30 = $over_q30 / $total;
		else
			$q30 = 0;
		// Set the computed metrics for the current lane
		$lane_metrics[$cur_read][$lane]["q30"] = $q30;
	}
}

// Parsing the error rates
if(is_readable($error_metrics))//skip if error metrics file is missing (e.g. because of no PhiX spike in)
{
	$error_metrics_handle = fopen($error_metrics, "rb");
	// Skip the file version number and the length of each record
	$data = fread($error_metrics_handle, 2);
	$error_metrics_data = array();
	while(!feof($error_metrics_handle)) {
		// 2 bytes lane number (uint16)
		$data = fread($error_metrics_handle, 2);
		if(feof($error_metrics_handle)) break;	// feof returns true if one has tried to read from an empty file
		$cur_lane = unpack("S", $data)[1];
		// 2 bytes tile number (uint16)
		$data = fread($error_metrics_handle, 2);
		$cur_tile = unpack("S", $data)[1];
		// 2 bytes cycle number (uint16)
		$data = fread($error_metrics_handle, 2);
		$cur_cycle = unpack("S", $data)[1];
		// 4 bytes error rate (float)
		$data = fread($error_metrics_handle, 4);
		$cur_error_rate = unpack("f", $data)[1];
		// 5 * 4 bytes data to skip (#perfect reads, reads with one, two, three, four errors) (uint32)
		$data = fread($error_metrics_handle, 20);
		// Save the quality scores
		$error_metrics_data[$cur_lane][$cur_tile][$cur_cycle] = $cur_error_rate;
	}
	fclose($error_metrics_handle);

	// Compute the mean error rates for each read and each lane
	$last_cycle = 1;
	for($cur_read=1; $cur_read<=$num_reads; $cur_read++) {
		$first_cycle = $last_cycle;												// first cycle of the current read
		$last_cycle = $first_cycle + $read_metrics[$cur_read]["num_cycles"];	// last cycle of the current read		
		// Iterate over the lanes
		foreach($error_metrics_data as $lane => $lane_data) {
			// Iterate over the tiles of the lane and the cycles that belong the the actual read
			$error_rate_sum = 0;
			$divisor = 0;
			foreach($lane_data as $tile_data) {
				for($cycle=$first_cycle; $cycle<($last_cycle-1); $cycle++) {
					if(isset($tile_data[$cycle])) {
						$error_rate_sum += $tile_data[$cycle];
						$divisor += 1;
					}
				}
			}
		
			// Compute mean error rate. Test if there are error rates given for at least one cycle to prevent division by 0. 
			if($divisor != 0)
				$error_rate = $error_rate_sum / $divisor;
			else
				$error_rate = 0;
			// Set the computed mean error rate for the current lane
			$lane_metrics[$cur_read][$lane]["error_rate"] = $error_rate;
		}
	}
}


// Append the lane metrics to the read_metrics and compute Q30 and mean error rate
for($i=1; $i<=$num_reads; $i++) {
	$read_metrics[$i]["lane_metrics"] = $lane_metrics[$i];
	// Q30 and mean error rate
	$q30 = 0;
	$error_rate = null;
	$total_cluster = 0;
	// Iterate with the foreach such that the number of lanes is known (there is no variable holding the number of lanes)
	foreach($read_metrics[$i]["lane_metrics"] as $lane => $lane_data) {
		$q30 += $read_metrics[$i]["lane_metrics"][$lane]["q30"] * $read_metrics[$i]["lane_metrics"][$lane]["num_cluster_pf"];
		if(array_key_exists("error_rate",$read_metrics[$i]["lane_metrics"][$lane]))
		{
			$error_rate += $read_metrics[$i]["lane_metrics"][$lane]["error_rate"] * $read_metrics[$i]["lane_metrics"][$lane]["num_cluster_pf"];
		}
		$total_cluster += $read_metrics[$i]["lane_metrics"][$lane]["num_cluster_pf"];
	}
	// Set the Q30 and the error rate
	$q30 = $q30 / $total_cluster;
	if ($error_rate!==null)
	{
		$error_rate = $error_rate / $total_cluster;
	}
	$read_metrics[$i]["q30"] = $q30;
	$read_metrics[$i]["error_rate"] = $error_rate;
}

// TODO just for testing
//print_r($read_metrics);

// Check if runqc data was already imported for this run
$hash = $db_connect->prepare("SELECT id  FROM runqc_read WHERE sequencing_run_id = :seq_run_id");
$db_connect->bind($hash, 'seq_run_id', $seq_run_id);
$db_connect->execute($hash, true); 
$result = $db_connect->fetch($hash);
if(count($result) != 0 && !$force) {
	trigger_error("QC for this run and read was already imported. Use the flag '-force' to overwrite them.", E_USER_ERROR);
}

// Remove old runqc
if(count($result) != 0 && $force) {
	foreach($result as $read) {
		$read_id = $read['id'];
		// Remove all lanes for that read
		$db_connect->executeStmt("DELETE FROM runqc_lane WHERE runqc_read_id = $read_id ;");
		// Remove entry for that read
		$db_connect->executeStmt("DELETE FROM runqc_read WHERE id = $read_id ;");
	}
}
	
// Insert statistics for reads and lanes into db
if(count($result) == 0 || $force) {
	foreach($read_metrics as $read_num => $read_metric) {
		// Insert read
		$hash = $db_connect->prepare("INSERT INTO runqc_read (read_num, cycles, is_index, q30_perc, error_rate, sequencing_run_id) VALUES (:read_num, :cycles, :is_index, :q30_perc, :error_rate, :seq_run_id);");
		$db_connect->bind($hash, "read_num", $read_num);
		$db_connect->bind($hash, "cycles", $read_metric['num_cycles']);
		$db_connect->bind($hash, "is_index", $read_metric['indexed']);
		$db_connect->bind($hash, "q30_perc", $read_metric['q30']*100);
		$db_connect->bind($hash, "error_rate", $read_metric['error_rate']);
		$db_connect->bind($hash, "seq_run_id", $seq_run_id);
		$db_connect->execute($hash, true);
		
		// ID of the inserted read
		$read_id = $db_connect->lastInsertId();
		
		// Insert lane metrics for that read
		foreach($read_metric['lane_metrics'] as $lane_num => $lane_metric) {
			// Insert this lane/tile/cycle
			$hash = $db_connect->prepare("INSERT INTO runqc_lane (lane_num, cluster_density, cluster_density_pf, yield, error_rate, q30_perc, runqc_read_id) VALUES (:lane_num, :cl_dens, :cl_dens_pf, :yield, :error_rate, :q30_perc, :runqc_read_id);");
			$db_connect->bind($hash, "lane_num", $lane_num);
			$db_connect->bind($hash, "cl_dens", $lane_metric['cluster_dens']);
			$db_connect->bind($hash, "cl_dens_pf", $lane_metric['cluster_dens_pf']);
			$db_connect->bind($hash, "yield", $lane_metric['yield']);
			if (array_key_exists('error_rate',$lane_metric))
			{
				$db_connect->bind($hash, "error_rate", $lane_metric['error_rate']);
			}
			else
			{
				$db_connect->bind($hash, "error_rate", null);
			}
			$db_connect->bind($hash, "q30_perc", $lane_metric['q30']*100);
			$db_connect->bind($hash, "runqc_read_id", $read_id);
			$db_connect->execute($hash, true);
		}
	}
}
?>