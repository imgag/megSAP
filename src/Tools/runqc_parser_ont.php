<?php
/**
 * @page runqc_parser_ont
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("runqc_parser_ont", "Parses and imports ONT run QC metrics into the NGSD.");
$parser->addString("name", "Name of the run.", false);
$parser->addString("run_dir", "Run directory.", false);
$parser->addFlag("force", "Overwrites already existing DB entries instead of throwing an error.");
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
$parser->addFlag("ignore_flowcell_id_check", "Ignores Errors due to mismatching flowcell IDs between NGSD and folder name (e.g. if flowcell was replaced during sequencing)");
$parser->addFlag("no_db", "No NGSD import, only report run QC to STDOUT.");
extract($parser->parse($argv));

//init
$db_connect = DB::getInstance($db);
$run_dir = realpath($run_dir);

//check that flowcell ID is set in the NGSD
$result = $db_connect->executeQuery("SELECT id, fcid FROM sequencing_run WHERE name = '$name'");
if(count($result) != 1)	trigger_error("Sequencing run with the name $name not found in DB.", E_USER_ERROR);
$seq_run_id = $result[0]['id'];
$flowcell_id = $result[0]['fcid'];

//check run folder for run name suffix
if (!ends_with($run_dir, substr($name, 1))) trigger_error("Run folder name ('".basename2($run_dir)."') doesn't match run name ({$name}!", E_USER_ERROR);

// check fcid and raw data folder
if ($ignore_flowcell_id_check) $fc_dirs = glob("{$run_dir}/2*_*_*");
else $fc_dirs = glob("{$run_dir}/2*_{$flowcell_id}_*");
if (count($fc_dirs) == 0) trigger_error("Run raw data folder missing or wrong FlowCell ID in NGSD!", E_USER_ERROR);
if (count($fc_dirs) > 1) trigger_error("Multiple raw data folders found!\n".implode(", ", $fc_dirs), E_USER_WARNING);

$run_qc_values = array();

foreach ($fc_dirs as $fc_dir) 
{
	
	// check report file
	if ($ignore_flowcell_id_check) $report_files = glob("{$fc_dir}/report_*_*.json");
	else $report_files = glob("{$fc_dir}/report_{$flowcell_id}_*.json");
	if (count($report_files) == 0) trigger_error("No report JSON file found!", E_USER_ERROR);
	if (count($report_files) > 1) trigger_error("Multiple report JSON files found!\n".implode(", ", $report_files) , E_USER_ERROR);

	//parse report file
	$report_content = (array) json_decode(file_get_contents($report_files[0], true));

	//extract run info
	$run_info = $report_content["protocol_run_info"];
	$protocol_id = $run_info->protocol_id;
	$software_args = (array) $run_info->args;
	$device_firmware_versions = array();
	foreach ((array) $run_info->device->firmware_version as $component) 
	{
		$device_firmware_versions[] = $component->component.":".$component->version;
	}
	$minknow_version = $run_info->software_versions->minknow->full;


	//use the last acquisition to extract runQC
	$final_acquisition = end($report_content["acquisitions"]);

	//check for correct acquisition
	$stop_reason = $final_acquisition->acquisition_run_info->stop_reason;
	if ($stop_reason != "STOPPED_PROTOCOL_ENDED") trigger_error("Invalid stop reason '{$stop_reason}'!", E_USER_WARNING);


	$read_num = $final_acquisition->acquisition_run_info->yield_summary->read_count;
	$yield = $final_acquisition->acquisition_run_info->yield_summary->basecalled_pass_bases + $final_acquisition->acquisition_run_info->yield_summary->basecalled_fail_bases;
	$basecalled_pass_bases = $final_acquisition->acquisition_run_info->yield_summary->basecalled_pass_bases;
	$passing_filter_perc = ($basecalled_pass_bases / $yield) * 100.0;

	// check for non-basecalled bases
	$fraction_skipped = 0.0;
	if (isset($final_acquisition->acquisition_run_info->yield_summary->fraction_skipped)) 
	{
		$fraction_skipped = $final_acquisition->acquisition_run_info->yield_summary->fraction_skipped;
	}
	if ($fraction_skipped > 1e-5) trigger_error(($fraction_skipped*100)."% of bases are not basecalled! Consider manual basecalling!", E_USER_WARNING);


	// get q30/q20 values
	$qscore_histograms_bases = ((array) $final_acquisition->qscore_histograms)[1];
	if ($qscore_histograms_bases->bucket_value_type != "QScore_BasecalledBases") trigger_error("QScore histogram not found!", E_USER_ERROR);

	//get qscores of buckets
	$qscores = array();
	$ranges = (array) $qscore_histograms_bases->bucket_ranges;
	foreach($ranges as $range)
	{
		if (isset($range->end)) $qscores[] = $range->end;
		else $qscores[] = PHP_INT_MAX; //last bucket has no end	 
	}

	//get combined counts of buckets
	$combined_counts = array_fill(0, count($qscores), 0);
	$histogram_data = (array) $qscore_histograms_bases->histogram_data;
	foreach($histogram_data as $data_object)
	{
		$counts = (array) $data_object->bucket_values;
		if (count($qscores) != count($counts)) trigger_error("Count size doesn't match number of buckets! (".count($counts)."vs".count($qscores).")", E_USER_ERROR);
		for ($fc_dir_name=0; $fc_dir_name < count($qscores); $fc_dir_name++) 
		{ 
			$combined_counts[$fc_dir_name] += $counts[$fc_dir_name];
		}
	}
	//calculate q30 value
	$all_bases = 0;
	$q30_bases = 0;
	$q20_bases = 0;
	for ($fc_dir_name=0; $fc_dir_name < count($qscores); $fc_dir_name++) 
	{ 
		$all_bases += $combined_counts[$fc_dir_name];
		if ($qscores[$fc_dir_name] >= 30) $q30_bases += $combined_counts[$fc_dir_name];
		if ($qscores[$fc_dir_name] >= 20) $q20_bases += $combined_counts[$fc_dir_name];
	}
	$q30_perc = ($q30_bases/$all_bases) * 100.0;
	$q20_perc = ($q20_bases/$all_bases) * 100.0;

	//get N50 value
	$n50 = -1;
	$read_length_histograms = (array) $final_acquisition->read_length_histogram;
	foreach ($read_length_histograms as $histogram)
	{
		if (isset($histogram->read_length_type) && ($histogram->read_length_type == "BasecalledBases") 
			&& isset($histogram->bucket_value_type) && ($histogram->bucket_value_type == "ReadLengths"))
		{
			// var_dump($histogram->plot->histogram_data);
			$n50 = $histogram->plot->histogram_data[0]->n50;
			break;
		}
	}

	if ($n50 < 0) trigger_error("N50 value not found!", E_USER_ERROR);

	// collapse arrays to single string
	$software_args = implode(" ", $software_args);
	$device_firmware_versions = implode(";",$device_firmware_versions);


	//store per-folder QCs:
	$folder_qc = array();
	$folder_qc["run"] = $name;
	$folder_qc["fc_id"] = $flowcell_id;
	$folder_qc["read_num"] = $read_num;
	$folder_qc["yield"] = $yield;
	$folder_qc["passing_filter_perc"] = $passing_filter_perc; 
	$folder_qc["fraction_skipped"] = $fraction_skipped;
	$folder_qc["q30_perc"] = $q30_perc;
	$folder_qc["q20_perc"] = $q20_perc;
	$folder_qc["N50"] = $n50;
	$folder_qc["protocol_id"] = $protocol_id; 
	$folder_qc["software_args"] = $software_args;
	$folder_qc["device_firmware_versions"] = $device_firmware_versions;
	$folder_qc["minknow_version"] = $minknow_version;

	$run_qc_values[basename($fc_dir)] = $folder_qc;
}

//get total/combined values
$flowcell_ids = array();
$total_read_num = 0;
$total_yield = 0;
$passing_filter_percentages = array();
$fraction_skipped_list = array();
$q30_percentages = array();
$q20_percentages = array();
$n50s = array();
$protocol_ids = array();
$software_arguments = array();
$device_firmware_version_list = array();
$minknow_versions = array();

foreach ($run_qc_values as $fc_dir_name => $folder_qc) 
{
	$flowcell_ids[] = $folder_qc["fc_id"];
	$total_read_num += $folder_qc["read_num"];
	$total_yield += $folder_qc["yield"];
	$passing_filter_percentages[] = $folder_qc["passing_filter_perc"];
	$fraction_skipped_list[] = $folder_qc["fraction_skipped"];
	$q30_percentages[] = $folder_qc["q30_perc"];
	$q20_percentages[] = $folder_qc["q20_perc"];
	$n50s[] = $folder_qc["N50"];
	$protocol_ids[] = $folder_qc["protocol_id"];
	$software_arguments[] = $folder_qc["software_args"];
	$device_firmware_version_list[] = $folder_qc["device_firmware_versions"];
	$minknow_versions[] = $folder_qc["minknow_version"];
}

//check argument fields
if (count(array_unique($flowcell_ids)) > 1) trigger_error("Flocell IDs of sub folders differ!\n".implode(", ", $flowcell_ids), E_USER_ERROR);
$flowcell_id = array_unique($flowcell_ids)[0];
if (count(array_unique($protocol_ids)) > 1) trigger_error("Protocol IDs of sub folders differ!\n".implode(", ", $protocol_ids), E_USER_ERROR);
$protocol_id = array_unique($protocol_ids)[0];
//workaround for restarted runs with FAST
if (count(array_unique($software_arguments)) > 1) trigger_error("Software arguments of sub folders differ!\n".implode(", ", $software_arguments), E_USER_WARNING); 
$software_args = array_unique($software_arguments)[0];
if (count(array_unique($device_firmware_version_list)) > 1) trigger_error("Device firmware versions of sub folders differ!\n".implode(", ", $device_firmware_version_list), E_USER_ERROR);
$device_firmware_versions = array_unique($device_firmware_version_list)[0];
if (count(array_unique($minknow_versions)) > 1) trigger_error("MinKNOW versions of sub folders differ!\n".implode(", ", $minknow_versions), E_USER_ERROR);
$minknow_version = array_unique($minknow_versions)[0];

//calculate weight factors for yield and read numbers
$wf_reads = array();
$wf_yield = array();
foreach ($run_qc_values as $fc_dir_name => $folder_qc) 
{
	$wf_reads[$fc_dir_name] = floatval($folder_qc["read_num"] / $total_read_num);
	$wf_yield[$fc_dir_name] = floatval($folder_qc["yield"] / $total_yield);
}

//calculate weighted QC values
$read_num = $total_read_num;
$yield = $total_yield;
$passing_filter_perc = 0.00;
$fraction_skipped = 0.00;
$q30_perc = 0.00;
$q20_perc = 0.00;
$n50 = 0;
foreach ($run_qc_values as $fc_dir_name => $folder_qc)
{ 
	$passing_filter_perc += $wf_yield[$fc_dir_name] * $run_qc_values[$fc_dir_name]["passing_filter_perc"];
	$fraction_skipped += $wf_yield[$fc_dir_name] * $run_qc_values[$fc_dir_name]["fraction_skipped"];
	$q30_perc += $wf_yield[$fc_dir_name] * $run_qc_values[$fc_dir_name]["q30_perc"];
	$q20_perc += $wf_yield[$fc_dir_name] * $run_qc_values[$fc_dir_name]["q20_perc"];
	$n50 += $wf_reads[$fc_dir_name] * $run_qc_values[$fc_dir_name]["N50"];
}

//print TSV to STDOUT
print "#run\tfc_id\tread_num\tyield\tpassing_filter_perc\tfraction_skipped\tq30_perc\tq20_perc\tN50\tprotocol_id\tsoftware_args\tdevice_firmware_versions\tminknow_version\tsub_folder\n";
foreach ($run_qc_values as $fc_dir_name => $folder_qc) 
{
	$parts = array($name);
	$parts[] = $folder_qc["fc_id"];
	$parts[] = $folder_qc["read_num"];
	$parts[] = $folder_qc["yield"];
	$parts[] = $folder_qc["passing_filter_perc"];
	$parts[] = $folder_qc["fraction_skipped"];
	$parts[] = $folder_qc["q30_perc"];
	$parts[] = $folder_qc["q20_perc"];
	$parts[] = $folder_qc["N50"];
	$parts[] = $folder_qc["protocol_id"];
	$parts[] = $folder_qc["software_args"];
	$parts[] = $folder_qc["device_firmware_versions"];
	$parts[] = $folder_qc["minknow_version"];
	$parts[] = $fc_dir_name;
	print substr(implode("\t", $parts), 1)."\n";
}
print substr($name, 1)."\t{$flowcell_id}\t{$read_num}\t{$yield}\t{$passing_filter_perc}\t{$fraction_skipped}\t{$q30_perc}\t{$q20_perc}\t{$n50}\t{$protocol_id}\t{$software_args}\t{$device_firmware_versions}\t{$minknow_version}\tcombined\n";


if (!$no_db)
{
	//check if runqc data was already imported for this run
	$result = $db_connect->executeQuery("SELECT id FROM runqc_ont WHERE sequencing_run_id = {$seq_run_id}");
	if(count($result)!=0)
	{
		if ($force) //remove data
		{
			foreach($result as $row)
			{
				$db_connect->executeStmt("DELETE FROM runqc_ont WHERE id = ".$row['id']);
			}
		}
		else //abort
		{
			trigger_error("QC data for this run and read was already imported. Use the flag '-force' to overwrite it.", E_USER_ERROR);
		}
	}

	//insert read
	$hash = $db_connect->prepare("INSERT INTO runqc_ont (sequencing_run_id, read_num, yield, passing_filter_perc, fraction_skipped, q30_perc, q20_perc, n50, protocol_id, software_args, device_firmware_versions, minknow_version) "
									."VALUES (:sequencing_run_id, :read_num, :yield, :passing_filter_perc, :fraction_skipped, :q30_perc, :q20_perc, :n50, :protocol_id, :software_args, :device_firmware_versions, :minknow_version);");
	$db_connect->bind($hash, "sequencing_run_id", $seq_run_id);
	$db_connect->bind($hash, "read_num", $read_num);
	$db_connect->bind($hash, "yield", $yield);
	$db_connect->bind($hash, "passing_filter_perc", $passing_filter_perc);
	$db_connect->bind($hash, "fraction_skipped", $fraction_skipped);
	$db_connect->bind($hash, "q30_perc", $q30_perc);
	$db_connect->bind($hash, "q20_perc", $q20_perc);
	$db_connect->bind($hash, "n50", $n50);
	$db_connect->bind($hash, "protocol_id", $protocol_id);
	$db_connect->bind($hash, "software_args", $software_args);
	$db_connect->bind($hash, "device_firmware_versions", $device_firmware_versions);
	$db_connect->bind($hash, "minknow_version", $minknow_version);
	$db_connect->execute($hash, true);
}

?>