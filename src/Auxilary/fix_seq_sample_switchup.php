<?php
/** 
	@page fix_seq_sample_switchup
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("fix_seq_sample_switchup", "Gather switched data for documentation. Then correct the values in NGSD and move the necessary files.");
$parser->addInfile("samples", "Table which sample is wrong and which sample contains the correct data.", false);
$parser->addString("backup_folder", "folder to save files and NGSD content in the 'before' state.", false);
$parser->addFlag("commit", "Commits the changes into NGSD.");
$parser->addFlag("backup_overwrite", "remove data in backup_folder and copy the backup data again.");

extract($parser->parse($argv));


function sqlTableLine(&$db, $query)
{
	$res = $db->executeQuery($query);
	
	if (count($res) != 1)
	{
		var_dump($res);
		trigger_error("sqlTableLine can only be used with queries that return a single row as result! Query: ".$query, E_USER_ERROR);
	}
	
	$headers = array_keys($res[0]);
	$values = array_values($res[0]);
	
	return [$headers, $values];
}

function sqlTableLineSwitchup(&$db, $query, $id_before, $id_data, $switchup_id)
{
	list($headers_before, $values_before) = sqlTableLine($db, $query.$id_before);
	list($headers_data,   $values_data)   = sqlTableLine($db, $query.$id_data);
	
	$headers = ["#switchup_id"];
	$headers = array_merge($headers, array_map(fn($h) => 'sample1_' . $h, $headers_before));
	$headers = array_merge($headers, array_map(fn($h) => 'sample2_' . $h, $headers_data));
	
	$values = [$switchup_id];
	$values = array_merge($values, $values_before);
	$values = array_merge($values, $values_data);
	
	//clean values data: remove '\n', '\t' -> would break output TSV
	$values = array_map(fn($v) => str_replace("\n", "<newLine>", str_replace("\t", "<tab>", $v)), $values);
	
	return [$headers, $values];
}

function sqlTableBackup(&$db, $query, $key, $switched_samples, $out, $backup_overwrite)
{
	$headers_set = false;
	$lines = [];
	
	if (! $backup_overwrite && is_file($out))
	{
		echo "Backup file already exists! Skipping $out\n";
		return;
	}
	
	echo "creating backup: ".basename($out)."...";
	foreach($switched_samples as $samples)
	{
		$switchup_id = $samples["before"]["ps_name"]."<-".$samples["after"]["ps_name"];
		$id_before = $samples["before"][$key];
		$id_after = $samples["after"][$key];
		
		list($res_headers, $res_values) = sqlTableLineSwitchup($db, $query, $id_before, $id_after, $switchup_id);
		
		if (! $headers_set)
		{
			$lines[] = implode("\t", $res_headers);
			$headers_set = true;
		}
		$lines[] = implode("\t", $res_values);
	}
	
	file_put_contents($out, implode("\n", $lines));
	echo "finished!\n";
}


function copySampleFiles($switched_samples, $backup_folder_files, $backup_overwrite)
{
	global $parser;
	global $file_change_done;
	global $file_backup_done;
	
	if (is_file($file_backup_done) && ! $backup_overwrite)
	{
		echo "Status file backup done exists: $file_backup_done\nFile backup was already completed! Skipping use 'backup_overwrite' to restart backup.\n";
		return;
	}
	
	if (is_file($file_change_done))
	{
		trigger_error("File changes were started cannot update the backup files anymore!", E_USER_ERROR);
	}
	
	foreach($switched_samples as $samples)
	{
		$sample_before = $samples["before"]["ps_folder"];
		$sample_after   = $samples["after"]["ps_folder"];
		
		foreach([$sample_before, $sample_after] as $folder)
		{
			$dest = $backup_folder_files.basename($folder);
			
			echo "copying sample dir: rsync    -r -a $folder $dest\n";
			$parser->exec(           "rsync", "-r -a $folder $dest");
		}
	}
	
	if (!is_file($file_backup_done)) touch($file_backup_done);
}

function gatherData(&$db, $switched_samples, $backup_folder, $backup_folder_files, $backup_overwrite)
{	
	//backup_files:
	global $ps_table_backup;
	global $sample_table_backup;
	global $run_table_backup;
	
	#gather NGSD data:
	$ps_table_query = "SELECT * FROM processed_sample WHERE id=";
	sqlTableBackup($db, $ps_table_query,     "ps_id",  $switched_samples, $ps_table_backup,     $backup_overwrite);
	$sample_table_query = "SELECT * FROM sample WHERE id=";
	sqlTableBackup($db, $sample_table_query, "s_id",   $switched_samples, $sample_table_backup, $backup_overwrite);
	$run_table_query = "SELECT * FROM sequencing_run WHERE id=";
	sqlTableBackup($db, $run_table_query,    "run_id", $switched_samples, $run_table_backup,    $backup_overwrite);
	
	#create backup copies of data:
	copySampleFiles($switched_samples, $backup_folder_files, $backup_overwrite);
}

function prepareQueryChanges($switched_samples, $backup_folder, $backup_folder_files)
{
	//to fix the switchup read the data from the backup files -> as they stay the same
	//two sets of changes: NGSD and files
	
	global $ps_table_backup;
	
	//NGSD - "rename" sample1 by changing sample id
	//fields to correct:
	//	processed_sample table: sample_id, molarity, input, project_id, process_id (add ps_comment?)
	$query_data_main = [];
	$query_data_fix_process_id = [];
	
	$headers = [];
	$lines = file($ps_table_backup);
	$headers = explode("\t", $lines[0]);
	$headers[count($headers)-1] = trim($headers[count($headers)-1]);
	
	$idx_switchup_id     = array_search("#switchup_id", $headers);
	$idx_ps_id      = array_search("sample1_id", $headers);
	$idx_sample_id  = array_search("sample2_sample_id", $headers);
	$idx_molarity   = array_search("sample2_molarity", $headers);
	$idx_input      = array_search("sample2_processing_input", $headers);
	$idx_project_id = array_search("sample2_project_id", $headers);
	$idx_process_id = array_search("sample2_process_id", $headers);
	$idx_urgent = array_search("sample2_urgent", $headers);
	$idx_comment = array_search("sample2_comment", $headers);
	
	if (in_array(false, [$idx_ps_id, $idx_sample_id, $idx_molarity, $idx_input, $idx_project_id, $idx_process_id], true))
	{
		trigger_error("Couldn't find all necessary columns to create the NGSD UPDATE queryies: Missing columns in file: '$ps_table_backup'", E_USER_ERROR);
	}
	
	//A temporary process id has to be used for the switch as there is a unique index on "sample_id + process_id"
	$tmp_process_id = 99;
	$query_str_main = "UPDATE processed_sample SET sample_id=:sample_id, molarity=:molarity, processing_input=:processing_input, project_id=:project_id, process_id=:process_id, urgent=:urgent, comment=:comment WHERE id=:id";
	$query_str_fix_process_id = "UPDATE processed_sample SET process_id=:process_id WHERE id=:id";
	
	for($switchup_idx=1; $switchup_idx<count($lines); $switchup_idx++)
	{
		$line = $lines[$switchup_idx];
		
		$parts = explode("\t", $line);
		$parts[count($parts)-1] = trim($parts[count($parts)-1]);
		
		$switchup_id = $parts[$idx_switchup_id];
		$sample1 = explode("<-", $switchup_id)[0];
		$sample2 = explode("<-", $switchup_id)[1];
		
		$ps_id      = $parts[$idx_ps_id];
		$sample_id  = $parts[$idx_sample_id];
		$molarity   = $parts[$idx_molarity];
		$input      = $parts[$idx_input];
		$project_id = $parts[$idx_project_id];
		$process_id = $parts[$idx_process_id];
		$urgent = $parts[$idx_urgent];
		$comment = $parts[$idx_comment];
		
		//add comment to document the change:
		$comment = $comment." Was part of a sample switchup. This was previously {$sample1} and was corrected to {$sample2} on the ".date("Y-m-d H:i")." by script.";
		
		if (intval($process_id) == $tmp_process_id)
		{
			trigger_error("Temporary process id clashes with actual process id! Switchup: ".$parts[0], E_USER_ERROR);
		}
		
		$query_data_main[] = ["query_str"=>$query_str_main, "binds"=>["id"=>$ps_id, "sample_id"=>$sample_id, "molarity"=>$molarity, "processing_input"=>$input, "project_id"=>$project_id, "process_id"=>$tmp_process_id, "urgent"=>$urgent, "comment"=>$comment]];
		$query_data_fix_process_id[] = ["query_str"=>$query_str_fix_process_id, "binds"=>["id"=>$ps_id, "process_id"=>$process_id]];
	}
	
	$query_data = array_merge($query_data_main, $query_data_fix_process_id);
	
	return $query_data;
}

function prepareFileChanges($sample_switches, $backup_folder, $backup_folder_files)
{
	//files:
	//move necessary files for reanalysis from "before" sample into "after" sample folder
	//	remove current files
	//	copy and generate fastqs from .cram
	
	global $parser;
	
	$file_changes = [];
	foreach($sample_switches as $samples)
	{
		$before_ps_info = $samples["before"];
		$after_ps_info = $samples["after"];
		
		
		//remove old data:
		$cmd = "rm";
		$args = [];
		$args[] = "-r";
		$args[] = $after_ps_info["ps_folder"]."/*";
		
		if (!is_dir($after_ps_info["ps_folder"]))
		{
			trigger_error("Destination processed sample doesn't exist!", E_USER_ERROR);
		}
		
		$file_changes[] = $cmd." ".implode(" ", $args);
		
		
		//source the data from the backup folders
		$source_folder_name = basename($before_ps_info["ps_folder"]);
		$source_folder = $backup_folder_files.$source_folder_name."/";
		
		$dest_folder = $after_ps_info["ps_folder"]."/";
		
		//convert cram to fastqs in correct folder: 
		
		if (!is_dir($source_folder) || !is_dir($dest_folder))
		{
			trigger_error("One of processed sample folders doesn't exist!\nSource folder in backup: '$source_folder'\nDestination folder: '$dest_folder'", E_USER_ERROR);
		}
		
		$source_bam_or_cram = $source_folder.basename($before_ps_info["ps_bam"]);
		
		if (!is_file($source_bam_or_cram))
		{
			trigger_error("The source BAM or CRAM file doesn't exist: $source_bam_or_cram", E_USER_ERROR);
		}
		
		$dest_fastq1 = $dest_folder.$after_ps_info["ps_name"]."_switchup_fix_R1_001.fastq.gz";
		$dest_fastq2 = $dest_folder.$after_ps_info["ps_name"]."_switchup_fix_R2_001.fastq.gz";
		
		//COMMAND ONLY!!
		$cmd = $parser->execApptainer("ngs-bits", "BamToFastq", "-in $source_bam_or_cram -out1 $dest_fastq1 -out2 $dest_fastq2", [$source_bam_or_cram], [$dest_fastq1, $dest_fastq2], TRUE);
		$file_changes[] = $cmd;
	}
	
	return $file_changes;
}

function printFileChanges($file_changes)
{
	global $parser;
	global $file_change_doc;
	
	$count = 0;
	foreach($file_changes as $change)
	{
		echo $count." print: ".$change."\n";
		$count++;
	}
	
	file_put_contents($file_change_doc, implode("\n", $file_changes));
}

function executeFileChanges($file_changes)
{
	global $parser;
	global $file_change_done;
	
	if (!is_file($file_change_done))
	{
		touch($file_change_done);
	}
	
	$commands_done = file($file_change_done);
	$commands_done = array_map("trim", $commands_done);

	foreach($file_changes as $change)
	{
		if (in_array(trim($change), $commands_done))
		{
			echo "already done: $change\n";
			continue;
		}
		echo "File change being done: $change\n";
		exec2($change);
		file_put_contents($file_change_done, $change."\n", FILE_APPEND);
	}
}

function printNGSDChanges($query_data)
{
	global $parser;
	global $queries_doc;
	
	$count = 0;
	foreach($query_data as $query)
	{
		$query_str = $query["query_str"];
			
		foreach($query["binds"] as $name => $value)
		{
			$query_str = str_replace(":".$name, "'".$value."'", $query_str);
		}
		
		echo $count."print: ".$query_str."\n";
		$lines[] = $query_str;
		$count++;
	}
	
	file_put_contents($queries_doc, implode("\n", $lines));
	
}

function executeNGSDChanges(&$db, $query_data)
{
	global $queries_COMMITTED;
	
	if (is_file($queries_COMMITTED))
	{
		echo "Queries were already committed. Skipping!\n";
		return;
	}
	
	$db->beginTransaction();
	
	try
	{
		foreach($query_data as $query)
		{
			$hash = $db->prepare($query["query_str"]);
			
			foreach($query["binds"] as $name => $value)
			{
				$db->bind($hash, $name, $value);
			}
		
			$db->execute($hash);
		}
	} 
	catch (Exception $e)
	{
		$db->rollBack();
		throw $e;
	}
	
	$db->endTransaction();
	
	touch($queries_COMMITTED);
}

function printChanges($query_data, $file_changes)
{
	printNGSDChanges($query_data);
	printFileChanges($file_changes);
}


$db = DB::getInstance("NGSD");
$switched_samples = [];
$idx_before = -1;
$idx_data = -1;

$msg_sample_file_error = "This scripts expects the samples input table to be a TSV file with exactly two columns: 'current_processed_sample_name' and 'processed_sample_containing_correct_data'. The header needs to be marked with '#'";

//create log file in backup folder if none is provided
$backup_folder_files = $backup_folder."/sample_data/";

#backup_files
$ps_table_backup     = $backup_folder."/processed_sample_table_data.tsv";
$sample_table_backup = $backup_folder."/sample_table_data.tsv";
$run_table_backup    = $backup_folder."/run_table_data.tsv";
#documentation files
$queries_doc         = $backup_folder."/queries.tsv";
$file_change_doc     = $backup_folder."/file_change_commands.tsv";
#status files
$file_backup_done     = $backup_folder."/file_backup_DONE.txt";
$queries_COMMITTED         = $backup_folder."/queries_COMITTED.txt";
$file_change_done     = $backup_folder."/file_change_commands_DONE.tsv";
$commit_file = $backup_folder."/ALL_COMMITTED.txt";

if (!is_dir($backup_folder)) mkdir($backup_folder);
if (!is_dir($backup_folder_files)) mkdir($backup_folder_files);

if (is_file($commit_file))
{
	trigger_error("Cannot modify the backup data in '$backup_folder'. The changes were already committed.\n The folder and files are now for documentation purposes only!", E_USER_ERROR);
}


if ($parser->getLogFile()=="") $parser->setLogFile($backup_folder."/fix_seq_sample_switchup.log");

//log server, user, etc.
$parser->logServerEnvronment();

$samples_before = [];
$samples_data = [];

echo "parsing input....";
foreach(file($samples) as $line)
{
	if (trim($line) == "") continue;
	
	$line = trim($line);
	
	if ($line[0] == "#")
	{
		#check header:
		$parts = explode("\t", $line);
		
		$idx_before = array_search("#current_processed_sample_name", $parts);
		$idx_data = array_search("processed_sample_it_should_be", $parts);
		
		if (count($parts) != 2 || $parts[0] != "#current_processed_sample_name" || $parts[1] != "processed_sample_it_should_be")
		{
			var_dump($parts);
			trigger_error("1".$msg_sample_file_error, E_USER_ERROR);
		}
		continue;
	}
	
	if ($idx_before == -1 || $idx_data == -1) trigger_error("2".$msg_sample_file_error, E_USER_ERROR);
	
	
	$parts = explode("\t", $line);
	$sample_before = $parts[$idx_before];
	$sample_data = $parts[$idx_data];
	
	if (array_key_exists($sample_before, $samples_before))
	{
		trigger_error("Sample '$sample_before' exists at  least twice as 'before' sample. This is not allowed!", E_USER_ERROR);
	}
	$samples_before[$sample_before] = true;
	
	if (array_key_exists($sample_data, $samples_data))
	{
		trigger_error("Sample '$sample_data' exists at least twice as 'after' sample. This is not allowed!", E_USER_ERROR);
	}
	$samples_data[$sample_data] = true;
	
	$switched_samples[] = ["before" => get_processed_sample_info($db, $parts[$idx_before]), "after" => get_processed_sample_info($db, $parts[$idx_data])];
}

echo "DONE\n";

//all samples need to be before and after in one case:
foreach($samples_before as $s)
{
	if (! array_key_exists($s, $samples_data))
	{
		trigger_error("The sample '$s' is used as 'current_processed_sample_name' but not as 'processed_sample_it_should_be'. This script expects the switchup to be balanced.", E_USER_ERROR);
	}
}

$p_systems = [];
foreach($switched_samples as $samples)
{
	foreach($samples as $s)
	{
		if (!array_key_exists($s["sys_id"], $p_systems))
		{
			$ps_systems[$s["sys_id"]] = true;
		}
	}
}

if (count($p_systems) > 1)
{
	trigger_error("Not all samples have the same processing system. This is not supported!", E_USER_ERROR);
}

if (!$commit)
{
	gatherData($db, $switched_samples, $backup_folder, $backup_folder_files, $backup_overwrite);
	$query_data = prepareQueryChanges($switched_samples, $backup_folder, $backup_folder_files);
	$file_changes = prepareFileChanges($switched_samples, $backup_folder, $backup_folder_files);
	printChanges($query_data, $file_changes);
}
else
{
	$query_data = prepareQueryChanges($switched_samples, $backup_folder, $backup_folder_files);
	$file_changes = prepareFileChanges($switched_samples, $backup_folder, $backup_folder_files);
	executeNGSDChanges($db, $query_data);
	executeFileChanges($file_changes);
	
	//mark this switchup complete -> any further calls to this data will be rejected
	file_put_contents($commit_file, "All commands executed: ".date("Y-m-d H:i:s"));
	#TODO add chmod to make all files read-only
}






?>
