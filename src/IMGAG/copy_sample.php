<?php 
/** 
	@page copy_sample
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");


error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("copy_sample", "Creates a Makefile to copy de-multiplexed sample data to projects and queue data analysis.");
//optional
$parser->addString("samplesheet",  "Input samplesheet that was used to create the data folder. (Has to be located in the run folder.)", true, "SampleSheet_bcl2fastq.csv");
$parser->addString("folder",  "Input data folder.", true, "Unaligned");
$parser->addString("runinfo" ,"Illumina RunInfo.xml file. Necessary for checking metadata in NGSD",true, "RunInfo.xml");
$parser->addOutfile("out",  "Output Makefile. Default: 'Makefile'.", true);
$parser->addFlag("high_priority", "Assign high priority to all queued samples.");
$parser->addFlag("overwrite", "Do not prompt before overwriting FASTQ files.");
$parser->addFlag("no_rename_r3", "Do not rename R2/R3 FASTQ files to index/R2.");

$parser->addFlag("no_genlab", "Do not query GENLAB for metadata.");
$parser->addFlag("ignore_analysis", "Use FASTQ files from BCLConvert.");
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
$parser->addInt("threads_ora", "Number of threads used to decompress ORA files during the copy.", true, 15);
$parser->addFlag("test", "Run in test mode, e.g. set the pipeline path and project folder to a fixed value.");

$parser->addFlag("use_dragen", "Use DRAGEN to analyze the data.");
$parser->addFlag("no_dragen", "Do not use DRAGEN to analyze the data.");

$parser->addFlag("merge_runs", "Automatically merge this run if a previous run was detected.");
$parser->addFlag("skip_run_merging", "Do not merge this run with a previous run.");

$parser->addFlag("manual_demux", "Ignore NovaSeqX Analysis results and use maunal demux.");
$parser->addFlag("no_queuing", "Do not include queuing commands in the default target of the Makefile.");

$parser->addString("ignore_samples" ,"Comma-separated list of samples which should not be copied/analyzed",true, "");

extract($parser->parse($argv));

//check if GenLab is available
if (get_path("genlab_host", false)=="" || get_path("genlab_name", false)=="" || get_path("genlab_user", false)=="" || get_path("genlab_pass", false)=="")
{
	trigger_error("GenLab SQL settings not present: skipping GenLab import!", E_USER_NOTICE);
	$no_genlab = true;
}

//extract samples names and sequencer type from sample sheet
function extract_sample_data(&$db_conn, $filename)
{
	$file = file($filename);
	//remove everything above [DATA] + [DATA]:
	for ($i=0; $i < count($file); $i++) 
	{ 
		if(starts_with($file[$i], "[Data]"))
		{
			$file = array_slice($file, $i+1); //section before [DATA], [DATA]
			break;
		}
	}
	array_shift($file);//skip header

	//extract sample data
	$sample_data = array();
	foreach($file as $line)
	{
		if (trim($line)=="") continue;
		
		list(, , $sample) = explode(",", $line);
		$sample = trim($sample);
		
		$sample_data[$sample] = get_processed_sample_info($db_conn, $sample);
	}
	return $sample_data;
}

function create_mail_command(&$db_conn, $project_name, $samples)
{	
	$result = $db_conn->executeQuery("SELECT u.name, u.email FROM user u, project p WHERE p.internal_coordinator_id=u.id AND p.name=:name", array("name"=>$project_name));
	$coordinator = $result[0]['name'];
	$email = $result[0]['email'];
	
	$mail_text = array();
	$mail_text[] = "Hallo $coordinator,";
	$mail_text[] = "";
	$mail_text[] = "die FASTQ-Dateien der folgenden Proben des Projekts $project_name liegen vor:";
	$mail_text[] = "";
	foreach(array_unique($samples) as $sample)
	{
		$mail_text[] = $sample;
	}
	$mail_text[] = "";
	$mail_text[] = "Viele Gruesse";
	$mail_text[] = "";
	$mail_text[] = "  die Bioinformatik";
	return "php -r 'mail(\"$email\",\"Neue Daten fuer $project_name\", \"".implode("\\n",$mail_text)."\",\"Reply-To: no-reply@med.uni-tuebingen.de\");'";
}

//get file size in MB
function filesize_mb($file)
{
	return number_format(filesize($file)/1024/1024, 3, ".", "");
}

//Checks whether there is at least one sample in "sample_sheet" that has the same number of lanes as specified Illumina RunInfo.xml. Returns false if RunInfo XML-file or samplesheet could not be parsed.
function check_number_of_lanes($run_info_xml_file, $sample_sheet)
{
	//RunInfo.xml contains attribute "LaneCount"
	if(!file_exists($run_info_xml_file))
	{
		return false;
	}
	$xml = simplexml_load_file($run_info_xml_file);
	if ($xml===false) trigger_error("Could not load XML file: $run_info_xml_file", E_USER_ERROR);
	if(empty($xml->Run->FlowcellLayout->attributes()->LaneCount))
	{
		return false;
	}
	$lane_count = (int) $xml->Run->FlowcellLayout->attributes()->LaneCount;
	
	
	$sample_sheet_data = Matrix::fromTSV($sample_sheet);
	$i_lane = $sample_sheet_data->getColumnIndex("Lane",false,false);
	if($i_lane == -1)
	if(empty($xml->Run->FlowcellLayout->attributes()->LaneCount))
	{
		return false;
	}
	
	//Check whether LaneCount from RunInfo.xml occurs at least one time in samplesheet
	
	//generate a cleaned sample sheet without adapter and sections
	$file = file($sample_sheet);
	for ($i=0; $i < count($file); $i++) 
	{ 
		if(starts_with($file[$i], "[Data]"))
		{
			$file = array_slice($file, $i+1); //section before [DATA], [DATA]
			break;
		}
	}
	$cleaned_sample_sheet = temp_file("sample_sheet.csv");
	file_put_contents($cleaned_sample_sheet, implode("\n", $file));
	$sample_sheet_content = Matrix::fromTSV($cleaned_sample_sheet, ",");
	$i_lane = array_search("Lane",  $sample_sheet_content->getRow(0)); // first row contains headers
	
	if($i_lane == -1) return false;
	
	$lane_count_found_in_samplesheet = false;
	for($i=1;$i<$sample_sheet_content->rows();++$i)
	{
		$lane = $sample_sheet_content->get($i,$i_lane);
		if($lane == $lane_count)
		{
			$lane_count_found_in_samplesheet = true;
			break;
		}
	}
	
	if(!$lane_count_found_in_samplesheet)
	{
		$message = array("",
		"*** !!!!!!!!!!!!!!! ***",
		"*** !!! WARNING !!! ***",
		"*** !!!!!!!!!!!!!!! ***",
		"*** Could not find any sample that uses $lane_count lanes as specified in {$run_info_xml_file}. ***",
		"*** CHECK NUMBER OF LANES IN DATABASE AND {$sample_sheet}! ***",
		""
		);
		trigger_error(implode("\n",$message),E_USER_WARNING);
	}
	return true;
}

//Checks whether there is former run containing more than 50% same sample names of current run
function check_former_run(&$db_conn, $run_name)
{
	$res = $db_conn->getValues("SELECT id FROM sequencing_run WHERE name='{$run_name}'");
	if(count($res) != 1)
	{
		trigger_error("Could not find info for run '{$run_name}'.\nNo check for former run to be merged...", E_USER_WARNING);
		return "";
	}
	$current_id = $res[0];
	
	$current_samples = get_processed_samples_from_run($db_conn, $run_name);
	for($i=0; $i<count($current_samples); ++$i)
	{
		$current_samples[$i] = explode("_",$current_samples[$i])[0];
	}
	
	//Check former 50 runs
	$res = $db_conn->executeQuery($query = "SELECT id,name FROM sequencing_run WHERE quality != 'bad' AND id < $current_id AND id > $current_id - 50");
	foreach($res as $res2)
	{
		$other_ps = get_processed_samples_from_run($db_conn, $res2["name"]);
		
		for($i=0; $i<count($other_ps); ++$i)
		{
			$other_ps[$i] = explode("_", $other_ps[$i])[0];
		}
		
		
		if(count(array_intersect($current_samples, $other_ps)) > count($current_samples)/2 )
		{
			return $res2["name"];
		}

	}
	return "";
}

//Returns father/mother of the sample
function get_trio_parents($ps)
{
	global $db_conn;
	
	$info = get_processed_sample_info($db_conn, $ps, false);
	if (is_null($info)) return NULL;
	$s_id = $info['s_id'];
	$sys_id = $info['sys_id'];
	
	//check if father/mother are defined
	$s_id_mother = $db_conn->getValue("SELECT s.id FROM sample_relations sr, sample s WHERE s.id=sr.sample1_id AND sr.relation='parent-child' AND sample2_id='{$s_id}' AND s.gender='female'", "");
	$s_id_father = $db_conn->getValue("SELECT s.id FROM sample_relations sr, sample s WHERE s.id=sr.sample1_id AND sr.relation='parent-child' AND sample2_id='{$s_id}' AND s.gender='male'", "");
	if ($s_id_mother=="" || $s_id_father=="") return NULL;
	
	//select latest processing of father/mother with same system
	$ps_id_mother = $db_conn->getValue("SELECT id FROM processed_sample WHERE sample_id={$s_id_mother} AND processing_system_id={$sys_id} ORDER BY id DESC LIMIT 1", "");
	$ps_id_father = $db_conn->getValue("SELECT id FROM processed_sample WHERE sample_id={$s_id_father} AND processing_system_id={$sys_id} ORDER BY id DESC LIMIT 1", "");
	if ($ps_id_mother=="" || $ps_id_father=="") return NULL;
	
	return array(processed_sample_name($db_conn, $ps_id_father), processed_sample_name($db_conn, $ps_id_mother));
}

//extract samples names and sequencer type from sample sheet
function get_sample_data_from_db(&$db_conn, $run_name)
{
	$samples = get_processed_samples_from_run($db_conn, $run_name);
	$sample_data = array();
	foreach ($samples as $sample)
	{
		$sample_data[$sample] = get_processed_sample_info($db_conn, $sample);
	}
	return $sample_data;
}

//get md5 checksum (NovaSeqX analysis)
function get_md5_line($file)
{
	$md5_file = $file.".md5sum";
	if (!file_exists($md5_file)) trigger_error("No MD5 file for {$file} found!", E_USER_ERROR);
		
	$md5 = trim(file($md5_file)[0]);

	return $md5." ".$file;
}

//init
if (!isset($out)) $out = "Makefile";
$db_conn = DB::getInstance($db);
if($test) $repo_folder = "/mnt/storage2/megSAP/pipeline"; //fixed absolute path to make the tests work for all users
else $repo_folder = repository_basedir(); //use repositories tool in production

//fallback for NovaSeq X default SampleSheet name
if (($samplesheet == "SampleSheet_bcl2fastq.csv") && !file_exists($samplesheet)) $samplesheet = "SampleSheet.csv";
//check SampleSheet
if(!file_exists($samplesheet)) trigger_error("ERROR: Provided SampleSheet '{$samplesheet}' not found!", E_USER_ERROR);

//get file names
$run_folder = dirname(realpath($samplesheet));

//get ignored samples
$ignored_samples = explode(",", $ignore_samples);

//get run id and check FlowCell id
$run_folder_parts = explode("_", basename($run_folder));
$run_name = $run_folder_parts[count($run_folder_parts) - 1];
if((strlen($run_name) != 5) || ((int) $run_name > 9999) || ((int)$run_name < 0)) trigger_error("ERROR: invalid run folder suffix '{$run_name}' provided!", E_USER_ERROR);
//add '#' prefix
$run_name = "#".$run_name;

//get flowcell id
if(!file_exists($runinfo)) trigger_error("ERROR: Run info file '{$runinfo}' not found!", E_USER_ERROR);
$xml = simplexml_load_file($runinfo);
if ($xml===false) trigger_error("Could not load XML file: $runinfo", E_USER_ERROR);
if(empty($xml->Run->Flowcell)) trigger_error("ERROR: Run info file doesn't contain flowcell id!", E_USER_ERROR);
$flowcell_id = $xml->Run->Flowcell;

//check flowcell id
$flowcell_id_ngsd = $db_conn->getValue("SELECT fcid FROM sequencing_run WHERE `name`='{$run_name}'");
if($flowcell_id != $flowcell_id_ngsd) trigger_error("ERROR: FlowCell id from run info doesn't match FlowCell id in NGSD! (Run info: {$flowcell_id} <-> NGSD: {$flowcell_id_ngsd})", E_USER_ERROR);

//get instrument type
$run_parameters_xml = $run_folder."/RunParameters.xml";
if(!file_exists($run_parameters_xml)) trigger_error("ERROR: Required RunParameters.xml file is missing in the run folder!", E_USER_ERROR);
$is_novaseq_x = is_novaseq_x_run($run_parameters_xml) && !$manual_demux;

//change default data folder for NovaSeqX 
if($is_novaseq_x && ($folder=="Unaligned")) 
{
	trigger_error("NovaSeqX run detected!", E_USER_NOTICE);
	$folder="Analysis";
}

//check that data folder exists
if (!file_exists($folder)) trigger_error("Data folder '$folder' does not exist!", E_USER_ERROR);

if(! $is_novaseq_x && !file_exists("Fq") && !check_number_of_lanes($runinfo,$samplesheet))
{
	trigger_error("***!!!WARNING!!!***\nCould not verify number of lanes used for Demultiplexing and actually used on Sequencer. Please check manually!", E_USER_WARNING);
}



//get analysis options

//DRAGEN:
if ($use_dragen && $no_dragen) trigger_error("Cannot use and not use DRAGEN!", E_USER_ERROR);
if (($use_dragen == false) && ($no_dragen == false) && (get_path("queues_dragen", false)!=""))
{
	// ask if DRAGEN should be used
	echo "Should the genome samples on this run mapped with DRAGEN mapping? (y/n)?\n";
	if(trim(fgets(STDIN)) == "y") $use_dragen = true;
}

//merge former run
if ($merge_runs && $skip_run_merging) trigger_error("Can either merge or skip merging runs!", E_USER_ERROR);
$merge_former_run = "";
if ($merge_runs) $merge_former_run = "y";
else if ($skip_run_merging) $merge_former_run = "n";


//get sample data
if($is_novaseq_x) $sample_data = get_sample_data_from_db($db_conn, $run_name);
else $sample_data = extract_sample_data($db_conn, $samplesheet);

//remove ignored samples
foreach ($ignored_samples as $sample) 
{
	unset($sample_data[$sample]);
}

//import data from Genlab
if (!$no_genlab)
{
	print "Importing information from GenLab...\n";

	foreach($sample_data as $sample => $data)
	{
		$args = [];
		$args[] = "-ps {$sample}";
		if ($db=="NGSD_TEST") $args[] = "-test";
		$parser->execApptainer("ngs-bits", "NGSDImportGenlab", implode(" ", $args));
	}

	//update sample data after importing sample relations
	if($is_novaseq_x) $sample_data = get_sample_data_from_db($db_conn, $run_name);
	else $sample_data = extract_sample_data($db_conn, $samplesheet);
}
//extract tumor-normal pair infos
$normal2tumor = array();
$tumor2normal = array();
foreach($sample_data as $sample => $sample_infos)
{
	$normal_name = $sample_infos['normal_name'];
	if ($normal_name!="")
	{
		$normal2tumor[$normal_name] = $sample;
		$tumor2normal[$sample] = $normal_name;
	}
	//check run name (the same for all samples)
	if($run_name != $sample_infos['run_name']) trigger_error("Sequencing run doesn't match sample info ('".$sample_infos['run_name']."')", E_USER_ERROR);
}

//Check for former run that contains more than 50% same samples and offer merging to user
$former_run_name = "";
$merge_files = array(); //contains merge commands, commands are inserted before queue

if($run_name != "")
{
	$former_run_name =  check_former_run($db_conn, $run_name);
	if($former_run_name != "")
	{
		// if no option is set by cli -> ask
		if ($merge_former_run == "")
		{
			echo "Former run '{$former_run_name}' detected. Merge (y/n)?\n";
			$merge_former_run = trim(fgets(STDIN));
		}
		
		if($merge_former_run == "y")
		{
			//get samples to merge
			$old_samples = array();
			foreach(get_processed_samples_from_run($db_conn, $former_run_name) as $ps)
			{
				$key = explode("_", $ps)[0];
				$old_samples[$key] = $ps;
			}
			$current_samples = array();
			
			foreach(get_processed_samples_from_run($db_conn, $run_name) as $ps)
			{
				$key = explode("_", $ps)[0];
				$current_samples[$key] = $ps;
			}
			
			$merge_files[] = "merge:";
			
			//match samples from old and new run
			foreach($current_samples as $current_key => $current_ps)
			{
				if(array_key_exists($current_key, $old_samples))
				{
					$merge_files[$current_ps] = "\tphp {$repo_folder}/src/Tools/merge_samples.php -ps ".$old_samples[$current_key]." -into $current_ps";
				}
			}
		}
		elseif($merge_former_run != "n")
		{
			trigger_error("Please choose whether files from both runs shall be merged.", E_USER_ERROR);
		}
	}
	else
	{
		//skip first run only on NovaSeq6000 
		if (!$is_novaseq_x && !$no_queuing)
		{
			// check if run is first genome run
			$processed_samples = get_processed_samples_from_run($db_conn, $run_name);
			$n_samples = count($processed_samples);
			$n_genomes = 0;
			foreach($processed_samples as $ps)
			{
				$processed_sample_info = get_processed_sample_info($db_conn,$ps);
				if ($processed_sample_info['sys_type'] == "WGS") $n_genomes++;
			}

			if (($n_genomes/$n_samples) >= 0.9)
			{
				// ask if queuing should be skipped
				echo "The current run seems to be the first of a genome run. Should the queuing of the samples be skipped? (y/n)?\n";
				if(trim(fgets(STDIN)) == "y")
				{
					$no_queuing = true;
				} 
			}
		}

	}
}
else
{
	trigger_error("Could not determine current run name. Skipping check for former run that could be merged into current run.", E_USER_WARNING);
}

//get analysis folder
$nsx_analysis_done = false;
if ($is_novaseq_x)
{
	$analyses = array_filter(glob($folder."/[0-9]"), 'is_dir');
	if(count($analyses) == 0) trigger_error("ERROR: No analysis found!", E_USER_ERROR);
	if(count($analyses) > 1) trigger_error("ERROR: Multiple analyses found!", E_USER_ERROR);
	$analysis_id = $analyses[0];

	//check for analysis
	$nsx_analysis_done = file_exists($analysis_id."/Data/DragenEnrichment/") || file_exists($analysis_id."/Data/DragenGermline/");
}



//parse input
$create_project_folder = array();
$target_to_copylines = array();
$target_to_queuelines = array();
$queue_trios = array();
$queue_somatic = array();
$project_to_fastqonly_samples = array();
$sample_to_newlocation = array();
$md5sum_buffer = array();

//get all processed samples which are scheduled for resequencing:
$samples_to_resequence = array();
$result = $db_conn->executeQuery("SELECT id, sample_id, processing_system_id FROM processed_sample ps WHERE scheduled_for_resequencing=1");
foreach($result as $row)
{
	$samples_to_resequence[$row["sample_id"]."_".$row["processing_system_id"]] = $row["id"];
}
$merge_notice = array("The following samples are scheduled for merging:");


foreach($sample_data as $sample => $sample_infos)
{
	$project_name = $sample_infos['project_name'];
	$project_type = $sample_infos['project_type'];
	$project_folder = $test ? "/test/project/folder/" : $sample_infos['project_folder'];
	$project_analysis = $sample_infos['project_analysis'];
	$sample_folder = $sample_infos['ps_folder'];
	$sample_is_tumor = $sample_infos['is_tumor'];
	$sys = $sample_infos['sys_name_short'];
	$sys_type = $sample_infos['sys_type'];
	$sys_target = $sample_infos['sys_target'];

	//check if sample should be merged
	$merge_sample = false;
	if (count($samples_to_resequence) > 0)
	{
		$key = $sample_infos["s_id"]."_".$sample_infos["sys_id"];
		if (isset($samples_to_resequence[$key]))
		{
			//merge sample
			$merge_sample = true;
			$source_ps_name = processed_sample_name($db_conn, $samples_to_resequence[$key]);
			$merge_files[$sample] = "\tphp {$repo_folder}/src/Tools/merge_samples.php -ps {$source_ps_name} -into {$sample}";
			$merge_notice[] = "\t{$source_ps_name} => {$sample}";
		}
	}

	//calculate current location of sample
	if ($is_novaseq_x)
	{
		$processing_system_ini = $parser->tempFile(".ini", $sys);
		store_system($db_conn, $sys, $processing_system_ini);
		$processing_system_info = parse_ini_file($processing_system_ini);
		$umi_type = $processing_system_info['umi_type'];
		$fastq_folder = $analysis_id."/Data/BCLConvert/fastq";

		$old_location = $analysis_id."/Data";
		if(($sys_type == "RNA") || ($sys_type == "Panel") || ($sys_type == "cfDNA (patient-specific)") || ($sys_type == "cfDNA") || ($sys_type == "WGS (shallow)"))
		{
			$old_location .= "/BCLConvert";
			if (file_exists($old_location."/ora_fastq"))
			{
				$fastq_folder = $old_location."/ora_fastq";
			}
			else
			{
				$fastq_folder = $old_location."/fastq";
			}
		}
		else if($sys_type == "WGS")
		{
			if (file_exists($old_location. "/DragenGermline")) $old_location .= "/DragenGermline";
			else $old_location .= "/BCLConvert"; // fallback to BCLConvert if no analysis is performed

			//check for ORA folder:
			if (file_exists($old_location."/ora_fastq"))
			{
				$fastq_folder = $old_location."/ora_fastq";
			}
			else
			{
				$fastq_folder = $old_location."/fastq";
			}

		}
		else if($sys_type == "WES")
		{
			if (file_exists($old_location. "/DragenEnrichment")) $old_location .= "/DragenEnrichment";
			else $old_location .= "/BCLConvert"; // fallback to BCLConvert if no analysis is performed
			if (file_exists($old_location."/ora_fastq"))
			{
				$fastq_folder = $old_location."/ora_fastq";
			}
			else
			{
				$fastq_folder = $old_location."/fastq";
			}
		}
		else trigger_error("ERROR: Invalid processing system type '{$sys_type}' for NovaSeq X analysis!", E_USER_ERROR);
	}
	else
	{
		$old_location = "{$folder}/{$project_name}";
	}
	
	//determine project 
	$sample_to_newlocation[$sample] = $sample_folder;
	
	//copy_target line if first sample of project in run, mkdir and chmod if first sample of project at all
	$tag = $sample;
	if (!is_dir($project_folder)) 
	{
		$create_project_folder[] = "\tmkdir -p {$project_folder}";
		$create_project_folder[] = "\tchmod 775 {$project_folder}";
	}
	
	//build copy line
	if($is_novaseq_x)
	{
		if(($sys_type == "RNA") || ($sys_type == "Panel") || ($sys_type == "cfDNA (patient-specific)") || ($sys_type == "cfDNA") || ($sys_type == "WGS (shallow)"))
		{
			//only copy FastQ files, analyzed by megSAP
			//get FastQs
			$fastq_files = glob($fastq_folder."/{$sample}*_L00[0-9]_R[123]_00[0-9].fastq.{gz,ora}", GLOB_BRACE);
			if($umi_type == "n/a" || $umi_type == "IDT-UDI-UMI" || $umi_type == "IDT-xGen-Prism" || $umi_type == "Twist")
			{
				//no index files: simply copy/convert FastQs
				
				//check count
				if(count($fastq_files) != count($sample_infos["ps_lanes"]) * 2) 
				{
					trigger_error("Number of FastQ files for sample {$sample} doesn't match number of lanes in run info! (expected: ".(count($sample_infos["ps_lanes"]) * 2).", found: ".count($fastq_files).")", E_USER_ERROR);
				}

				//copy files
				$target_to_copylines[$tag][] = "\tmkdir -p {$project_folder}Sample_{$sample}";
				foreach ($fastq_files as $fastq_file) 
				{
					if(ends_with(strtolower($fastq_file), ".fastq.ora"))
					{
						//convert to fastq.gz
						$orad_files = [
							"{$project_folder}/Sample_{$sample}/",
							$fastq_file,
							get_path("data_folder")."/dbs/oradata/"
						];
						$out_folder = ["{$project_folder}/Sample_{$sample}/"];
						$orad_command = $parser->execApptainer("orad", "orad", "--ora-reference ".get_path("data_folder")."/dbs/oradata/ ".($overwrite ? " -f" : "")." -t {$threads_ora} -P {$project_folder}/Sample_{$sample}/ ".realpath($fastq_file), $orad_files, $out_folder, true);
						$target_to_copylines[$tag][] = "\t".$orad_command;
					}
					else
					{
						$target_to_copylines[$tag][] = "\tcp ".($overwrite ? "-f " : "")."{$fastq_file} {$project_folder}/Sample_{$sample}/";
					}	
				}

			}
			else
			{
				trigger_error("ERROR: Currently unsupported UMI type '{$umi_type}' provided!", E_USER_ERROR);
			}
		}
		else if(($sys_type == "WGS") || ($sys_type == "WES"))
		{
			//create folder
			$target_to_copylines[$tag][] = "\tmkdir -p {$project_folder}Sample_{$sample}";

			//get FastQs
			$fastq_files = glob($fastq_folder."/{$sample}*_L00[0-9]_R[12]_00[0-9].fastq.{gz,ora}", GLOB_BRACE);
			//check count
			if(count($fastq_files) != count($sample_infos["ps_lanes"]) * 2) 
			{
				trigger_error("Number of FastQ files for sample {$sample} doesn't match number of lanes in run info! (expected: ".(count($sample_infos["ps_lanes"]) * 2).", found in {$fastq_folder}: ".count($fastq_files).")", E_USER_ERROR);
			}

			//check & copy analyzed data
			if ($nsx_analysis_done && !$merge_sample && !$ignore_analysis)
			{
				//check if all files are present
				$tmp = glob("$old_location/{$sample}/*_seq");
				if(count($tmp) == 0) trigger_error("ERROR: Source sample folder '$old_location/{$sample}/*_seq' not found!", E_USER_ERROR);
				if(count($tmp) > 1) trigger_error("ERROR: Multiple folders in '$old_location/{$sample}/'!", E_USER_ERROR);
				$source_folder = $tmp[0];

				//logs & report
				$log_folder = "$old_location/{$sample}/logs";
				if(!file_exists($log_folder)) trigger_error("ERROR: Sample log folder '$log_folder' not found!", E_USER_ERROR);

				$source_mapping_file = "{$source_folder}/{$sample}.bam";
				if(!file_exists($source_mapping_file))
				{
					//No BAM file -> check for CRAM
					$source_mapping_file = "{$source_folder}/{$sample}.cram";
					if(!file_exists($source_mapping_file))
					{
						//No mapping file found
						trigger_error("ERROR: BAM/CRAM file '{$source_mapping_file}/.bam' is missing!", E_USER_ERROR);
					}
					//else: CRAM file was found -> check index
					if(!file_exists($source_mapping_file.".crai")) trigger_error("ERROR: CRAM index file '{$source_mapping_file}.crai' is missing!", E_USER_ERROR);
					
				}
				else
				{
					//BAM file was created -> check index
					if(!file_exists($source_mapping_file.".bai")) trigger_error("ERROR: BAM index file '{$source_mapping_file}.bai' is missing!", E_USER_ERROR);
				} 
				
				if(!$sample_is_tumor) //VC is not done on tumor samples
				{
					$source_vcf_file = "{$source_folder}/{$sample}.hard-filtered.vcf.gz";
					if(!file_exists($source_vcf_file)) trigger_error("ERROR: VCF file '{$source_vcf_file}' is missing!", E_USER_ERROR);
					if(!file_exists($source_vcf_file.".tbi")) trigger_error("ERROR: VCF index file '{$source_vcf_file}.tbi' is missing!", E_USER_ERROR);
					$source_gvcf_file = "{$source_folder}/{$sample}.hard-filtered.gvcf.gz";
					if(!file_exists($source_gvcf_file)) trigger_error("ERROR: gVCF file '{$source_gvcf_file}' is missing!", E_USER_ERROR);
					if(!file_exists($source_gvcf_file.".tbi")) trigger_error("ERROR: gVCF index file '{$source_gvcf_file}.tbi' is missing!", E_USER_ERROR);
					$source_sv_vcf_file = "{$source_folder}/{$sample}.sv.vcf.gz";
					if(!file_exists($source_sv_vcf_file)) trigger_error("ERROR: SV VCF file '{$source_sv_vcf_file}' is missing!", E_USER_ERROR);
					if(!file_exists($source_sv_vcf_file.".tbi")) trigger_error("ERROR: SV VCF index file '{$source_sv_vcf_file}.tbi' is missing!", E_USER_ERROR);
					$source_cnv_vcf_file = "{$source_folder}/{$sample}.cnv.vcf.gz";
					if(($sys_type == "WGS") && !file_exists($source_cnv_vcf_file)) trigger_error("ERROR: CNV VCF file '{$source_cnv_vcf_file}' is missing!", E_USER_ERROR);
					if(($sys_type == "WGS") && !file_exists($source_cnv_vcf_file.".tbi")) trigger_error("ERROR: CNV VCF index file '{$source_cnv_vcf_file}.tbi' is missing!", E_USER_ERROR);
				}

				//get md5sum 
				$md5sum_buffer[] = get_md5_line($source_mapping_file);
				if(!$sample_is_tumor) //VC is not done on tumor samples
				{
					$md5sum_buffer[] = get_md5_line($source_vcf_file);
					$md5sum_buffer[] = get_md5_line($source_gvcf_file);
					// $md5sum_buffer[] = get_md5_line($source_sv_vcf_file);
					if (file_exists($source_cnv_vcf_file)) $md5sum_buffer[] = get_md5_line($source_cnv_vcf_file);
				}


				//*********************** copy analyzed data *********************************************
				$target_to_copylines[$tag][] = "\tmkdir -p {$project_folder}Sample_{$sample}/dragen";
				//copy logs
				$target_to_copylines[$tag][] = "\tcp -r ".($overwrite ? "-f " : "")."{$log_folder} {$project_folder}Sample_{$sample}/dragen/";

				//copy analzed data to dragen folder
				$target_to_copylines[$tag][] = "\tcp -r ".($overwrite ? "-f " : "")."{$source_folder}/* {$project_folder}Sample_{$sample}/dragen/";

				//touch all indices (to prevent warnings)
				$target_to_copylines[$tag][] = "\ttouch {$project_folder}Sample_{$sample}/dragen/*.bai";
				$target_to_copylines[$tag][] = "\ttouch {$project_folder}Sample_{$sample}/dragen/*.crai";
				$target_to_copylines[$tag][] = "\ttouch {$project_folder}Sample_{$sample}/dragen/*.tbi";
				$target_to_copylines[$tag][] = "\ttouch {$project_folder}Sample_{$sample}/dragen/*.csi";

			}
			
			// copy fastq/ORA
			if(($sample_infos["preserve_fastqs"] == 1) || ($project_analysis=="fastq") || !$nsx_analysis_done || $merge_sample || $ignore_analysis)
			{
				foreach ($fastq_files as $fastq_file) 
				{
					if ((ends_with(strtolower($fastq_file), ".fastq.ora") && !$use_dragen && !$merge_sample) || ($project_analysis=="fastq"))
					{
						//convert to fastq.gz
						$orad_files = [
							"{$project_folder}/Sample_{$sample}/",
							$fastq_file,
							get_path("data_folder")."/dbs/oradata/"
						];
						$out_folder = ["{$project_folder}/Sample_{$sample}/"];
						$orad_command = $parser->execApptainer("orad", "orad", "--ora-reference ".get_path("data_folder")."/dbs/oradata/ ".($overwrite ? " -f" : "")." -t {$threads_ora} -P {$project_folder}/Sample_{$sample}/ ".realpath($fastq_file), $orad_files, $out_folder, true);
						$target_to_copylines[$tag][] = "\t".$orad_command;
					}
					else
					{
						$target_to_copylines[$tag][] = "\tcp ".($overwrite ? "-f " : "")."{$fastq_file} {$project_folder}/Sample_{$sample}/";
					}
				}
			}

		}
		else
		{
			trigger_error("ERROR: Analysis other than WES, WGS, cfDNA, Panel or RNA are currently not supported on NovaSeq X!", E_USER_ERROR);
		}
		
	}
	else
	{
		$fastqgz_files = glob($old_location."/Sample_{$sample}/*.fastq.gz");
		$r3_count = 0;
		foreach($fastqgz_files as $file)
		{
			$r3_count += contains($file, "_R3_");
		}
		if (count($fastqgz_files)>=3 && $r3_count==count($fastqgz_files)/3 && !$no_rename_r3) //handling of molecular barcode in index read 2 (HaloPlex HS, Swift, ...)
		{
			//create target folder
			$target_to_copylines[$tag][] = "\tmkdir -p {$sample_folder}";
			$target_to_copylines[$tag][] = "\tchmod 775 {$sample_folder}";
				
			//copy fastq.gz files and change names
			foreach($fastqgz_files as $file) 
			{
				$old_name = basename($file);
				$new_name = strtr($old_name, array("_R2_"=>"_index_", "_R3_"=>"_R2_"));
				if ($sample_infos["sys_umi_type"] === "n/a" && contains($new_name, "_index_"))
				{
					//do not move index/UMI FASTQ files for processing systems with UMI type n/a
					//useful for flow cells with mixed UMI/non-UMI protocols
					trigger_error("Processing system of sample '{$sample}' has UMI type n/a, skipping index bases FASTQ!", E_USER_NOTICE);
					continue;
				}
				$target_to_copylines[$tag][] = "\tmv ".($overwrite ? "-f " : "")."$old_location/Sample_{$sample}/$old_name {$sample_folder}/$new_name";				
			}
		}
		else
		{
			$target_to_copylines[$tag][] = "\tmv ".($overwrite ? "-f " : "")."$old_location/Sample_{$sample}/ {$project_folder}";		
		}
	}

	// make sure all files are accessible by bioinf
	$target_to_copylines[$tag][] = "\tchmod -R 775 {$project_folder}Sample_{$sample}";
	$target_to_copylines[$tag][] = "\tchgrp -R f_ad_bi_l_medgen_access_storages {$project_folder}Sample_{$sample}";

	//skip normal samples which have an associated tumor sample on the same run
	$is_normal_with_tumor = !$sample_is_tumor && isset($normal2tumor[$sample]);
	
	//additional arguments for db_queue_analysis
	$args = array();

	if($project_analysis!="fastq") //if more than FASTQ creation should be done for samples's project
	{	
		//queue analysis
		if ($sample_is_tumor && $sys_type!="RNA")
		{
			//queue tumor, with somatic specific options
			//add variant calling for diagnostic normal samples
			$outputline = "php {$repo_folder}/src/Tools/db_queue_analysis.php -type 'single sample' -samples {$sample}";
			if($is_novaseq_x && $nsx_analysis_done)
			{
				$outputline .= " -args '-steps db -somatic'";
			} 
			else 
			{
				$outputline .= " -args '-steps ma,db -somatic'";
				//use DRAGEN mapping:
				$outputline .= ($use_dragen ? " -use_dragen": "");
			}
			
			if (isset($tumor2normal[$sample]))
			{
				$normal = $tumor2normal[$sample];				
				$normal_processed_info = get_processed_sample_info($db_conn,$normal);
				$n_dir = $normal_processed_info["ps_folder"];
				//queue somatic analysis: only if normal sample is included in this run or its already folder exists
				if(in_array($normal, array_keys($sample_data)) || is_dir($n_dir))
				{
					$queue_somatic[] = "\tphp {$repo_folder}/src/Tools/db_queue_analysis.php -type 'somatic' -samples {$sample} {$normal} -info tumor normal ".($use_dragen ? "-use_dragen": "");
				}
				else
				{
					trigger_error("Skipping {$sample}-{$normal} somatic analysis because normal folder $n_dir does not exist!",E_USER_WARNING);
				}
			}
			else
			{
				//queue tumor-only somatic analysis
				$queue_somatic[] = "\tphp {$repo_folder}/src/Tools/db_queue_analysis.php -type 'somatic' -samples {$sample} -info tumor";
			}
		}
		else
		{
			//default analysis
			$outputline = "php {$repo_folder}/src/Tools/db_queue_analysis.php -type 'single sample' -samples {$sample}";

			//determine analysis steps from project
			$args = [];
			if ($project_analysis=="mapping")
			{
				if ($nsx_analysis_done && (($sys_type == "WGS") || ($sys_type == "WES") && !$merge_sample)) $args[] = "-steps db"; //mapping already done on device
				else $args[] = "-steps ma,db";
			}
			else if ($project_analysis=="variants")
			{
				if (($sys_type == "WGS") || ($sys_type == "WES"))
				{
					if ($nsx_analysis_done && !$merge_sample && !$ignore_analysis)
					{
						$args[] = "-steps vc,cn,sv,re,db";
					}
					else
					{
						if ($use_dragen)
						{
							//perform DRAGEN analysis at front:
							$outputline .= " -use_dragen";
						}
						//complete analysis => no '-steps' parameter 
					}
				}
				//other sample types > use all default steps
			}
			
			//add somatic specific options if it is the normal to a tumor sample
			if ($is_normal_with_tumor)
			{
				$args[] = "-somatic";
			}
			
			//check if trio needs to be queued
			$parents = get_trio_parents($sample);
			if (!is_null($parents))
			{
				list($father, $mother) = $parents;
				$command = "php {$repo_folder}/src/Tools/db_queue_analysis.php -type 'trio' -samples {$sample} {$father} {$mother} -info child father mother";
				if ($high_priority) $outputline .= " -high_priority";
				$queue_trios[] = "\t".$command;
			}
		}
		if ($high_priority || $sample_infos['urgent']==1)
		{
			$outputline .= " -high_priority";
		}
		
		if(count($args)>0)
		{
			$outputline .= " -args '".implode(" ", $args)."'";
		}
		
		$target_to_queuelines[$tag][]="\t".$outputline;
	}
	else if(!$is_normal_with_tumor)
	{		
		$project_to_fastqonly_samples[$project_name][] = $sample;
	}
}

//create Makefile
$output = array();
$all_parts = array();
$output[] = "all: ".(($is_novaseq_x)?"":"chmod ")."import_runqc import_read_counts ";
$output[] = "";

if (!$is_novaseq_x) //not required for NovaSeqX run
{
	//target 'chmod'
	$output[] = "chmod:";
	$output[] = "\tchmod -R 775 {$folder}";
	$output[] = "";
}

//target 'import_runqc'
$output[] = "import_runqc:";
$output[] = "\tphp {$repo_folder}/src/Tools/runqc_parser.php -name \"{$run_name}\" -run_dir $folder/../ -force";
$output[] = "";

//target 'import_read_counts'
$output[] = "import_read_counts:";
if ($is_novaseq_x)
{
	$output[] = "\tphp {$repo_folder}/src/IMGAG/import_sample_read_counts.php -csv_mode -stats {$analysis_id}/Data/Demux/Demultiplex_Stats.csv -db $db ";
}
else
{
	$output[] = "\tphp {$repo_folder}/src/IMGAG/import_sample_read_counts.php -stats $folder/Stats/Stats.json -db $db ";
}
$output[] = "";

if ($is_novaseq_x)
{
	//check md5 checksums of analysis data
	if (count($md5sum_buffer) > 0)
	{
		//write checksums to file
		file_put_contents("analysis_checksums.txt", implode("\n", $md5sum_buffer));

		$output[] = "check_md5_analysis_data:";
		$output[] = "\tmd5sum -c analysis_checksums.txt";

		$output[] = "";

		$all_parts[] = "check_md5_analysis_data";
	}


	//check ORA files
	$ora_files = glob($folder."/[0-9]/Data/*/ora_fastq/*.fastq.ora");
	if (count($ora_files) > 0)
	{
		$output[] = "check_ora_data:";
		//add check for ORA files
		foreach($ora_files as $ora_file)
		{
			$ora_file = realpath($ora_file);
			$in_files = array($ora_file, get_path("data_folder")."/dbs/oradata/");
			$orad_command = $parser->execApptainer("orad", "orad", "--ora-reference ".get_path("data_folder")."/dbs/oradata/ -t {$threads_ora} -C ".$ora_file, $in_files, [], true);
			$output[] = "\t".$orad_command;
		}

		$output[] = "";

		$all_parts[] = "check_ora_data";
	}
}

//target 'create_project_folder'
if (count($create_project_folder) > 0)
{
	$create_project_folder = array_unique($create_project_folder);
	$output[] = "create_project_folder:";
	foreach ($create_project_folder as $line) 
	{
		$output[] = $line;
	}
	$output[] = "";
}



//group commands by sample
$urgent_sample_buffer = array();
$normal_sample_buffer = array();
$skipped_queuing_commands = array();
$urgent_samples = array();
$normal_samples = array();

foreach($sample_data as $sample => $sample_infos)
{
	$sample_target = array();

	$sample_target[] = "Sample_{$sample}:";

	//copy commands
	$sample_target[] = "\t#copy";
	$copy_cmds = $target_to_copylines[$sample];
	if (count($copy_cmds) < 1) trigger_error("No data found to copy for Sample {$sample}!", E_USER_ERROR);
	foreach ($copy_cmds as $line) 
	{
		$sample_target[] = $line;
	}

	//merge
	if (isset($merge_files[$sample]))
	{
		$sample_target[] = "\t#merge";
		$sample_target[] = $merge_files[$sample];
	}

	//queue
	if (isset($target_to_queuelines[$sample]))
	{
		if (!$no_queuing)
		{
			$sample_target[] = "\t#queue";
			foreach ($target_to_queuelines[$sample] as $line) 
			{
				$sample_target[] = $line;
			}
		}
		else
		{
			foreach ($target_to_queuelines[$sample] as $line) 
			{
				$skipped_queuing_commands[] = $line;
			}
		}
	}

	$sample_target[] = "";

	//divide samples in urgent and non urgent
	if ($sample_infos['urgent']==1)	
	{
		$urgent_samples[] = "Sample_".$sample;
		$urgent_sample_buffer[] = implode("\n", $sample_target);
	}
	else 
	{
		$normal_samples[] = "Sample_".$sample;
		$normal_sample_buffer[] = implode("\n", $sample_target);
	}
	
}
//add to Makefile
$all_parts = array_merge($all_parts, $urgent_samples, $normal_samples);
$output = array_merge($output, $urgent_sample_buffer, $normal_sample_buffer);

//target 'queue_somatic'
if (count($queue_somatic)>0)
{
	if (!$no_queuing) $all_parts[] = "queue_somatic";
	
	$output = array_merge($output, ["queue_somatic:"], array_unique($queue_somatic));
	$output[] = "";	
}

//target 'queue_trios'
if (count($queue_trios)>0)
{
	if (!$no_queuing) $all_parts[] = "queue_trios";
	
	$output = array_merge($output, ["queue_trios:"], array_unique($queue_trios));
	$output[] = "";	
}


//target(s) 'email_...'
foreach ($project_to_fastqonly_samples as $target => $samples)
{
	$all_parts[] = "email_{$target}";
	
	$output[] = "email_{$target}:";
	$output[] = "\t".create_mail_command($db_conn, $target, $samples);
	$output[] = "";
}

//add targets to 'all' target
$output[0] .= implode(" ", $all_parts);

//write data
file_put_contents($out, implode("\n", $output));

?>
