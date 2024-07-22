<?php 
/** 
	@page copy_sample
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");


error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("copy_sample", "Creates a Makefile to copy de-multiplexed sample data to projects and queue data analysis.");
$parser->addString("samplesheet",  "Input samplesheet that was used to create the data folder. (Has to be located in the run folder.)", true, "SampleSheet_bcl2fastq.csv");
$parser->addString("folder",  "Input data folder.", true, "Unaligned");
$parser->addString("runinfo" ,"Illumina RunInfo.xml file. Necessary for checking metadata in NGSD",true, "RunInfo.xml");
$parser->addOutfile("out",  "Output Makefile. Default: 'Makefile'.", true);
$parser->addFlag("high_priority", "Assign high priority to all queued samples.");
$parser->addFlag("overwrite", "Do not prompt before overwriting FASTQ files.");
$parser->addFlag("no_rename_r3", "Do not rename R2/R3 FASTQ files to index/R2.");
$parser->addFlag("ignore_nsx_analysis", "Ignore NovaSeqX Analysis results.");
$parser->addFlag("no_genlab", "Do not query GENLAB for metadata.");
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
$parser->addInt("threads_ora", "Number of threads used to decompress ORA files during the copy.", true, 8);
$parser->addFlag("test", "Run in test mode, e.g. set the pipeline path to a fixed value.");
extract($parser->parse($argv));

//extract samples names and sequencer type from sample sheet
function extract_sample_data(&$db_conn, $filename)
{
	$file = file($filename);
	if(starts_with($file[0], "[Data]")) array_shift($file); //NovaSeq 6000: skip "[Data]"
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
	return "php -r 'mail(\"$email\",\"Neue Daten fuer $project_name\", \"".implode("\\n",$mail_text)."\",\"Reply-To: ".get_path("queue_email")."\");'";
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
	$sample_sheet_content = Matrix::fromTSV($sample_sheet);
	$i_lane = $sample_sheet_content->getColumnIndex("Lane",false,false);
	if($i_lane == -1) return false;
	
	$lane_count_found_in_samplesheet = false;
	for($i=0;$i<$sample_sheet_content->rows();++$i)
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

//init
if (!isset($out)) $out = "Makefile";
$db_conn = DB::getInstance($db);
$ngsbits = get_path("ngs-bits");
if($test) $repo_folder = "/mnt/storage2/megSAP/pipeline"; //fixed absolute path to make the tests work for all users
else $repo_folder = repository_basedir(); //use repositories tool in production

$file_acccess_group = get_path("file_access_group", false);
if ($file_acccess_group == "") trigger_error("File access group not set in the settings.ini. File group will not be set!", E_USER_WARNING);



//fallback for NovaSeq X default SampleSheet name
if (($samplesheet == "SampleSheet_bcl2fastq.csv") && !file_exists($samplesheet)) $samplesheet = "SampleSheet.csv";
//check SampleSheet
if(!file_exists($samplesheet)) trigger_error("ERROR: Provided SampleSheet '{$samplesheet}' not found!", E_USER_ERROR);

//get file names
$run_folder = dirname(realpath($samplesheet));

//get run id and check FlowCell id
$run_folder_parts = explode("_", basename($run_folder));
$run_name = $run_folder_parts[count($run_folder_parts) - 1];
if((strlen($run_name) != 5) || ((int) $run_name > 9999) || ((int)$run_name < 0)) trigger_error("ERROR: invalid run folder suffix '{$run_name}' provided!", E_USER_ERROR);
//add '#' prefix
$run_name = "#".$run_name;

//get flowcell id
if(!file_exists($runinfo)) trigger_error("ERROR: Run info file '{$runinfo}' not found!", E_USER_ERROR);
$xml = simplexml_load_file($runinfo);
if(empty($xml->Run->Flowcell)) trigger_error("ERROR: Run info file doesn't contain flowcell id!", E_USER_ERROR);
$flowcell_id = $xml->Run->Flowcell;

//check flowcell id
$flowcell_id_ngsd = $db_conn->getValue("SELECT fcid FROM sequencing_run WHERE `name`='{$run_name}'");
if($flowcell_id != $flowcell_id_ngsd) trigger_error("ERROR: FlowCell id from run info doesn't match FlowCell id in NGSD! (Run info: {$flowcell_id} <-> NGSD: {$flowcell_id_ngsd})", E_USER_ERROR);

//get instrument type
$run_parameters_xml = $run_folder."/RunParameters.xml";
if(!file_exists($run_parameters_xml)) trigger_error("ERROR: Required RunParameters.xml file is missing in the run folder!", E_USER_ERROR);
$is_novaseq_x = is_novaseq_x_run($run_parameters_xml) && !$ignore_nsx_analysis;

//change default data folder for NovaSeqX 
if($is_novaseq_x && ($folder=="Unaligned")) 
{
	trigger_error("NovaSeqX run detected!", E_USER_NOTICE);
	$folder="Analysis";
}

//check that data folder exists
if (!file_exists($folder)) trigger_error("Data folder '$folder' does not exist!", E_USER_ERROR);

if(!file_exists("Fq") && !check_number_of_lanes($runinfo,$samplesheet))
{
	trigger_error("***!!!WARNING!!!***\nCould not verify number of lanes used for Demultiplexing and actually used on Sequencer. Please check manually!", E_USER_WARNING);
}

//get sample data
if($is_novaseq_x) $sample_data = get_sample_data_from_db($db_conn, $run_name);
else $sample_data = extract_sample_data($db_conn, $samplesheet);


//import data from Genlab
if (! $no_genlab)
{
	print "Importing information from GenLab...\n";

	foreach($sample_data as $sample => $data)
	{
		$args = [];
		$args[] = "-ps {$sample}";
		if ($db=="NGSD_TEST") $args[] = "-test";
		$parser->exec("{$ngsbits}/NGSDImportGenlab", implode(" ", $args), true);
	}
}


//update sample data after importing sample relations
if($is_novaseq_x) $sample_data = get_sample_data_from_db($db_conn, $run_name);
else $sample_data = extract_sample_data($db_conn, $samplesheet);

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
	if($run_name != $sample_infos['run_name']) trigger_error("ERROR: Sequencing run doesn't match sample info ('".$sample_infos['run_name']."')");
}
$queued_normal_samples = [];



//Check for former run that contains more than 50% same samples and offer merging to user
$former_run_name = "";
$merge_files = array(); //contains merge commands, commands are inserted before queue
$use_dragen = false; //parameter to enable DRAGEN mapping
$skip_queuing = false; //parameter to skip queuing if run is first genome run
$wgs_use_dragen_data = false; //parameter to use created Dragen data for analysis

if($run_name != "")
{
	$former_run_name =  check_former_run($db_conn, $run_name);
	if($former_run_name != "")
	{
		//skip user interaction on test runs:
		$answer = "n";
		if($db != "NGSD_TEST")
		{
			echo "Former run '{$former_run_name}' detected. Merge (y/n)?\n";
			$answer = trim(fgets(STDIN));
		}
		
		if($answer == "y")
		{
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
					$merge_files[] = "\tphp {$repo_folder}/src/Tools/merge_samples.php -ps ".$old_samples[$current_key]." -into $current_ps";
				}
			}

			// check if DRAGEN can be used
			if ((get_path("dragen_user", false) != "") && (get_path("queues_dragen", false) != ""))
			{
				// ask if DRAGEN should be used
				echo "Should the genome samples on this run mapped with DRAGEN mapping? (y/n)?\n";
				if(trim(fgets(STDIN)) == "y") $use_dragen = true;
				
			}

		}
		elseif($answer != "n")
		{
			trigger_error("Please choose whether files from both runs shall be merged.", E_USER_ERROR);
		}
	}
	else
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
				$skip_queuing = true;
			} 
			else
			{
				if ($is_novaseq_x)
				{
					echo "The genomes are already analyzed on the device. Should this data be used for the analysis? (y/n)?\n";
					if(trim(fgets(STDIN)) == "y") $wgs_use_dragen_data = true;
				}
			}
		}
		elseif (($n_genomes > 0) && $is_novaseq_x)
		{
			echo "There are genomes on this run which are already analyzed on the device. Should this data be used for the analysis? (y/n)?\n";
			if(trim(fgets(STDIN)) == "y") $wgs_use_dragen_data = true;
		}

	}
}
else
{
	trigger_error("Could not determine current run name. Skipping check for former run that could be merged into current run.", E_USER_WARNING);
}

//get analysis folder
if ($is_novaseq_x)
{
	$analyses = array_filter(glob($folder."/[0-9]"), 'is_dir');
	if(count($analyses) == 0) trigger_error("ERROR: No analysis found!", E_USER_ERROR);
	if(count($analyses) > 1) trigger_error("ERROR: Multiple analyses found!", E_USER_ERROR);
	$analysis_id = $analyses[0];
}

//parse input
$target_to_copylines = array();
$target_to_queuelines = array();
$queue_trios = array();
$project_to_fastqonly_samples = array();
$sample_to_newlocation = array();

if ((get_path("dragen_user", false) != "") && (get_path("queues_dragen", false) != "") && !$test)
{
	$use_dragen_somatic = null;
}
else
{
	$use_dragen_somatic = false;
}

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
	$project_folder = $sample_infos['project_folder'];
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
			//create function name for first sample
			if (count($merge_files) == 0) $merge_files[] = "merge:";
			$merge_files[] = "\tphp {$repo_folder}/src/Tools/merge_samples.php -ps {$source_ps_name} -into {$sample}";
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
		if($sys_type == "WGS")
		{
			$old_location .= "/DragenGermline";
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
			$old_location .= "/DragenEnrichment";
			if (file_exists($old_location."/ora_fastq"))
			{
				$fastq_folder = $old_location."/ora_fastq";
			}
			else
			{
				$fastq_folder = $old_location."/fastq";
			}
		}
		else if(($sys_type == "RNA") || ($sys_type == "Panel") || ($sys_type == "cfDNA (patient-specific)") || ($sys_type == "cfDNA") || ($sys_type == "WGS (shallow)"))
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
		else trigger_error("ERROR: Invalid processing system type '{$sys_type}' for NovaSeq X analysis!", E_USER_ERROR);
	}
	else
	{
		$old_location = "{$folder}/{$project_name}";
	}
	
	//determine project 
	$sample_to_newlocation[$sample] = $sample_folder;
	
	//copy_target line if first sample of project in run, mkdir and chmod if first sample of project at all
	$tag = $project_name."_".$project_type;
	if (!isset($target_to_copylines[$tag]))
	{
		$target_to_copylines[$tag] = array("copy_{$tag}:");
		if (!is_dir($project_folder)) 
		{
			$target_to_copylines[$tag][] = "\tmkdir -p {$project_folder}";
			$target_to_copylines[$tag][] = "\tchmod 775 {$project_folder}";
		}
	}
	
	//build copy line
	if($is_novaseq_x)
	{
		if(($sys_type == "RNA") || ($sys_type == "Panel") || ($sys_type == "cfDNA (patient-specific)") || ($sys_type == "cfDNA") || ($sys_type == "WGS (shallow)"))
		{
			//only copy FastQ files
			//get FastQs
			$fastq_files = glob($fastq_folder."/{$sample}*_L00[0-9]_R[123]_00[0-9].fastq.{gz,ora}", GLOB_BRACE);
			if($umi_type == "n/a" || $umi_type == "IDT-UDI-UMI" || $umi_type == "IDT-xGen-Prism" || $umi_type == "Twist")
			{
				//no index files: simply copy/convert FastQs
				
				//check count
				if(count($fastq_files) != count($sample_infos["ps_lanes"]) * 2) 
				{
					trigger_error("ERROR: Number of FastQ files for sample {$sample} doesn't match number of lanes in run info! (expected: ".(count($sample_infos["ps_lanes"]) * 2).", found: ".count($fastq_files).")", E_USER_ERROR);
				}

				//copy files
				$target_to_copylines[$tag][] = "\tmkdir -p {$project_folder}Sample_{$sample}";
				foreach ($fastq_files as $fastq_file) 
				{
					if(ends_with(strtolower($fastq_file), ".fastq.ora"))
					{
						//convert to fastq.gz
						$target_to_copylines[$tag][] = "\t".get_path("orad")." --ora-reference ".dirname(get_path("orad"))."/oradata/ ".($overwrite ? " -f" : "")." -t {$threads_ora} -P {$project_folder}/Sample_{$sample}/ {$fastq_file}";
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
			//WES/WGS sample ==> use full analysis (ma,vc,sv)

			//if not manually set, ignore analysis of WGS samples since the WGS samples will be mapped by the DRAGEN server
			if($sys_type != "WGS" || $wgs_use_dragen_data)
			{
				//check if all files are present
				$tmp = glob("$old_location/{$sample}/*_seq");
				if(count($tmp) == 0) trigger_error("ERROR: Source sample folder '$old_location/{$sample}/*_seq' not found!", E_USER_ERROR);
				if(count($tmp) > 1) trigger_error("ERROR: Multiple folders in '$old_location/{$sample}/'!", E_USER_ERROR);
				$source_folder = $tmp[0];

				//logs & report
				$log_folder = "$old_location/{$sample}/logs";
				if(!file_exists($log_folder)) trigger_error("ERROR: Sample log folder '$log_folder' not found!", E_USER_ERROR);
				$report_file = $source_folder."/report.html";
				if(!file_exists($report_file)) trigger_error("ERROR: Report HTML file '{$report_file}' not found!", E_USER_ERROR);

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
				
				$source_vcf_file = "{$source_folder}/{$sample}.hard-filtered.vcf.gz";
				if(!file_exists($source_vcf_file)) trigger_error("ERROR: VCF file '{$source_vcf_file}' is missing!", E_USER_ERROR);
				if(!file_exists($source_vcf_file.".tbi")) trigger_error("ERROR: VCF index file '{$source_vcf_file}.tbi' is missing!", E_USER_ERROR);
				$source_gvcf_file = "{$source_folder}/{$sample}.hard-filtered.gvcf.gz";
				if(!file_exists($source_gvcf_file)) trigger_error("ERROR: gVCF file '{$source_gvcf_file}' is missing!", E_USER_ERROR);
				if(!file_exists($source_gvcf_file.".tbi")) trigger_error("ERROR: gVCF index file '{$source_gvcf_file}.tbi' is missing!", E_USER_ERROR);
				$source_sv_vcf_file = "{$source_folder}/{$sample}.sv.vcf.gz";
				if(!file_exists($source_sv_vcf_file)) trigger_error("ERROR: SV VCF file '{$source_sv_vcf_file}' is missing!", E_USER_ERROR);
				if(!file_exists($source_sv_vcf_file.".tbi")) trigger_error("ERROR: SV VCF index file '{$source_sv_vcf_file}.tbi' is missing!", E_USER_ERROR);
			}
			

			//get FastQs
			$fastq_files = glob($fastq_folder."/{$sample}*_L00[0-9]_R[12]_00[0-9].fastq.{gz,ora}", GLOB_BRACE);
			//check count
			if(count($fastq_files) != count($sample_infos["ps_lanes"]) * 2) 
			{
				trigger_error("ERROR: Number of FastQ files for sample {$sample} doesn't match number of lanes in run info! (expected: ".(count($sample_infos["ps_lanes"]) * 2).", found: ".count($fastq_files).")", E_USER_ERROR);
			}

			$move_cmd = "mv ".($overwrite ? "-f " : "");
			
			//create folder
			$target_to_copylines[$tag][] = "\tmkdir -p {$project_folder}Sample_{$sample}";
			//ignore analysis of WGS samples for 2-run WGS samples
			if(!$merge_sample && ($sys_type != "WGS" || $wgs_use_dragen_data))
			{
				$target_to_copylines[$tag][] = "\tmkdir -p {$project_folder}Sample_{$sample}/dragen_variant_calls";
				//copy logs
				$target_to_copylines[$tag][] = "\tcp -r {$log_folder} {$project_folder}Sample_{$sample}/dragen_variant_calls/";
				$target_to_copylines[$tag][] = "\tcp {$report_file} {$project_folder}Sample_{$sample}/dragen_variant_calls/logs";

				//move BAM
				$target_to_copylines[$tag][] = "\t{$move_cmd} {$source_mapping_file} {$project_folder}Sample_{$sample}/";
				if (file_exists($source_mapping_file.".bai"))
				{
					$target_to_copylines[$tag][] = "\t{$move_cmd} {$source_mapping_file}.bai {$project_folder}Sample_{$sample}/";
				}
				else
				{
					$target_to_copylines[$tag][] = "\t{$move_cmd} {$source_mapping_file}.crai {$project_folder}Sample_{$sample}/";
				}
				
				if(!$sample_is_tumor)
				{
					//move (g)VCFs
					$target_to_copylines[$tag][] = "\t{$move_cmd} {$source_vcf_file} {$project_folder}Sample_{$sample}/dragen_variant_calls/{$sample}_dragen.vcf.gz";
					$target_to_copylines[$tag][] = "\t{$move_cmd} {$source_vcf_file}.tbi {$project_folder}Sample_{$sample}/dragen_variant_calls/{$sample}_dragen.vcf.gz.tbi";
					$target_to_copylines[$tag][] = "\t{$move_cmd} {$source_gvcf_file} {$project_folder}Sample_{$sample}/dragen_variant_calls/{$sample}_dragen.gvcf.gz";
					$target_to_copylines[$tag][] = "\t{$move_cmd} {$source_gvcf_file}.tbi {$project_folder}Sample_{$sample}/dragen_variant_calls/{$sample}_dragen.gvcf.gz.tbi";
					//move SV VCFs
					$target_to_copylines[$tag][] = "\t{$move_cmd} {$source_sv_vcf_file} {$project_folder}Sample_{$sample}/dragen_variant_calls/{$sample}_dragen_svs.vcf.gz";
					$target_to_copylines[$tag][] = "\t{$move_cmd} {$source_sv_vcf_file}.tbi {$project_folder}Sample_{$sample}/dragen_variant_calls/{$sample}_dragen_svs.vcf.gz.tbi";
				}
			}

			if(($sample_infos["preserve_fastqs"] == 1) || ($project_analysis=="fastq") || ($sys_type == "WGS" && !$wgs_use_dragen_data) || $merge_sample)
			{
				foreach ($fastq_files as $fastq_file) 
				{
					if(ends_with(strtolower($fastq_file), ".fastq.ora"))
					{
						//convert to fastq.gz
						$target_to_copylines[$tag][] = "\t".get_path("orad")." --ora-reference ".dirname(get_path("orad"))."/oradata/".($overwrite ? " -f" : "")." -t {$threads_ora} -P {$project_folder}/Sample_{$sample}/ {$fastq_file}";
					}
					else
					{
						$target_to_copylines[$tag][] = "\tcp ".($overwrite ? "-f " : "")."{$fastq_file} {$project_folder}/Sample_{$sample}/";
					}
					// // make sure all files are accessible by bioinf
					// $target_to_copylines[$tag][] = "\tchmod 777 {$project_folder}Sample_{$sample}/*.fastq.gz";	
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
	if ($file_acccess_group != "") $target_to_copylines[$tag][] = "\tchgrp -R {$file_acccess_group} {$project_folder}Sample_{$sample}";

	//skip normal samples which have an associated tumor sample on the same run
	$is_normal_with_tumor = !$sample_is_tumor && isset($normal2tumor[$sample]);
	
	//additional arguments for db_queue_analysis
	$args = array();

	if($project_analysis!="fastq" && !$is_normal_with_tumor) //if more than FASTQ creation should be done for samples's project
	{	
		//build target lines for analysis using Sungrid Engine's queues if first sample of project on this run
		if (!array_key_exists($tag, $target_to_queuelines))
		{
			$target_to_queuelines[$tag] = array("queue_{$tag}:");
		}
		
		//queue analysis
		if ($sample_is_tumor && $sys_type!="RNA")
		{
			//queue tumor, with somatic specific options
			//add variant calling for diagnostic normal samples
			$outputline = "php {$repo_folder}/src/NGS/db_queue_analysis.php -type 'single sample' -samples {$sample}";
			if($is_novaseq_x)
			{
				$outputline .= " -args '-steps db -somatic'\n\t";
			} 
			else 
			{
				$outputline .= " -args '-steps ma,db -somatic'\n\t";
			}
			if (isset($tumor2normal[$sample]))
			{
				$normal = $tumor2normal[$sample];
				//queue normal if on same run, with somatic specific options
				if (!in_array($normal, $queued_normal_samples)  && array_key_exists($normal,$sample_data) && $sample_data[$normal]["run_name"] === $sample_infos["run_name"])
				{
					$steps = ($is_novaseq_x ? "vc,cn,db" : "ma,vc,cn,db");
					$outputline .= "php {$repo_folder}/src/NGS/db_queue_analysis.php -type 'single sample' -samples {$normal}";
					if($is_novaseq_x)
					{
						$outputline .= " -args '-steps {$steps} -use_dragen -somatic'\n\t";
					}
					else
					{
						$outputline .= " -args '-steps {$steps} -somatic'\n\t";
					}
					//track that normal sample is queued
					$queued_normal_samples[] = $normal;
				}
				
				$normal_processed_info = get_processed_sample_info($db_conn,$normal);
				$n_dir = $normal_processed_info["ps_folder"];
				//queue somatic analysis: only if normal sample is included in this run or its folder exists
				if(in_array($normal,$queued_normal_samples) || is_dir($n_dir))
				{
					if ($use_dragen_somatic === null)
					{
						echo "Somatic anylsis about to be started. Use dragen for somatic analysis on this run (y/n)?\n";
						$answer = trim(fgets(STDIN));
						
						$use_dragen_somatic = ($answer == "y");
					}
					
					if ($use_dragen_somatic)
					{
						$outputline .= "php {$repo_folder}/src/NGS/db_queue_analysis.php -type 'somatic' -samples {$sample} {$normal} -info tumor normal -args '-use_dragen'";
					}
					else
					{
						$outputline .= "php {$repo_folder}/src/NGS/db_queue_analysis.php -type 'somatic' -samples {$sample} {$normal} -info tumor normal";
					}
				}
				else
				{
					trigger_error("Skipping {$sample}-{$normal} somatic analysis because normal folder $n_dir does not exist!",E_USER_WARNING);
				}
			}
			else
			{
				//queue tumor-only somatic analysis
				$outputline .= "php {$repo_folder}/src/NGS/db_queue_analysis.php -type 'somatic' -samples {$sample} -info tumor";
			}
		}
		else
		{
			//TODO: check if sample will be merged => redo (mapping & ) VC
			$outputline = "php {$repo_folder}/src/NGS/db_queue_analysis.php -type 'single sample' -samples {$sample}";

			// add DRAGEN parameter
			if ($use_dragen && !$is_novaseq_x)
			{
				$processed_sample_info = get_processed_sample_info($db_conn,$sample);
				if ($processed_sample_info['sys_type'] == "WGS") $args[] = "-use_dragen -no_abra";
			}

			//determine analysis steps from project
			if ($project_analysis=="mapping")
			{
				$args[] = (($is_novaseq_x && ($sys_type != "RNA") && ($sys_type != "cfDNA (patient-specific)") && ($sys_type != "cfDNA") && ($sys_type != "WGS (shallow)") && ($sys_type != "Panel") && !$merge_sample) ? "-steps db -use_dragen" : "-steps ma,db");
			}
			if ($project_analysis=="variants")
			{
				//no steps parameter > use all default steps
				if($is_novaseq_x && ($sys_type != "RNA") && ($sys_type != "cfDNA (patient-specific)") && ($sys_type != "cfDNA") && ($sys_type != "WGS (shallow)") && ($sys_type != "Panel")) // do  complete analysis for RNA/cfDNA samples 
				{
					//do mapping for WGS samples
					if (($sys_type == "WGS" && !$wgs_use_dragen_data) || $merge_sample) $args[] = "-steps ma,vc,cn,sv,re,db";
					else $args[] = "-steps vc,cn,sv,re,db";
					$args[] = "-use_dragen";
					$args[] = "-no_abra";
				}
			}
			
			//check if trio needs to be queued
			$parents = get_trio_parents($sample);
			if (!is_null($parents))
			{
				list($father, $mother) = $parents;
				$command = "php {$repo_folder}/src/NGS/db_queue_analysis.php -type 'trio' -samples {$sample} {$father} {$mother} -info child father mother";
				if ($high_priority) $outputline .= " -high_priority";
				$queue_trios[] = "\t".$command;
			}
		}
		if ($high_priority || contains(strtolower(implode(" ", $sample_infos['ps_comments'])), "eilig"))
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
$output[] = "\tphp {$repo_folder}/src/NGS/runqc_parser.php -name \"{$run_name}\" -run_dir $folder/../ -force";
$output[] = "";

//target 'import_read_counts'
$output[] = "import_read_counts:";
if ($is_novaseq_x)
{
	$output[] = "\tphp {$repo_folder}/src/NGS/import_sample_read_counts.php -csv_mode -stats ${analysis_id}/Data/Demux/Demultiplex_Stats.csv -db $db ";
}
else
{
	$output[] = "\tphp {$repo_folder}/src/NGS/import_sample_read_counts.php -stats $folder/Stats/Stats.json -db $db ";
}
$output[] = "";

//target(s) 'copy_...'
$all_parts = array();
foreach ($target_to_copylines as $target => $lines)
{
	$all_parts[] = "copy_{$target}";
	
	$output = array_merge($output, array_unique($lines));
	$output[] = "";
}

//target 'merge'
if(count($merge_files) > 0)
{
	$all_parts[] = "merge";
	$output = array_merge($output, $merge_files);
	$output[] = "";

	//report merged samples
	print(implode("\n", $merge_notice));
}
	
//target(s) 'queue_...'
if (!$skip_queuing)
{
	foreach ($target_to_queuelines as $target => $lines)
	{
		$all_parts[] = "queue_{$target}";
		
		$output = array_merge($output, array_unique($lines));
		$output[] = "";	
	}
}


//target(s) 'queue_trios'
if (!$skip_queuing)
{
	if (count($queue_trios)>0)
	{
		$all_parts[] = "queue_trios";
		
		$output = array_merge($output, ["queue_trios:"], array_unique($queue_trios));
		$output[] = "";	
	}
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
