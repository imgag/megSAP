<?php 
/** 
	@page copy_sample
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");
$repo_folder = "/mnt/users/all/megSAP"; //fixed absolute path to make the tests work for all users

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("copy_sample", "Creates a Makefile to copy de-multiplexed sample data to projects and queue data analysis.");
$parser->addString("samplesheet",  "Input samplesheet that was used to create the data folder.", true, "SampleSheet_bcl2fastq.csv");
$parser->addString("folder",  "Input data folder.", true, "Unaligned");
$parser->addString("runinfo" ,"Illumina RunInfo.xml file. Necessary for checking metadata in NGSD",true, "RunInfo.xml");
$parser->addOutfile("out",  "Output Makefile. Default: 'Makefile'.", true);
$parser->addFlag("high_priority", "Assign high priority to all queued samples.");
$parser->addFlag("overwrite", "Do not prompt before overwriting FASTQ files.");
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
extract($parser->parse($argv));

//extract samples names and sequencer type from sample sheet
function extract_sample_data(&$db_conn, $filename)
{
	//determine type
	$file = file($filename);
	if (trim($file[0])=="[Data]")//NxtSeq run
	{
		$is_nextseq = true;
		array_shift($file);//skip [Data]
		array_shift($file);//skip header
	}
	else
	{
		$is_nextseq = false;
		array_shift($file);//skip header
	}
	
	//extract sample data
	$sample_data = array();
	foreach($file as $line)
	{
		if (trim($line)=="") continue;
		
		if ($is_nextseq)
		{
			list(, $name_with_prefix) = explode(",", $line);
			$sample = substr($name_with_prefix, 7); //remove "Sample_"
		}
		else
		{
			list(, , $sample) = explode(",", $line);
		}
		$sample = trim($sample);
		
		$sample_data[$sample] = get_processed_sample_info($db_conn, $sample);
	}
	
	return array($sample_data, $is_nextseq);
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
		trigger_error("Could not find any sample that uses $lane_count lanes as specified in {$run_info_xml_file}. Please check number of lanes in database and {$sample_sheet}!",E_USER_WARNING);
	}
	return true;
}

//init
if (!isset($out)) $out = "Makefile";
$db_conn = DB::getInstance($db);

if(!check_number_of_lanes($runinfo,$samplesheet))
{
	trigger_error("Could not verify number of lanes used for Demultiplexing and actually used on Sequencer. Please check manually!",E_USER_WARNING);
}

//get sample data from samplesheet
list($sample_data, $is_nextseq) = extract_sample_data($db_conn, $samplesheet);

//create Illumina-like folder for MGI data
if (file_exists("Fq"))
{
	print "Detected that this is an MGI run folder!\n";
	if (file_exists($folder))
	{
		print "Data folder '$folder' exists => skipping creation of Illumina-like folder at '$folder'.\n";
	}
	else
	{
		print "Creating Illumina-like folder at '$folder'...\n";
		$copied_fastqs = [];
		foreach($sample_data as $sample => $data)
		{
			$sample_folder = "{$folder}/".$data['project_name']."/Sample_{$sample}/";
			print "$sample ($sample_folder)\n";
			
			//create sample folder
			exec2("mkdir -p $sample_folder");
			
			//determine data for the next step
			$fcid = $data["run_fcid"];
			$lanes = $data["ps_lanes"];
			$mids = array($data["ps_mid1"]);
			$reads = array("1");
			if (count(explode("+", $data['run_recipe']))>2) $reads[] = "2";
			foreach($data['ps_comments'] as $comment_line)
			{
				if (starts_with($comment_line, "add_mids:"))
				{
					$add_mids = explode(",", substr($comment_line, strpos($comment_line, ":")+1));
					foreach($add_mids as $mid)
					{
						$mids[] = trim($mid);
					}
				}
			}
			for ($i=0; $i<count($mids); ++$i)
			{
				$mids[$i] = trim(strtr($mids[$i], array("MGI"=>"")));
			}
			
			//copy FASTQs
			foreach($lanes as $lane)
			{
				foreach($mids as $mid)
				{
					foreach($reads as $read)
					{
						$from = "Fq/L0{$lane}/{$fcid}_L0{$lane}_{$mid}_{$read}.fq.gz";
						$to = "{$sample_folder}/{$sample}_L00{$lane}_MID{$mid}_R{$read}_001.fastq.gz";
						print "  Copying $from (".filesize_mb($from)." MB)\n";
						$parser->copyFile($from, $to);
						$copied_fastqs[] = $from;
					}
				}
			}
		}
		
		//check for large FASTQs that were not copied
		list($big_fastqs) = exec2("find Fq -name '*.fq.gz' -size +50M | sort");
		$unprocessed_fastqs = array_diff($big_fastqs, $copied_fastqs);
		if (count($unprocessed_fastqs)>0)
		{
			print "Unprocessed FASTQ files larger than 50MB:\n";
			foreach($unprocessed_fastqs as $fastq)
			{
				print "  $fastq (".filesize_mb($fastq)." MB)\n";
			}
		}
	}
}

//check that data folder exists
if (!file_exists($folder))
{
	trigger_error("Data folder '$folder' does not exist!", E_USER_ERROR);
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
	
	//get run name (the same for all samples)
	$run_name = $sample_infos['run_name'];
}
$queued_normal_samples = [];

//parse input
$target_to_copylines = array();
$target_to_queuelines = array();
$project_to_fastqonly_samples = array();
$sample_to_newlocation = array();
foreach($sample_data as $sample => $sample_infos)
{
	$project_name = $sample_infos['project_name'];
	$project_type = $sample_infos['project_type'];
	$project_folder = $sample_infos['project_folder'];
	$project_analysis = $sample_infos['project_analysis'];
	$sample_folder = $sample_infos['ps_folder'];
	$sample_is_tumor = $sample_infos['is_tumor'];
	$sys_type = $sample_infos['sys_type'];
	$sys_target = $sample_infos['sys_target'];

	//calculate current location of sample
	if ($is_nextseq)
	{
		$old_location = "{$folder}/".$project_name;
	}
	else
	{
		$old_location = "{$folder}/Project_".$project_name;
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
	$fastqgz_files = glob($old_location."/Sample_{$sample}/*.fastq.gz");
	$r3_count = 0;
	foreach($fastqgz_files as $file)
	{
		$r3_count += contains($file, "_R3_");
	}
	if (count($fastqgz_files)>=3 && $r3_count==count($fastqgz_files)/3 ) //handling of molecular barcode in index read 2 (HaloPlex HS, Swift, ...)
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
			$outputline = "php {$repo_folder}/src/NGS/db_queue_analysis.php -type 'single sample' -samples {$sample} -args '-steps ma,db -somatic'";
			$outputline .= "\n\t";

			if (isset($tumor2normal[$sample]))
			{
				$normal = $tumor2normal[$sample];
				//queue normal if on same run, with somatic specific options
				if (!in_array($normal, $queued_normal_samples)  && array_key_exists($normal,$sample_data) && $sample_data[$normal]["run_name"] === $sample_infos["run_name"])
				{
					$steps = ($project_type=="diagnostic") ? "ma,vc,an,cn,db" : "ma,db";
					$outputline .= "php {$repo_folder}/src/NGS/db_queue_analysis.php -type 'single sample' -samples {$normal} -args '-steps $steps -somatic'";
					$outputline .= "\n\t";
					//track that normal sample is queued
					$queued_normal_samples[] = $normal;
				}
				
				$normal_processed_info = get_processed_sample_info($db_conn,$normal);
				$n_dir = $normal_processed_info["ps_folder"];
				//queue somatic analysis: only if normal sample is included in this run or its folder exists
				if(in_array($normal,$queued_normal_samples) || is_dir($n_dir))
				{
					$outputline .= "php {$repo_folder}/src/NGS/db_queue_analysis.php -type 'somatic' -samples {$sample} {$normal} -info tumor normal";
				}
				else
				{
					trigger_error("Skipping {$tumor}-{$normal} somatic analysis because normal folder $n_dir does not exist!",E_USER_WARNING);
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
			$outputline = "php {$repo_folder}/src/NGS/db_queue_analysis.php -type 'single sample' -samples {$sample}";

			//determine analysis steps from project
			if ($project_analysis=="mapping")
			{
				$args[] = "-steps ma,db";
			}
			if ($project_analysis=="variant calling")
			{
				$args[] = "-steps ma,vc,db,cn";
			}
		}
		
		if ($high_priority)
		{
			$args[] = "-high_priority";
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
$output[] = "all: chmod import_runqc ";
$output[] = "";

//target 'chmod'
$output[] = "chmod:";
$output[] = "\tchmod -R 775 {$folder}";
$output[] = "";

//target 'import_runqc'
$output[] = "import_runqc:";
$output[] = "\tphp {$repo_folder}/src/NGS/runqc_parser.php -name \"{$run_name}\" -run_dir $folder/../ -force";
$output[] = "";

//target(s) 'copy_...'
$all_parts = array();
foreach ($target_to_copylines as $target => $lines)
{
	$all_parts[] = "copy_{$target}";
	
	$output = array_merge($output, array_unique($lines));
	$output[] = "";
}
	
//target(s) 'queue_...'
foreach ($target_to_queuelines as $target => $lines)
{
	$all_parts[] = "queue_{$target}";
	
	$output = array_merge($output, array_unique($lines));
	$output[] = "";	
}

//target(s) 'qc_...'
foreach ($project_to_fastqonly_samples as $target => $samples)
{
	$all_parts[] = "qc_{$target}";
	
	$output[] = "qc_{$target}:";
	foreach (array_unique($samples) as $sample)
	{
		$sample_location = $sample_to_newlocation[$sample];
		$output[] = "\tphp {$repo_folder}/src/NGS/qc_fastq.php -folder {$sample_location} -import";
	}
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