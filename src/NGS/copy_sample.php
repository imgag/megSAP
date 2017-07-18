<?php 
/** 
	@page copy_sample
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");
$repo_folder = "/mnt/users/all/megSAP"; //fixed absolute path to make the tests work

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("copy_sample", "Creates a Makefile to copy de-multiplexed sample data to projects and queue data analysis.");
$parser->addString("samplesheet",  "Input samplesheet that was used to create the data folder.", true, "SampleSheet_bcl2fastq.csv");
$parser->addString("folder",  "Input data folder.", true, "Unaligned");
$parser->addOutfile("out",  "Output Makefile. Default: 'Makefile'.", true);
extract($parser->parse($argv));

//finds and returns the SampleIDs from a Sample Sheet
function extract_IDs_from_samplesheet($samplesheet_name)
{
	$samplesheet_file=file($samplesheet_name);
	if (trim($samplesheet_file[0])=="[Data]")//if NxtSeq 500 run
	{
		array_shift($samplesheet_file);//skip [Data]
		array_shift($samplesheet_file);//skip header
		foreach($samplesheet_file as $line)
		{
			if (trim($line)=="") continue;
			list(, $nameWithPrefix) = explode(",", $line);
			$sampleIDs[] = substr($nameWithPrefix, strlen("Sample_"));//remove "Sample_"
		}
		return array($sampleIDs,true);//return sampleID and information that it's a NxtSeq
	}
	else
	{
		array_shift($samplesheet_file);//skip header
		$sampleIDs=array();
		
		foreach($samplesheet_file as $line)
		{
			list(, , $sampleIDs[]) = explode(",",$line);
		}
		return array($sampleIDs,false);//return sampleID and information that it's not a NxtSeq
	}
}

//find the project name and process sample type for a sample ID, find and check run number
function get_parameters($processed_sample_name)
{	
	//get sample name and process ID
	if (strpos($processed_sample_name,'_') !== false)
	{
		list($sample_name,$process_id) = explode("_",$processed_sample_name);
	}
	else
	{
		print "Warning: processed sample name \"$processed_sample_name\" is not of form SAMPLENAME_PROCESSID, Skipped it\n";
		return array("SKIPPED", "SKIPPED", "SKIPPED", "SKIPPED");
	}
	$process_id = ltrim($process_id,'0');
	$db_connect = DB::getInstance('NGSD');
	
	//get tumour status
	$result = $db_connect->executeQuery("SELECT tumor FROM  sample WHERE name='".$sample_name."'");
	$tumor_status = $result[0]['tumor'];
	
	//get sample id
	$sample_id = $db_connect->getId("sample", "name", $sample_name);
	
	//get processed sample processing_system and project ID
	$result = $db_connect->executeQuery("SELECT id,project_id,processing_system_id FROM  processed_sample WHERE sample_id=".$sample_id." AND process_id= ".$process_id);
	$psample_id = $result[0]["id"];
	$projectID = $result[0]["project_id"];
	$processing_systemID = $result[0]["processing_system_id"];
	
	//get project name, type, coordinator_id and analysis step
	$result = $db_connect->executeQuery("SELECT  name,type, analysis,internal_coordinator_id FROM project WHERE id=".$projectID);
	$projectname = $result[0]["name"];
	$project_type = $result[0]["type"];
	$project_analysis = $result[0]["analysis"];
	$project_coord_id = $result[0]["internal_coordinator_id"];
	
	//get coordinator name and email
	$result = $db_connect->executeQuery("SELECT name,email FROM user WHERE id=".$project_coord_id);
	$internal_coord=array("name" => $result[0]["name"], "email" => $result[0]["email"]);
	//get run ID
	$result = $db_connect->executeQuery("SELECT sequencing_run_id FROM  processed_sample WHERE sample_id=".$sample_id." AND process_id= ".$process_id);
	$run_ID = $result[0]["sequencing_run_id"];
	
	//get run name
	$result = $db_connect->executeQuery("SELECT name FROM  sequencing_run WHERE id=".$run_ID);
	$run_name = $result[0]["name"];
	
	//check for associated tumor sample(s)
	$result = $db_connect->executeQuery("SELECT id FROM processed_sample WHERE normal_id=".$psample_id);
	$assoc_tumor = array_column($result, "id");
	
	return array($projectname, $project_type, $run_name, $tumor_status, $project_analysis, $internal_coord, $assoc_tumor);
}

function create_mail_command($coordinator, $email_ad, $samples, $project_name)
{
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
	return "php -r 'mail(\"$email_ad\",\"Neue Daten fuer $project_name\", \"".implode("\\n",$mail_text)."\",\"Reply-To: ".get_path("queue_email")."\");'";
}

//write makefile lines for file and folder operations and stores them in a dictionary
function build_makefile($folder, $sample_IDs, $sample_projectname_map, $sample_projecttype_map, $project_coord_map,  $sample_tumor_status_map, $sample_assoc_tumor_map, $runnumber, $makefile_name, $repo_folder, $nxtSeq, $sample_analysis_step_map)
{
	if (!file_exists($folder))
	{
		trigger_error("Folder '$folder' does not exist!", E_USER_ERROR);
	}
	
	//init
	$output = array();
	$target_to_copylines = array();
	$target_to_queuelines = array();
	$skipped_lines = array();//stores lines to print the sample names which are not parseable and was not processed
	$skipped_lines[] = "skipped:";
	$project_to_fastqonly_samples = array();
	$sample_to_newlocation = array();
	
	//parse input
	foreach($sample_IDs as $sample_ID)
	{
		if ($sample_projectname_map[$sample_ID]=="SKIPPED")//if the correct parameters are not available
		{
			$skipped_lines[]="\t"."@echo '$sample_ID was skipped!'";
			continue;
		}
		$project_name=$sample_projectname_map[$sample_ID];
		$project_type=$sample_projecttype_map[$sample_ID];
		//build line to copy sample
		//calculate current location of sample
		if ($nxtSeq)
		{
			$old_location= "{$folder}/".$project_name;
		}
		else
		{
			$old_location= "{$folder}/Project_".$project_name;
		}
		
		//determine project 
		$new_location = get_path("project_folder").$project_type."/".$project_name;
		$sample_to_newlocation[$sample_ID] = $new_location."/Sample_".$sample_ID;
		
		//copy_target line if first sample of project in run, mkdir and chmod  if first sample of project at all
		$tag = $project_name."_".$project_type;
		if (!(array_key_exists($tag, $target_to_copylines)))
		{
			$outputline="copy_{$tag}:";
			$target_to_copylines[$tag]=array($outputline);
			if (!is_dir($new_location)) 
			{
				$outputline="mkdir -p ".$new_location."/";
				$target_to_copylines[$tag][]="\t".$outputline;
				$outputline="chmod 775 ".$new_location."/";
				$target_to_copylines[$tag][]="\t".$outputline;
			}
		}
		
		//bcl2fastq2 changes "-" to "_" in sample name, revert that
		if ($nxtSeq)
		{
			$sample_ID_modified = strtr($sample_ID, "_", "-");
			$old_files = glob("{$old_location}/Sample_{$sample_ID}/{$sample_ID_modified}*.fastq.gz");
			foreach($old_files as $file)
			{
				$file_corrected = strtr($file, array($sample_ID_modified => $sample_ID));
				rename($file, $file_corrected);
			}
		}
		
		//build copy line
		$fastqgz_files = glob($old_location."/Sample_".$sample_ID."/*.fastq.gz");
		$r3_count = 0;
		foreach($fastqgz_files as $file)
		{
			$r3_count += contains($file, "_R3_");
		}
		if (count($fastqgz_files)>=3 && $r3_count==count($fastqgz_files)/3 ) //handling of molecular barcode in index read 2 (Haloplex HS, Swift, ...)
		{
			//create target folder
			$target_to_copylines[$tag][]="\tmkdir -p ".$new_location."/Sample_".$sample_ID."/";
			$target_to_copylines[$tag][]="\tchmod 775 ".$new_location."/Sample_".$sample_ID."/";
				
			//copy fastq.gz files and change names
			foreach($fastqgz_files as $file) 
			{
				$old_name = basename($file);
				$new_name = strtr($old_name, array("_R2_"=>"_index_", "_R3_"=>"_R2_"));
				$target_to_copylines[$tag][]="\tcp -i ".$old_location."/Sample_".$sample_ID."/$old_name ".$new_location."/Sample_".$sample_ID."/$new_name";				
			}
		}
		else
		{
			$target_to_copylines[$tag][]="\tcp -i -r ".$old_location."/Sample_".$sample_ID."/ ".$new_location."/";
		}

		//skip normal samples for SomaticAndTreatment project
		$is_normal_with_tumor= ($sample_tumor_status_map[$sample_ID] == 0) && !empty($sample_assoc_tumor_map[$sample_ID]);
		
		if($sample_analysis_step_map[$sample_ID]!="fastq" && !$is_normal_with_tumor) //if more than FASTQ creation should be done for samples's project
		{					
			//build target lines for analysis using Sungrid Engine's queues if first sample of project on this run
			if (!(array_key_exists($tag, $target_to_queuelines)))
			{
				$outputline="queue_{$tag}:";//build target in the output makefile
				$target_to_queuelines[$tag]=array($outputline);
			}
			
			//build  first part of line for analysis using Sungrid Engine's queues,
			$outputline= "php {$repo_folder}/src/NGS/queue_sample.php -sample ".$sample_ID;

			//stop at mapping if analysis for project is set to mapping
			if ($sample_analysis_step_map[$sample_ID]=="mapping")
			{
				$outputline.=" -steps ma,db";
			}
			
			//stop at variant calling if analysis for project is set to variant calling
			if ($sample_analysis_step_map[$sample_ID]=="variant calling")
			{
				$outputline.=" -steps ma,vc,db,cn";
			}
						
			$target_to_queuelines[$tag][]="\t".$outputline." ";
		}
		elseif(!$is_normal_with_tumor)
		{		
			$project_to_fastqonly_samples[$project_name][] = $sample_ID;
		}
	}
	
	//target 'all'
	$all_line="all: chmod import_runqc ";
	foreach ($target_to_copylines as $target => $lines)
	{
		$all_line .= "copy_".$target." ";
	}
	foreach ($target_to_queuelines as $target => $lines)
	{
		$all_line .= "queue_".$target." ";	
	}
	foreach ($project_to_fastqonly_samples as $target => $samples)
	{
		$all_line .= "qc_".$target." ";
	}
	foreach ($project_to_fastqonly_samples as $target => $samples)
	{
		$all_line .= "email_".$target." ";
	}
	$all_line .= "skipped";
	$output[] = $all_line;
	$output[] = "";
	
	//target 'chmod'
	$output[] = "chmod:";
	$output[] = "\tchmod -R 775 {$folder}";
	$output[] = "";
	
	//target 'import_runqc'
	$output[] = "import_runqc:";
	$output[] = "\tphp {$repo_folder}/src/NGS/runqc_parser.php -name \"$runnumber\" -run_dir $folder/../ -force -db NGSD";
	$output[] = "";
	
	//target(s) 'copy_...'
	foreach ($target_to_copylines as $target => $lines)
	{
		$output = array_merge($output, array_unique($lines));
		$output[] = "";
	}
		
	//target(s) 'queue_...'
	foreach ($target_to_queuelines as $target => $lines)
	{
		$output = array_merge($output, array_unique($lines));
		$output[] = "";	
	}

	//target(s) 'qc_...'
	foreach ($project_to_fastqonly_samples as $project => $samples)
	{
		$output[] = "qc_".$project.":";
		foreach (array_unique($samples) as $sample) {
			$sample_location = $sample_to_newlocation[$sample];
			$output[] = "\tphp {$repo_folder}/src/NGS/qc_fastq.php -folder {$sample_location} -import";
		}
		$output[] = "";
	}
	
	//target(s) 'email_...'
	foreach ($project_to_fastqonly_samples as $project => $samples)
	{
		$output[] = "email_".$project.":";
		$output[] = "\t".create_mail_command($project_coord_map[$project]["name"], $project_coord_map[$project]["email"], $samples, $project);
		$output[] = "";
	}
	
	//target 'skipped'
	$output = array_merge($output, $skipped_lines);
	$output[] = "";
	
	//write data
	file_put_contents($makefile_name, implode("\n", $output));
}

if (!isset($out)) $out = "Makefile";
list($sample_IDs, $nxtSeq) = extract_IDs_from_samplesheet($samplesheet);
$sample_projectname_map = array();
$sample_projecttype_map = array();
$sample_tumor_status_map = array();
$sample_assoc_tumor_map = array();
$old_runnumber = -1;
$project_coord_map = array();
foreach($sample_IDs as $sampleID)
{	
	list($sample_projectname_map[$sampleID], $sample_projecttype_map[$sampleID], $runnumber, $sample_tumor_status_map[$sampleID], $sample_analysis_step_map[$sampleID], $internal_coord, $sample_assoc_tumor_map[$sampleID]) = get_parameters($sampleID);
	if ($sample_projectname_map[$sampleID]=="SKIPPED")//unable to extract all information due to malformed Sample ID
	{
		$runnumber = $old_runnumber;//set run number (which is "SKIPPED" now) back to last value
		continue;
	}
	if ($old_runnumber!=$runnumber && $old_runnumber!=-1)
	{
		trigger_error("Inconsistent run numbers within samples found: ".$old_runnumber." and ".$runnumber, E_USER_ERROR);
	}
	$old_runnumber = $runnumber;
	$project_coord_map[$sample_projectname_map[$sampleID]] = $internal_coord;
}

build_makefile($folder, $sample_IDs, $sample_projectname_map, $sample_projecttype_map, $project_coord_map, $sample_tumor_status_map, $sample_assoc_tumor_map, $runnumber, $out, $repo_folder, $nxtSeq, $sample_analysis_step_map);

?>