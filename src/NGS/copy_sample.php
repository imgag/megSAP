<?php 
/** 
	@page copy_sample
*/

$basedir = dirname($_SERVER['SCRIPT_FILENAME'])."/../";
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("copy_sample", "\$Rev: 899 $", "Creates a Makefile to copy de-multiplexed sample data to projects and queue data analysis.");
$parser->addInfile("ss",  "Input samplesheet that was used to create the data folder.", false);
$parser->addString("fo",  "Input data folder.", true, "Unaligned");
$parser->addOutfile("out",  "Output Makefile. Default: 'Makefile'.", true);
extract($parser->parse($argv));

function array_not_unique($input) 
{
	$duplicates=array();
	$processed=array();
	foreach($input as $i) 
	{
		if(in_array($i,$processed)) 
		{
			$duplicates[]=$i;
		} 
		else 
		{
			$processed[]=$i;
		}
	}
	return $duplicates;
}

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
		print ("Warning: processed sample name \"$processed_sample_name\" is not of form SAMPLENAME_PROCESSID, Skipped it\n");
		return array("SKIPPED", "SKIPPED", "SKIPPED", "SKIPPED");
	}
	$process_id=ltrim($process_id,'0');
	$db_connect = DB::getInstance('NGSD');
	
	//get tumour status
	$result = $db_connect->executeQuery("SELECT tumor FROM  sample WHERE name='".$sample_name."'");
	$tumor_status=$result[0]['tumor'];
	
	//get sample id
	$result = $db_connect->executeQuery("SELECT id FROM  sample WHERE name='".$sample_name."'");
	$sample_id=$result[0]["id"];
	
	//get processed sample processing_system and project ID
	$result = $db_connect->executeQuery("SELECT project_id,processing_system_id FROM  processed_sample WHERE sample_id=".$sample_id." AND process_id= ".$process_id);
	$projectID=$result[0]["project_id"];
	$processing_systemID=$result[0]["processing_system_id"];
	
	//check if processing system type is haloplex hs
	$result = $db_connect->executeQuery("SELECT type FROM  processing_system WHERE id=".$processing_systemID);
	$haloplex_hs=($result[0]["type"]=="Panel Haloplex HS");	
	
	//get project name, type, coordinator_id and analysis step
	$result = $db_connect->executeQuery("SELECT  name,type, analysis,internal_coordinator_id FROM project WHERE id=".$projectID);
	$projectname=$result[0]["name"];
	$project_type=$result[0]["type"];
	$project_analysis=$result[0]["analysis"];
	$project_coord_id=$result[0]["internal_coordinator_id"];
	
	//get coordinator name and email
	$result = $db_connect->executeQuery("SELECT name,email FROM user WHERE id=".$project_coord_id);
	$internal_coord=array("name" => $result[0]["name"], "email" => $result[0]["email"]);
	//get run ID
	$result = $db_connect->executeQuery("SELECT sequencing_run_id FROM  processed_sample WHERE sample_id=".$sample_id." AND process_id= ".$process_id);
	$run_ID=$result[0]["sequencing_run_id"];
	
	//get run name
	$result = $db_connect->executeQuery("SELECT name FROM  sequencing_run WHERE id=".$run_ID);
	$run_name=$result[0]["name"];
	return array($projectname, $project_type, $run_name, $tumor_status, $project_analysis, $internal_coord, $haloplex_hs);
}

function create_mail_command($coordinator, $email_ad, $samples, $project_name)
{
	$mail_text=array();
	$mail_text[]="Hallo $coordinator,";
	$mail_text[]="";
	$mail_text[]="die FASTQ-Dateien der folgenden Proben Deines Projekts $project_name liegen vor:";
	$mail_text[]="";
	foreach(array_unique($samples) as $sample)
	{
		$mail_text[]=$sample;
	}
	$mail_text[]="";
	$mail_text[]="Viele Gruesse";
	$mail_text[]="";
	$mail_text[]="der automatischer Mailversender";
	return "php -r 'mail(\"$email_ad\",\"Neue Daten fuer $project_name\", \"".implode("\\n",$mail_text)."\",\"Reply-To: Florian.Lenz@med.uni-tuebingen.de\");'";
}

function create_import_run_qc_command($run_name, $runfolder)
{
	return "import_runqc:\n\tphp /mnt/users/all/php/src/NGS/runqc_parser.php -name \"$run_name\" -run_dir $runfolder/ -force -db NGSD";
}

//write makefile lines for file and folder operations and stores them in a dictionary
function build_makefile($folder, $sample_IDs, $sample_projectname_map, $sample_projecttype_map, $project_coord_map,  $sample_tumor_status_map, $runnumber, $makefile_name, $basedir,$nxtSeq, $sample_analysis_step_map,$haloplex_hs_map,$import_run_qc_command)
{
	if (!file_exists($folder))
	{
		trigger_error("Folder '$folder' does not exist!", E_USER_ERROR);
	}
	
	$samples_to_merge=array_not_unique($sample_IDs); //samples that appears more than once in run should be proccesed using the "-merge" flag
	$outfile_lines=array();
	$outfile_lines[]="make:\n";
	$target_to_copylines=array();
	$target_to_queuelines=array();
	$skipped_lines=array();//stores lines to print the sample names which are not parseable and was not processed
	$skipped_lines[]="skipped:";
	$project_to_fastqonly_samples=array();
	
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
		
		//copy_target line if first sample of project in run, mkdir and chmod  if first sample of project at all
		if (!(array_key_exists($project_name."_".$project_type, $target_to_copylines)))
		{
			$outputline="\n"."\n"."copy_".$project_name."_".$project_type.":";
			$target_to_copylines[$project_name."_".$project_type]=array($outputline);
			if (!is_dir($new_location)) 
			{
				$outputline="mkdir -p ".$new_location."/";
				$target_to_copylines[$project_name."_".$project_type][]="\t".$outputline;
				$outputline="chmod 775 ".$new_location."/";
				$target_to_copylines[$project_name."_".$project_type][]="\t".$outputline;
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
		if ($haloplex_hs_map[$sample_ID])
		{
			//names of the FASTQ-Files need to be changed
			$files = glob($old_location."/Sample_".$sample_ID."/*.fastq.gz");
			foreach($files as $file) 
			{
				$target_to_copylines[$project_name."_".$project_type][]="\tmkdir -p ".$new_location."/Sample_".$sample_ID."/";
				$target_to_copylines[$project_name."_".$project_type][]="\tchmod 775 ".$new_location."/Sample_".$sample_ID."/";
				$old_name=basename($file);
				$replace_count=0;
				$new_name=str_replace("_R2_","_index_",$old_name,$replace_count);
				if ($replace_count)	
				{	
					$target_to_copylines[$project_name."_".$project_type][]="\tcp -i -p ".$old_location."/Sample_".$sample_ID."/$old_name ".$new_location."/Sample_".$sample_ID."/$new_name";
				}
				else
				{
					$new_name=str_replace("_R3_","_R2_",$old_name,$replace_count);
					if ($replace_count)
					{
						$target_to_copylines[$project_name."_".$project_type][]="\tcp -i -p ".$old_location."/Sample_".$sample_ID."/$old_name ".$new_location."/Sample_".$sample_ID."/$new_name";
					}
					else
					{
						$target_to_copylines[$project_name."_".$project_type][]="\tcp -i -p ".$old_location."/Sample_".$sample_ID."/$old_name ".$new_location."/Sample_".$sample_ID."/$old_name";
					}
				}
				
			}
		}
		else
		{
			$target_to_copylines[$project_name."_".$project_type][]="\tcp -i -r -p ".$old_location."/Sample_".$sample_ID."/ ".$new_location."/";
		}

		if($sample_analysis_step_map[$sample_ID]!="fastq") //if more than FASTQ creation should be done for samples's project
		{					
			//build target lines for analysis using Sungrid Engine's queues if first sample of project on this run
			if (!(array_key_exists($project_name."_".$project_type, $target_to_queuelines)))
			{
				$outputline="\n"."\n"."queue_".$project_name."_".$project_type.":";//build target in the output makefile
				$target_to_queuelines[$project_name."_".$project_type]=array($outputline);
			}
			
			//build  first part of line for analysis using Sungrid Engine's queues,
			$outputline= "php /mnt/users/all/php/src/NGS/queue_sample.php -sample ".$sample_ID;

			//set high_priority if diagnostic sample
			if ($project_type=='diagnostic') $outputline.=" -high_priority ";
			
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
						
			$target_to_queuelines[$project_name."_".$project_type][]="\t".$outputline." ";
		}
		else
		{		
			$project_to_fastqonly_samples[$project_name][]=$sample_ID;
		}
	}
	
	//build "all"-target line and chmod target
	$all_line="all: chmod_target import_runqc ";
	foreach (array_keys($target_to_copylines) as $target)
	{
		$all_line .= "copy_".$target." ";
		
	}
	foreach (array_keys($target_to_queuelines) as $target)
	{
		$all_line .= "queue_".$target." ";	
	}
	foreach (array_keys($project_to_fastqonly_samples) as $target)
	{
		$all_line .= "email_".$target." ";	
	}
	$all_line.="skipped";
	
	//save the "all" line to a new Makefile
	file_put_contents($makefile_name, implode ("\n",array($all_line, "","chmod_target:", "\tchmod -R 775 {$folder}/")));
	
	//add the other lines to the Makefile, minus the duplicates that may occur due to identical lines in Sample Sheet
	foreach ($target_to_copylines as $target => $lines)
	{
		//using "array_unique" to remove duplicates caused  by duplicate samples in SampleSheet which got the "-merge" flag
		file_put_contents($makefile_name, implode ("\n", array_unique($target_to_copylines[$target])), FILE_APPEND);
		if (isset($target_to_queuelines[$target]))
		{
			file_put_contents($makefile_name, implode ("\n", array_unique($target_to_queuelines[$target])), FILE_APPEND);
		}		
	}
	
	foreach ($project_to_fastqonly_samples as $project => $samples)
	{
		file_put_contents($makefile_name, "\nemail_".$project.":\n", FILE_APPEND);
		file_put_contents($makefile_name, "\t".create_mail_command($project_coord_map[$project]["name"], $project_coord_map[$project]["email"], $samples, $project)."\n", FILE_APPEND);
	}
	file_put_contents($makefile_name,"\n", FILE_APPEND);
	file_put_contents($makefile_name, $import_run_qc_command, FILE_APPEND);
	file_put_contents($makefile_name,"\n", FILE_APPEND);
	file_put_contents($makefile_name, implode ("\n", $skipped_lines), FILE_APPEND);
	file_put_contents($makefile_name,"\n", FILE_APPEND);
}

if (!isset($out)) $out = "Makefile";
list($sample_IDs, $nxtSeq)=extract_IDs_from_samplesheet($ss);
$sample_projectname_map=array();
$sample_projecttype_map=array();
$sample_tumor_status_map = array();
$old_runnumber=-1;
$project_coord_map=array();
$haloplex_hs_map=array();
foreach($sample_IDs as $sampleID)
{	
	list($sample_projectname_map[$sampleID],$sample_projecttype_map[$sampleID],$runnumber,$sample_tumor_status_map[$sampleID], $sample_analysis_step_map[$sampleID], $internal_coord, $haloplex_hs_map[$sampleID])=get_parameters($sampleID);
	if ($sample_projectname_map[$sampleID]=="SKIPPED")//unable to extract all information due to malformed Sample ID
	{
		$runnumber=$old_runnumber;//set run number (which is "SKIPPED" now) back to last value
		continue;
	}
	if (($old_runnumber!=$runnumber) and ($old_runnumber!=-1))
	{
		trigger_error("Inconsistent run numbers within samples found: ".$old_runnumber." and ".$runnumber, E_USER_ERROR);
	}
	$old_runnumber=$runnumber;
	$project_coord_map[$sample_projectname_map[$sampleID]]=$internal_coord;
}
$import_run_qc_command=create_import_run_qc_command($runnumber,realpath(dirname($ss)));
build_makefile($fo, $sample_IDs, $sample_projectname_map, $sample_projecttype_map, $project_coord_map, $sample_tumor_status_map, $runnumber,$out,$basedir,$nxtSeq, $sample_analysis_step_map,$haloplex_hs_map,$import_run_qc_command);

?>
