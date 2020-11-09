<?php
/** 
	@page qbic_copy
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("qbic_copy", "Copies QBIC data into the QBIC datamover folder.");
$parser->addFlag("upload", "Enable real upload (otherwise a dry run is performed: dummy data is written to the temporary folder, but it is not copied to the datamover folder).");
$parser->addString("project", "Restrict upload to a project.", true, "");
$parser->addStringArray("samples", "Restrict upload to a list of processed sample.", true, "");
$parser->addFlag("force_reupload", "Upload files even if already uploaded.", true, "");
$parser->addFlag("ignore_quality", "Upload all samples regardless of quality.");
extract($parser->parse($argv));

//check project
function project_names()
{
	$db = DB::getInstance("NGSD");
	$res = $db->executeQuery("SELECT name FROM project");
	$output = array();
	foreach($res as $r)
	{
		$output[] = $r['name'];
	}
	
	sort($output);
	return $output;
}
if ($project!="")
{
	$projects = project_names();
	if (!in_array($project, $projects))
	{
		trigger_error("Invalid project name '$project'! Valid project names are: ".implode(", ", $projects), E_USER_ERROR);
	}
}

//check that datamover is running
print "##Checking that datamover is running...\n";
print "##\n";
list($jobs) = exec2("ps aux | grep datamover | grep bioinf");
function contains_java($str) { return contains($str, "java"); };
$jobs = array_filter($jobs, "contains_java");
if (count($jobs)!=1)
{
	print "##\n";
	print "##============================ ERROR ============================\n";
	print "##=    No datamover of the user bioinf is running! Execute:     =\n";
	print "##=    > sudo su bioinf                                         =\n";
	print "##=    > /mnt/share/to_qbic/datamover.sh start                  =\n";
	print "##===============================================================\n";
	die(1);
}

//print documentation
$GLOBALS["datamover_done"] = "/mnt/share/to_qbic/data/done.txt";
$GLOBALS["datamover_path"] = "/mnt/share/to_qbic/data/to_qbic/";
$GLOBALS["datamover_tmp"] = "/mnt/share/to_qbic/data/tmp/";
print "##Datamover folder : ".$GLOBALS["datamover_path"]."\n";
print "##Datemover log    : /mnt/share/to_qbic/log/datamover_log.txt\n";
print "##Done samples file: ".$GLOBALS["datamover_done"]."\n";
print "##Temporary folder : ".$GLOBALS["datamover_tmp"]."\n";
if (!$upload)
{
	print "##\n";
	print "##========================== ATTENTION ==========================\n";
	print "##=         This is a dry run. No data is uploaded!             =\n";
	print "##=       Use the -upload flag to really update data!           =\n";
	print "##===============================================================\n";
}

//determine QBIC/FO ID from Probeneingang
function checkProbeneingang($name, $table)
{
	$name = strtr($name, array("FO-"=>"FO"));
	
	$db = DB::getInstance("Probeneingang");
	$res = $db->executeQuery("SELECT identifier, identifier_external FROM $table WHERE identifier LIKE '%$name%'");
	foreach($res as $row)
	{
		$ext = trim($row['identifier_external']);
		if (strlen($ext)==10 && $ext[0]=="Q")
		{
			return array(trim($row['identifier']), $ext);
		}
	}
	
	return null;
}

//determine QBIC/FO ID from sample name (and external sample name)
function getExternalNames($name, $name_ex)
{
	$qbic_name = "";
	$fo_name = "";
	
	$names = explode(" ", strtr($name_ex, ",;", "  "));
	foreach($names as $name)
	{
		$name = trim($name);
		if (strlen($name)==10 && $name[0]=="Q")
		{
			$qbic_name = $name;
		}
		
		//Check (derived) samples in Probeneingang
		if (starts_with($name, "FO") && strlen($name)>6)
		{
			$tmp = checkProbeneingang($name, "sample");
			if (!is_null($tmp))
			{
				list($fo_name, $qbic_name) = $tmp;			
			}
			$tmp = checkProbeneingang($name, "derived_sample");
			if (!is_null($tmp))
			{
				list($fo_name, $qbic_name) = $tmp;			
			}
		}
	}
	
	if ($qbic_name=="") return null;
	
	return array($qbic_name, $fo_name);
}

//determine samples meta data
function getSampleInfo($ps_id)
{
	$db = DB::getInstance("NGSD");
	$res = $db->executeQuery("SELECT ps.process_id, sys.name_manufacturer, gen.build, ps.normal_id, ps.quality, s.sample_type, s.quality as quality2, s.tumor, s.name as name2, s.name_external FROM processed_sample ps, processing_system sys, genome gen, sequencing_run run, sample s WHERE ps.id='$ps_id' AND ps.sample_id=s.id AND ps.processing_system_id=sys.id AND sys.genome_id=gen.id");
	list ($ps_num, $ps_sys, $ps_genome, $ps_normal, $ps_qual, $s_type, $s_qual, $s_tumor, $s_name, $s_name_ex) = array_values($res[0]);
	$ps_name = $s_name."_".str_pad($ps_num, 2, '0', STR_PAD_LEFT);
	
	$output = array();
	$output['ngsd_processedsample_id'] = $ps_id;
	$output['id_genetics'] = $ps_name;
	$names = getExternalNames($s_name, $s_name_ex);
	$output['id_qbic'] = $names[0];
	$output['processing_system'] = $ps_sys;
	$output['tumor'] = $s_tumor ? "yes" : "no";
	$output['genome'] = $ps_genome;
	$output['quality_sample'] = $s_qual;
	$output['quality_processed_sample'] = $ps_qual;
	$output['qc'] = array();
	$output['normal_id'] = $ps_normal;
	
	//get run quality (if available, normally not for RNA)
	$run_qual = "n/a";
	$res_rq = $db->executeQuery("SELECT r.quality FROM processed_sample ps, sequencing_run r WHERE ps.id='$ps_id' AND r.id=ps.sequencing_run_id");
	if (count($res_rq))
	{
		$run_qual = $res_rq[0]['quality'];
	}
	$output['quality_run'] = $run_qual;
	
	//determine experiment type
	$output['experiment_type'] = ($s_type=="RNA" ? 'rna_seq' : 'dna_seq');

	//determine QC terms
	$qc_names = array("QC:2000005", "QC:2000007", "QC:2000008", "QC:2000025", "QC:2000027", "QC:2000014");
	$res_qc = $db->executeQuery("SELECT n.qcml_id, n.name, nm.value FROM qc_terms n, processed_sample_qc nm WHERE nm.processed_sample_id='$ps_id' AND nm.qc_terms_id=n.id");
	foreach($res_qc as $row_qc)
	{
		if (in_array($row_qc['qcml_id'], $qc_names))
		{
			$output['qc'][] = $row_qc;
		}
	}
	
	return $output;	
}

//removes values that are not exported to JSON
function prepareForExport($array)
{
	unset($array['ngsd_processedsample_id']);
	unset($array['normal_id']);
	unset($array['quality_sample']);
	unset($array['quality_processed_sample']);
	unset($array['quality_run']);
	unset($array['experiment_type']);
	
	return $array;
}

//marks a sample as already uploaded
function markAsUploaded($sample1, $sample2=null, $files)
{
	//in done.txt file
	$name = $sample1['id_genetics'];
	if (isset($sample2)) $name .= "-".$sample2['id_genetics'];
	$date = get_timestamp(false);
	$files = array_map("basename", $files);
	$files = implode(" ", $files);
	file_put_contents($GLOBALS["datamover_done"], "\n{$name}\t{$date}\t{$files}", FILE_APPEND);
	
	//in NGSD
	$db = DB::getInstance("NGSD");
	$text = "\nUploaded $name to QBIC with QBIC-ID ".$sample1['id_qbic'];
	if (isset($sample2)) $name .= "/".$sample2['id_qbic'];
	$db->executeStmt("UPDATE processed_sample SET comment = CONCAT(COALESCE(comment,''), '$text') WHERE id='".$sample1['ngsd_processedsample_id']."'");
}


//returns the date/files, or 'false' if the sample was not uploaded yes
function alreadyUploaded($name)
{
	//make sure the file exists
	$dm_done = $GLOBALS["datamover_done"];
	if (!file_exists($dm_done)) touch($dm_done);
	
	//load first column of file
	$tmp = file($dm_done);
	foreach($tmp as $line)
	{
		$line = trim($line);
		if ($line=="") continue;
		
		list($file, $date, $files) = explode("\t", $line."\t\t");
		if ($file==$name)
		{
			return array($date, $files);
		}
	}
	
	return false;
}

//stores meta data 
function storeMetaData($path, $name, $data)
{
	//print json_encode($data, JSON_PRETTY_PRINT);
	file_put_contents($path."/".$name.".metadata", json_encode($data));
}

//copies files to a folder (or touches output file in debug mode)
function copyFiles($files, $to_folder, $upload)
{
	global $parser;

	foreach($files as $file)
	{
		$outfile = $to_folder."/".basename($file);

		if (!$upload)
		{
			if (!touch($outfile))
			{
				trigger_error("Could not touch '$outfile'!", E_USER_ERROR);
			}
		}
		else
		{
			$parser->copyFile($file, $outfile);
		
			//determine SHA256
			list($stdout) = $parser->exec("sha256sum", $outfile, true);
			list($checksum) = explode(" ", $stdout[0]);
			print "##sha256sum ".basename($outfile)." {$checksum}\n";
		}
	}
}

//prints a TSV output line
function printTSV($output, $state, $state_comment)
{	
	$header = $GLOBALS['tsv_header'];
	$expected = count($header);
	
	//fill missing fields with 'n/a'
	while(count($output)<$expected-2)
	{
		$output[] = "n/a";
	}
	
	//append state
	$output[] = $state;
	$output[] = $state_comment;
	
	//error when too many cols are given
	if (count($output)>$expected)
	{
		for($i=0; $i<count($output); ++$i)
		{
			print ($i<count($header) ? $header[$i] : "[extra]").": ".$output[$i]."\n";
		}
		
		trigger_error("Too many columns in output!", E_USER_ERROR);
	}
	
	//print
	print implode("\t", $output)."\n";
}


//print TSV header
$GLOBALS['tsv_header'] = array("sample_name", "qbic_name", "fo_name", "quality(sample,ps,run)", "ps_name(s)", "ps_project", "ps_system", "type", "quality_normal(sample,ps,run)", "state", "state_comment");
print "#".implode("\t", $GLOBALS['tsv_header'])."\n";

//init
$db = DB::getInstance("NGSD");
$conditions = array();
if (empty($samples))
{
	$conditions[] = "name_external LIKE '%Q%'";
	$conditions[] = "name_external LIKE '%FO%'";
}
else
{
	foreach($samples as $s)
	{
		list($name, $ps_num) = explode("_", $s);
		$conditions[] = "name='$name'";
	}
}

$res = $db->executeQuery("SELECT id, name, name_external, quality, tumor FROM sample WHERE ".implode(" || ", $conditions)." ORDER BY name");
foreach($res as $row)
{
	list ($s_id, $s_name, $s_name_ex, $s_qual) = array_values($row);
	//check if we have a QBIC name
	$names = getExternalNames($s_name, $s_name_ex);
	if (is_null($names)) continue;
	list($qbic_name, $fo_name) = $names;
	
	
	$output = array();
	$output[] = $s_name;
	$output[] = $qbic_name;
	$output[] = $fo_name;
	$output[] = $s_qual;
	
	//check all processed samples
	$output_sample = $output;
	$res2 = $db->executeQuery("SELECT ps.id, p.name FROM processed_sample ps, project p WHERE sample_id='$s_id' AND ps.project_id=p.id AND ps.id NOT IN (SELECT processed_sample_id FROM merged_processed_samples)");
	foreach($res2 as $row2)
	{
		//truncate output to sample info
		$output = $output_sample;
		
		$ps_id = $row2['id'];
		$pro = $row2['name'];
		$sample1 = getSampleInfo($row2['id']);
		$ps_name = $sample1['id_genetics'];
		$sys_name = $sample1['processing_system'];
		
		
		$output[] = $ps_name;
		$output[] = $pro;
		$output[] = $sys_name;
		$output[3] .= ",".$sample1['quality_processed_sample'].",".$sample1['quality_run'];
		$output[] = $sample1['experiment_type'];

		//skip wrong projects
		if ($project!="" && $pro!=$project)
		{
			//printTSV($output, "RESTRICTED" ,"wrong project");
			continue;
		}

		//skip wrong samples
		if (!empty($samples) && !in_array($ps_name,$samples))
		{
			//printTSV($output, "RESTRICTED" ,"wrong processed sample");
			continue;
		}

		//skip bad processed samples
		if ($sample1['quality_processed_sample']=="bad" && !$ignore_quality)
		{
			printTSV($output, "SKIPPED" ,"bad processed sample quality");
			continue;
		}
		
		//skip bad runs
		if ($sample1['quality_run']=="bad" && !$ignore_quality)
		{
			printTSV($output, "SKIPPED" ,"bad run quality");
			continue;
		}

		//get folder and files
		$skipped = false;
		$missing = false;

		//skip samples where data is not at default location
		$info = get_processed_sample_info($db,$ps_name);
		$project_folder = $info["project_folder"];
		$data_folder = $info["ps_folder"];
		if (!$skipped && !file_exists($data_folder))
		{
			printTSV($output, "ERROR" ,"folder does not exist");
			$missing = true;
			continue;
		}
		
		//determine files to transfer
		$files = array();
		$paths  = glob($data_folder.$s_name."*.*");
		if(is_dir($data_folder."+original")) $paths = glob($data_folder."+original/".$ps_name."*.*");	//for backward compatibility
		
		//get FASTQ/VCF
		$tumor_normal_pair = $info["is_tumor"] == 1 && !empty($info["normal_id"]);
		$is_single_sample_dna = !$tumor_normal_pair && $sample1["experiment_type"] != "rna_seq";
		$fastqs_present = false;
		foreach($paths as $file)
		{
			if (!$tumor_normal_pair && ends_with($file, ".fastq.gz") &&
					!($sample1["experiment_type"] == "rna_seq" && $sample1["tumor"] == 1))
			{
				$files[] = $file;
				$fastqs_present= true;
			}
			if (ends_with($file, "_var_annotated.vcf.gz")) $files[] = $file;
		}

		// Special treatment for PacBio project (Currently works for unmapped CLR Reads)
		if($info['device_type']=="SequelII")
		{	
			$files[] = "{$data_folder}{$ps_name}.bam";
			$files[] = "{$data_folder}{$ps_name}.bam.pbi";
			$files[] = "{$data_folder}{$ps_name}.subreadset.xml";
		} 
	
		//generate FASTQs from BAM is not available
		if ($is_single_sample_dna && !$fastqs_present && !$info['device_type']=="SequelII")
		{
			$bam = "{$data_folder}/{$ps_name}.bam";
			if (file_exists($bam))
			{
				$tmp_folder = $parser->tempFolder();	
				$fq1 = $tmp_folder."/{$ps_name}_BamToFastq_R1_001.fastq.gz";
				$fq2 = $tmp_folder."/{$ps_name}_BamToFastq_R2_001.fastq.gz";
				if ($upload) //skip generating FASTQs in dry run
				{
					$parser->exec(get_path("ngs-bits")."BamToFastq", "-in {$bam} -out1 $fq1 -out2 $fq2", true);
				}
				$files[] = $fq1;
				$files[] = $fq2;
			}
		}
		
		//Special treatment for tumor-normal fastqs
		$merged_fastq_files = array();
		if($tumor_normal_pair) //merge fastqs in case of normal-tumor pairs 
		{
			if($row['tumor'] == 1)
			{
				exec2("cat {$data_folder}*R1*.fastq.gz > {$data_folder}{$qbic_name}_tumor.1.fastq.gz");
				$merged_fastq_files[] = "{$data_folder}{$qbic_name}_tumor.1.fastq.gz";
				exec2("cat {$data_folder}*R2*.fastq.gz > {$data_folder}{$qbic_name}_tumor.2.fastq.gz");
				$merged_fastq_files[] = "{$data_folder}{$qbic_name}_tumor.2.fastq.gz";
			}
			//parse related normal file if found
			$normal_id = $sample1["normal_id"];
			$normal_sample = getSampleInfo($normal_id);
			if($normal_id != "" && $normal_sample['quality_sample'] != "bad" && $normal_sample['quality_processed_sample'] != "bad"
				&& $normal_sample['quality_run'] != "bad")
			{
				$data_folder2 = $project_folder."/Sample_".$normal_sample['id_genetics']."/";
				exec2("cat {$data_folder2}*R1*.fastq.gz > {$data_folder2}{$qbic_name}_normal.1.fastq.gz");
				$merged_fastq_files[] = "{$data_folder2}{$qbic_name}_normal.1.fastq.gz";
				exec2("cat {$data_folder2}*R2*.fastq.gz > {$data_folder2}{$qbic_name}_normal.2.fastq.gz");
				$merged_fastq_files[] = "{$data_folder2}{$qbic_name}_normal.2.fastq.gz";
			}
		}
		elseif($sample1["experiment_type"] == "rna_seq" && $sample1["tumor"] == 1)
		{
			exec2("cat {$data_folder}*R1*.fastq.gz > {$data_folder}{$qbic_name}_tumor_rna.1.fastq.gz");
			$merged_fastq_files[] = "{$data_folder}{$qbic_name}_tumor_rna.1.fastq.gz";
			exec2("cat {$data_folder}*R2*.fastq.gz > {$data_folder}{$qbic_name}_tumor_rna.2.fastq.gz");
			$merged_fastq_files[] = "{$data_folder}{$qbic_name}_tumor_rna.2.fastq.gz";
		}
		
		
		//skip already uploaded
		$skipped = false;
		$uploaded = false;
		if(!$force_reupload)	$uploaded = alreadyUploaded($ps_name);
		if ($uploaded!==false)
		{
			printTSV($output, "SKIPPED" ,"already uploaded on ".$uploaded[0].": ".$uploaded[1]."; found files for upload: ".implode(", ",array_map("basename",$files)));
			$skipped = true;
		}
		
		//skip if no files were found
		if (!$skipped)
		{
			if (count($files)==0 && count($merged_fastq_files) == 0)
			{
				printTSV($output, "ERROR" ,"no files to transfer in data folder '$data_folder'");
				$skipped = true;
			}
		}

		//determine/create subfolder
		if (!$skipped)
		{
			$folder_name = $qbic_name."_".$ps_name;
			$tmpfolder = $GLOBALS["datamover_tmp"]."/".$folder_name;
			if (file_exists($tmpfolder)) exec2("rm -rf $tmpfolder");
			mkdir($tmpfolder);
			
			//copy files to folder
			copyFiles($files, $tmpfolder, $upload);

			foreach($merged_fastq_files as $file)
			{
				$parser->moveFile($file,$tmpfolder."/".basename($file));
			}
			
			//store meta data file
			storeMetaData($tmpfolder, $ps_name, array("type"=>$sample1["experiment_type"], "sample1"=>prepareForExport($sample1), "files"=>array_map("basename", $files)));
			exec2("chmod -R 777 $tmpfolder");
			
			//upload data
			if ($upload)
			{
				$parser->moveFile($tmpfolder, $GLOBALS["datamover_path"]."/".$folder_name);
				$files = array_merge($files, $merged_fastq_files);
				markAsUploaded($sample1, null, $files);
			}
			printTSV($output, $upload ? "UPLOADED" : "TO_UPLOAD" , implode(", ", $files));
		}
	
		//handle tumor-normal pairs
		$somatic_data_uploaded = false;
		$normal_id = $sample1["normal_id"];
		if ($normal_id!="")
		{
			$sample2 = getSampleInfo($normal_id);
			$ps_name2 = $sample2['id_genetics'];
			$output[4] .= "-".$ps_name2;
			$output[] = $sample2['quality_sample'].",".$sample2['quality_processed_sample'].",".$sample2['quality_run'];
			
			//check sample quality	
			if ($sample2['quality_sample']=="bad" && !$ignore_quality)
			{
				printTSV($output, "SKIPPED" , "bad sample quality");
				continue;
			}
			if ($sample2['quality_processed_sample']=="bad" && !$ignore_quality)
			{ 
				printTSV($output, "SKIPPED" , "bad processed sample quality");
				continue;
			}
			if ($sample2['quality_run']=="bad" && !$ignore_quality)
			{ 
				printTSV($output, "SKIPPED" , "bad run quality");
				continue;
			}
				
			//skip tumor-normal pairs without unique data folder
			$folders = glob("{".$project_folder."Somatic_{$ps_name}-{$ps_name2}*}",GLOB_BRACE);
			if (count($folders)==0)
			{
				printTSV($output, "ERROR" , "somatic folder does not exist");
				continue;
			}
			else if (count($folders)>1)
			{
				printTSV($output, "ERROR" , "somatic folder exists several times: ".implode(" ", $folders));
				continue;
			}
			$data_folder = $folders[0];
			
			//determine files to transfer
			$patterns =  array();
			$patterns[] = "{$ps_name}-{$ps_name2}*_var_annotated.vcf.gz";
			$patterns[] = "{$ps_name}-{$ps_name2}*.GSvar";
			$patterns[] = "{$ps_name2}*_var_annotated.vcf.gz";
			$patterns[] = "{$ps_name2}*.GSvar";
			$patterns[] = "*_counts.tsv,*_var_fusions.tsv";
			$files = glob("{$data_folder}/{".implode(",", $patterns)."}", GLOB_BRACE);
						
			//skip already uploaded
			$uploaded = false;
			if(!$force_reupload)	$uploaded = alreadyUploaded($ps_name."-".$ps_name2);
			if ($uploaded!==false)
			{
				printTSV($output, "SKIPPED" ,"already uploaded on ".$uploaded[0].": ".$uploaded[1]."; found files for upload: ".implode(", ",array_map('basename',$files)).".");
				continue;
			}

			if (count($files)==0)
			{
				printTSV($output, "ERROR" , "no files to transfer in data folder '$data_folder'");
				continue;
			}
			
			//determine/create subfolder
			$folder_name = $qbic_name."_".$ps_name."-".$ps_name2;
			$tmpfolder = $GLOBALS["datamover_tmp"]."/".$folder_name;
			if (file_exists($tmpfolder)) exec2("rm -rf $tmpfolder");
			mkdir($tmpfolder);
			
			//copy files to transfer folder
			copyFiles($files, $tmpfolder, $upload);
						
			//store meta data
			storeMetaData($tmpfolder, $ps_name."-".$ps_name2, array("type"=>"dna_seq_somatic", "sample1"=>prepareForExport($sample1), "sample2"=>prepareForExport($sample2), "files"=>array_map("basename", $files)));
			exec2("chmod -R 777 $tmpfolder");
			
			//upload data
			if ($upload)
			{
				$parser->moveFile($tmpfolder, $GLOBALS["datamover_path"]."/".$folder_name);
				markAsUploaded($sample1, $sample2, $files);
			}
			printTSV($output, $upload ? "UPLOADED" : "TO_UPLOAD" , implode(" ", $files));
			
			$somatic_data_uploaded = true;
		}
		
		//Upload somatic data in special format for MTB (Molecular Tumor Board)
		if ($somatic_data_uploaded && file_exists("{$data_folder}/QBIC_files/"))
		{
			$output[4] = "{$ps_name}-{$ps_name2} (MTB)";
			
			//determine files to zip
			$files = glob("{$data_folder}/QBIC_files/*.tsv");
			
			//determine/create subfolder
			$folder_name = "{$qbic_name}_{$ps_name}-{$ps_name2}-MTB";
			$tmpfolder = $GLOBALS["datamover_tmp"]."/".$folder_name;
			if (file_exists($tmpfolder)) exec2("rm -rf $tmpfolder");
			mkdir($tmpfolder);
			
			//copy and rename files
			foreach($files as $file)
			{
				$outfile = basename($file);
				$outfile =  $qbic_name.substr($outfile, 4);
				$parser->copyFile($file, $tmpfolder."/".$outfile);
			}
			
			//zip files
			$zip = "{$tmpfolder}/{$qbic_name}_{$ps_name}-{$ps_name2}.zip";
			exec2("cd {$tmpfolder} && zip {$zip} *.tsv");
			
			//upload data
			if ($upload)
			{
				$parser->moveFile($zip, $GLOBALS["datamover_path"]."/".basename($zip));
			}
			printTSV($output, $upload ? "UPLOADED" : "TO_UPLOAD" , implode(" ", $files));
		}
	}
}


?>
