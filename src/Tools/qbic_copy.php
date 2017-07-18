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
$parser->addFlag("force_sap", "Makes the SAP ID mandatory (and the QBIC ID optional).", true, "");
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
	die;
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

//cache samples folders list to speed up RNA samples search
print "##\n";
print "##Caching sample folders...\n";
list($sample_folders) = exec2("find ".get_path("project_folder")." -maxdepth 3 -type d -name 'Sample_*'");

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

//determine SAP identifier from GenLab
function getSapId($name)
{	
	$db = DB::getInstance("GL8");
	$res = $db->executeQuery("SELECT identnr FROM v_ngs_sap WHERE labornummer LIKE '$name'");
	
	//no hit
	if (count($res)==0) return null;
	
	//more than one hit => error
	if (count($res)>1)
	{
		print_r($res);
		trigger_error("Found more than one sample '$name' in GenLab!", E_USER_ERROR);
	}
	
	//one hit
	$sap_id = trim($res[0]['identnr']);
	if ($sap_id=="") return null;
	
	return $sap_id;
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
	if (is_null($names)) //this can only be null when 'force_sap' is used
	{
		$output['id_qbic'] = "";
	}
	else
	{
		$output['id_qbic'] = $names[0];
	}
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
	$db->executeStmt("UPDATE processed_sample SET comment = CONCAT(comment, '$text') WHERE id='".$sample1['ngsd_processedsample_id']."'");
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
			copy2($file, $outfile);
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

$res = $db->executeQuery("SELECT id, name, name_external, quality FROM sample WHERE ".implode(" || ", $conditions)." ORDER BY name");
foreach($res as $row)
{
	list ($s_id, $s_name, $s_name_ex, $s_qual) = array_values($row);
	
	//check if we have a QBIC name
	$names = getExternalNames($s_name, $s_name_ex);
	if (is_null($names) && !$force_sap) continue;
	list($qbic_name, $fo_name) = $names;
	
	$output = array();
	$output[] = $s_name;
	$output[] = $qbic_name;
	$output[] = $fo_name;
	$output[] = $s_qual;
	
	//skip bad samples
	if ($s_qual=="bad")
	{
		printTSV($output, "SKIPPED", "bad sample quality");
		continue;
	}
	
	//check all processed samples
	$output_sample = $output;
	$res2 = $db->executeQuery("SELECT ps.id, p.name FROM processed_sample ps, project p WHERE sample_id='$s_id' AND ps.project_id=p.id");
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
		if ($sample1['quality_processed_sample']=="bad")
		{
			printTSV($output, "SKIPPED" ,"bad processed sample quality");
			continue;
		}
		
		//skip bad runs
		if ($sample1['quality_run']=="bad")
		{
			printTSV($output, "SKIPPED" ,"bad run quality");
			continue;
		}
		
		//get folder and files
		$skipped = false;
		$missing = false;
		if($sample1['experiment_type']=="dna_seq")
		{
			//skip samples where data is not at default location
			$res3 = $db->executeQuery("SELECT p.type, p.name FROM processed_sample ps, project p WHERE ps.id='$ps_id' AND ps.project_id=p.id");
			$project_folder = get_path("project_folder")."/".$res3[0]['type']."/".$res3[0]['name']."/";
			$data_folder = $project_folder."Sample_".$ps_name."/";
			if (!$skipped && !file_exists($data_folder))
			{
				printTSV($output, "ERROR" ,"folder does not exist");
				$missing = true;
				continue;
			}
		}
		else //rna_seq
		{
			if (!$skipped)
			{
				$hits = array();
				foreach($sample_folders as $f)
				{
					if (contains($f, "Sample_".$ps_name))
					{
						$hits[] = $f;
					}
				}
				
				if (count($hits)>1)
				{
					printTSV($output, "ERROR" ,"several folders found: ".implode(", ", $hits));
					continue;
				}
				if (count($hits)==0)
				{
					printTSV($output, "ERROR" ,"folder does not exist");
					continue;
				}
				$data_folder = $hits[0]."/";
			}
		}
		
		//determine files to transfer
		$files = array();
		$paths  = glob($data_folder.$s_name."*.*");
		if(is_dir($data_folder."+original")) $paths = glob($data_folder."+original/".$ps_name."*.*");	//for backward compatibility
		foreach($paths as $file)
		{
			if (ends_with($file, ".fastq.gz")) $files[] = $file;
			if (ends_with($file, "_var_annotated.vcf.gz")) $files[] = $file;
			if (ends_with($file, "_counts_raw.tsv")) $files[] = $file;
			if (ends_with($file, "_counts_fpkm.tsv")) $files[] = $file;
			if (ends_with($file, ".bam")) $files[] = $file;
			if (ends_with($file, ".bam.bai")) $files[] = $file;
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
			if (count($files)==0)
			{
				printTSV($output, "ERROR" ,"no files to transfer in data folder '$data_folder'");
				$skipped = true;
			}
		}
		
		//check SAP identifier
		if (!$skipped)
		{	
			$sample1['id_sap'] = getSapId($s_name);
			if ($force_sap && is_null($sample1['id_sap']))
			{
				printTSV($output, "ERROR" ,"No SAP identifier found for '$ps_id'!");
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
			
			//store meta data file
			storeMetaData($tmpfolder, $ps_name, array("type"=>$sample1["experiment_type"], "sample1"=>prepareForExport($sample1), "files"=>array_map("basename", $files)));
			exec2("chmod -R 777 $tmpfolder");
			
			//upload data
			if ($upload)
			{
				rename($tmpfolder, $GLOBALS["datamover_path"]."/".$folder_name);
				markAsUploaded($sample1, null, $files);
			}
			printTSV($output, $upload ? "UPLOADED" : "TO_UPLOAD" , implode(", ", $files));
		}
	
		//handle tumor-normal pairs
		$normal_id = $sample1["normal_id"];
		if ($normal_id!="")
		{
			$sample2 = getSampleInfo($normal_id);
			$ps_name2 = $sample2['id_genetics'];
			$output[4] .= "-".$ps_name2;
			$output[] = $sample2['quality_sample'].",".$sample2['quality_processed_sample'].",".$sample2['quality_run'];
			
			//check sample quality	
			if ($sample2['quality_sample']=="bad")
			{
				printTSV($output, "SKIPPED" , "bad sample quality");
				continue;
			}
			if ($sample2['quality_processed_sample']=="bad")
			{ 
				printTSV($output, "SKIPPED" , "bad processed sample quality");
				continue;
			}
			if ($sample2['quality_run']=="bad")
			{ 
				printTSV($output, "SKIPPED" , "bad run quality");
				continue;
			}
				
			//skip tumor-normal pairs without unique data folder
			$folders = glob("{".$project_folder."Somatic_{$ps_name}-{$ps_name2}*,".$project_folder."Sample_{$ps_name}-{$ps_name2}*}",GLOB_BRACE);	//backward compatibility
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
			
			$old_pipeline = false;
			//get pipeline version
			if(strpos($data_folder,$project_folder."Sample_")!==FALSE)	$old_pipeline = true;

			//determine files to transfer
			$pattern = null;
			if($old_pipeline)
			{
				$pattern = "{".$data_folder."/{$ps_name}-{$ps_name2}*_vc_strelka.vcf,".$data_folder."/{$ps_name}-{$ps_name2}*.GSvar}"; #,".$data_folder."/{$ps_name2}.GSvar,".$data_folder."/{$ps_name2}_var.vcf
			}
			else
			{
				$pattern = "{$data_folder}/{{$ps_name}-{$ps_name2}*_var_annotated.vcf.gz,{$ps_name}-{$ps_name2}*.GSvar,{$ps_name2}*_var_annotated.vcf.gz,{$ps_name2}*.GSvar}";
			}
			
			$files = glob($pattern, GLOB_BRACE);
						
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
			
			//get SAP identifier
			list($s_name2) = explode("_", $ps_name2);
			$sample2['id_sap'] = getSapId($s_name);			
			if ($force_sap && is_null($sample2['id_sap']))
			{
				printTSV($output, "ERROR" ,"No SAP identifier found for '$ps_name2'!");
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
				rename($tmpfolder, $GLOBALS["datamover_path"]."/".$folder_name);
				markAsUploaded($sample1, $sample2, $files);
			}
			printTSV($output, $upload ? "UPLOADED" : "TO_UPLOAD" , implode(" ", $files));
		}
	}
}

?>
