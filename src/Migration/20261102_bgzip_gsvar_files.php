<?php 

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("20261102_bgzip_gsvar_files", "Use bgzip to compress GSvar files.");
$parser->addInfile("samples", "Samples to process. Create with 'NGSDExportSamples -add_path SAMPLE_FOLDER'.", false);
$parser->addFlag("delete_unzipped", "Delete unzipped GSvar file after creating bgzipped version.");
extract($parser->parse($argv));

$db = DB::getInstance("NGSD");

function exists_and_not_empty($filename)
{
	return file_exists($filename) && filesize($filename)>0;
}

print "#ps\tfolder\tstatus\n";

//single-sample analyses
print "### processing single sample analyses ###\n";
$c_no_folder = 0;
$c_no_gsvar = 0;
list($stdout) = exec2("TsvSlice -in $samples -cols name,path");
print "##PS overall in NGSD: ".count($stdout)."\n";
foreach($stdout as $line)
{
	$line = nl_trim($line);
	if ($line=="" || $line[0]=="#") continue;
	list($ps, $folder) = explode("\t", $line);
	
	//check if folder exists
	if (!file_exists($folder))
	{
		++$c_no_folder;
		continue;
	}
	
	//check if GSvar exists
	$gsvar = $folder."/".$ps.".GSvar";
	if (!exists_and_not_empty($gsvar))
	{
		++$c_no_gsvar;
		continue;
	}
	
	$gsvar_gz = $folder."/".$ps.".GSvar.gz";
	$gsvar_gz_exists = exists_and_not_empty($gsvar_gz);
	$order_ok = false;
	$gz_ok = false;
	if ($gsvar_gz_exists)
	{
		$order_ok = filemtime($gsvar)<filemtime($gsvar_gz);
		list($stdout2, $stderr2, $exit_code2) = exec2("gzip -t $gsvar_gz", false);
		$gz_ok = $exit_code2==0;
	}
	if ($order_ok && $gz_ok)
	{
		print $ps."\t".$folder."\tGZ already exists\n";
		if ($delete_unzipped) exec2("rm -f $gsvar");
	}
	else
	{
		list($stdout2, $stderr2, $exit_code2) = exec2("bgzip -@ 8 -c $gsvar > $gsvar_gz", false);
		if ($exit_code2==0)
		{
			print $ps."\t".$folder."\tGZ created\n";
			if ($delete_unzipped) exec2("rm -f $gsvar");
		}
		else
		{
			print $ps."\t".$folder."\tGZ creation failed\n";
		}
	}
	fflush(STDOUT);
}
print "##PS without folder: $c_no_folder\n";
print "##PS without GSvar file: $c_no_gsvar\n";

//secondary analyses (from NGSD)
$res = $db->executeQuery("SELECT * FROM secondary_analysis");
foreach($res as $row)
{
	$id = $row['id'];
	$type = $row['type'];
	$gsvar = $row['gsvar_file'];
	$folder = dirname($gsvar);
	
	//code to fix paths of moved secondary analyses (commented out because it was needed only once at IMGAG)
	/*
	$parts = array_slice(explode("/", $gsvar), 0, 6);
	$project_folder = implode("/", $parts)."/";
	if (!file_exists($project_folder))
	{
		$project_name = end($parts);
		$to = get_project_folder($db, $project_name);
		print "$project_folder => $to ($id)\n";
		if ($project_folder!=$to)
		{
			$gsvar2 = strtr($gsvar, [$project_folder=>$to]);
			$id2 = $db->getValue("SELECT id FROM secondary_analysis WHERE gsvar_file='$gsvar2'", -1);
			if ($id2==-1)
			{
				$db->executeStmt("UPDATE secondary_analysis SET gsvar_file='$gsvar2' WHERE id=$id");
				print "  UPDATED\n";
			}
			else
			{
				$db->executeStmt("DELETE FROM secondary_analysis WHERE id=$id");
				print "  DELETED (exists as $gsvar2 with ID $id2)\n";
			}
		}
	}
	*/
	
	//already as GZ in NGSD => nothing to do
	if (ends_with($gsvar, ".GSvar.gz")) continue;
	
	//check if GSvar exists
	if (!file_exists($gsvar))
	{
		if (!file_exists($folder))
		{
			print "secondary:$type\t$folder\tFOLDER missing\n";
		}
		else
		{
			print "secondary:$type\t$folder\tGSvar missing\n";
		}
		continue;
	}
	
	$gsvar_gz = $gsvar.".gz";
	$gsvar_gz_exists = exists_and_not_empty($gsvar_gz);
	$order_ok = false;
	$gz_ok = false;
	if ($gsvar_gz_exists)
	{
		$order_ok = filemtime($gsvar)<filemtime($gsvar_gz);
		list($stdout2, $stderr2, $exit_code2) = exec2("gzip -t $gsvar_gz", false);
		$gz_ok = $exit_code2==0;
	}
	if ($order_ok && $gz_ok)
	{
		print "secondary:$type\t".$folder."\tGZ already exists\n";
		if ($delete_unzipped)
		{
			//TODO update NGSD
			exec2("rm -f $gsvar");
		}
	}
	else
	{
		list($stdout2, $stderr2, $exit_code2) = exec2("bgzip -@ 8 -c $gsvar > $gsvar_gz", false);
		if ($exit_code2==0)
		{
			print "secondary:$type\t".$folder."\tGZ created\n";
			if ($delete_unzipped)
			{
				//TODO update NGSD
				exec2("rm -f $gsvar");
			}
		}
		else
		{
			print "secondary:$type\t".$folder."\tGZ creation failed\n";
		}
	}
	fflush(STDOUT);
}

?>