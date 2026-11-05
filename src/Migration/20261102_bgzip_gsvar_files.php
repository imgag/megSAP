<?php 

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$db = DB::getInstance("NGSD");

$parser = new ToolBase("20261102_bgzip_gsvar_files", "Use bgzip to compress GSvar files.");
$parser->addInfile("samples", "Samples to process. Create with 'NGSDExportSamples -add_path SAMPLE_FOLDER'.", true, ["prepare", "deploy"], "prepare");  
$parser->addEnum("mode", "Execution mode.", true, ["prepare", "deploy"], "prepare");  
extract($parser->parse($argv));

function exists_and_not_empty($filename)
{
	return file_exists($filename) && filesize($filename)>0;
}

//single-sample analyses
print "#ps\tfolder\tstatus\n";
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
	$gsvar_exists = exists_and_not_empty($gsvar);
	if (!$gsvar_exists)
	{
		++$c_no_gsvar;
		continue;
	}
	
	if ($mode=="prepare")
	{
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
		if ($gsvar_exists && $order_ok && $gz_ok)
		{
			print $ps."\t".$folder."\tGZ already exists\n";
		}
		else
		{
			list($stdout2, $stderr2, $exit_code2) = exec2("bgzip -@ 8 -c $gsvar > $gsvar_gz", false);
			if ($exit_code2==0)
			{
				print $ps."\t".$folder."\tGZ created\n";
			}
			else
			{
				print $ps."\t".$folder."\tGZ creation failed\n";
			}
		}
	}
	else
	{
		//TODO
	}
	fflush(STDOUT);
}
print "##PS without folder: $c_no_folder\n";
print "##PS without GSvar file: $c_no_gsvar\n";

//secondary analyses (from NGSD)
//TODO

//TODO: clean up left-over TMP files


?>