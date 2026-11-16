<?php 

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("db_cleanup_seconary_analyses", "Clean up entries in 'secondary_analysis' table of NGSD.");
$parser->addInfile("samples", "Samples to process. Create with 'NGSDExportSamples -add_path SAMPLE_FOLDER'.", false);
extract($parser->parse($argv));


function log_files_only($folder)
{
    foreach (new DirectoryIterator($folder) as $file)
    {
        if ($file->isDot()) continue;

        if (!$file->isFile() || strtolower($file->getExtension()) !== 'log')
        {
            return false;
        }
    }

    return true;
}

//init
$db = DB::getInstance("NGSD");

//get processed sample names for which the sample folder still exists
$ps_folder_exists = [];
$cache = $samples.".cache";
if (!file_exists($cache) || filemtime($cache)<filemtime($samples))
{
	list($stdout) = exec2("TsvSlice -in $samples -cols name,path");
	print "##PS in NGSD: ".count($stdout)."\n";
	foreach($stdout as $line)
	{
		$line = nl_trim($line);
		if ($line=="" || $line[0]=="#") continue;
		
		list($ps, $folder) = explode("\t", $line);
		
		//check if folder exists - we cannot check for GSvar because e.g. somatic tumor samples contain no GSvar file
		if (file_exists($folder))
		{
			$ps_folder_exists[$ps] = true;
		}
	}
	print "##PS in NGSD with existing folder: ".count($ps_folder_exists)."\n";
	
	//store cache file
	file_put_contents($cache, implode("\n", array_keys($ps_folder_exists)));
}
else
{
	foreach(file($cache) as $ps)
	{
		$ps = trim($ps);
		if ($ps=="") continue;
		
		$ps_folder_exists[$ps] = true;
	}
	print "##PS from cache: ".count($ps_folder_exists)."\n";
}

//process secondary analysis entries
print "#id\tgsvar\tstatus\n";
foreach($db->executeQuery("SELECT * FROM secondary_analysis") as $row)
{
	$id = $row['id'];
	$type = $row['type'];
	$gsvar = trim($row['gsvar_file']);
	$folder = dirname($gsvar);
	
	//folder does not exist > delete
	if (!file_exists($folder))
	{
		print "$id\t$gsvar\tDELETE (folder missing)\n";
		continue;
	}
	
	//folder exists but contains only log files
	if (log_files_only($folder))
	{
		print "$id\t$gsvar\tDELETE (folder contains only log files)\n";
		continue;
	}	
	
	//GSvar does not exist > ???
	if (!file_exists($gsvar))
	{
		print "$id\t$gsvar\tDELETE? (GSvar missing)\n";
		continue;
	}
	
	//GSvar exists > check if PS folders exist
	$foldername = basename($folder);
	$foldername = preg_replace('/(_[0-9]{2})-/', '$1_', $foldername); //somatic T/N paths contain '-' as separator, multi/trio paths contain '_' > normalize to '_'
	$foldername = strtr($foldername, ["Somatic_"=>"", "Trio_"=>"", "Multi_"=>"", "Sample_"=>""]);

	$invalid_ps = [];
	$parts = explode("_", $foldername);
	if (count($parts) % 2 !== 0)
	{
		print "##Invalid part count: ".implode(", ", $parts)."\n";
		continue;
	}
	for($i=0; $i<count($parts); $i+=2)
	{
		$ps = $parts[$i]."_".$parts[$i+1];
		if (!isset($ps_folder_exists[$ps]))
		{
			$invalid_ps[] = $ps;
		}
	}
	if (count($invalid_ps)>0)
	{
		print "$id\t$gsvar\tDELETE? (invalid PS: ".implode(", ", $invalid_ps).")\n";
		continue;
	}
	fflush(STDOUT);
}

?>