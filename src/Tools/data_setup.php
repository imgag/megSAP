<?php
/** 
	@page data_setup
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("data_setup", "Creates a local copy of static NGS data that is heavily used during data analysis and should not be accessed via the network.");
$parser->addString("build", "Genome build.", false);
extract($parser->parse($argv));

$local_folder = get_path("local_data");
$genome_folder = get_path("data_folder")."/genomes/";

//skip if they are the same
if (realpath($local_folder)==realpath($genome_folder))
{
	print "Skipping, because local and global genome folders are the same:\n";
	print "local: ".realpath($local_folder)."\n";
	print "global: ".realpath($genome_folder)."\n";
	return;
}

//make local folder
if (!file_exists($local_folder))
{
	if (!mkdir($local_folder))
	{
		trigger_error("Could not create local data folder '{$local_folder}'!", E_USER_ERROR);
	}
	if (!chmod($local_folder, 0777))
	{
		trigger_error("Could not change privileges of local data folder '{$local_folder}'!", E_USER_ERROR);
	}
}

//check md5 sum of genome
if (file_exists("{$local_folder}/{$build}.md5"))
{
	exec("diff {$local_folder}/{$build}.md5 {$genome_folder}/{$build}.md5", $output, $code);
	if ($code==0)
	{
		print "MD5 checksums of genome '$build' match.\n";
	}
	else
	{
		print "MD5 checksums of genome '$build' differ. Deleting old data!\n";
		exec2("rm -rf {$local_folder}/{$build}.*");
	}
}

//copy reference genome and index files
list($files) = exec2("ls {$genome_folder}/{$build}.*");
foreach($files as $file)
{
	$base = basename($file);
	if (!file_exists($local_folder.$base))
	{
		print "Copying genome file '$base' to $local_folder\n";
		$parser->copyFile($genome_folder.$base, $local_folder.$base);
		if (!chmod($local_folder.$base, 0777))
		{
			trigger_error("Could not change privileges of local data folder '{$local_folder}{$base}'!", E_USER_ERROR);
		}
	}
}

?>