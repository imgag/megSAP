<?php
/** 
	@page data_setup
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("data_setup", "Creates a local copy of the reference genome and annotation data. They are heavily used during data analysis and should not be accessed via the network!.");
$parser->addString("build", "Genome build.", false);
extract($parser->parse($argv));

$local_folder = get_path("local_data");
$genome_folder = get_path("data_folder")."/genomes/";

######################### reference genome #########################
print "### Reference genome ###\n";
print "from: {$genome_folder}\n";
print "to  : {$local_folder}\n";
print "\n";

//skip if they are the same
if (realpath($local_folder)==realpath($genome_folder))
{
	print "Skipping, because local and remote genome folders are the same!\n";
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


//remove outdated genome FASTA files
if (file_exists("{$local_folder}/{$build}.fa.md5"))
{
	exec("diff {$local_folder}/{$build}.fa.md5 {$genome_folder}/{$build}.fa.md5", $output, $code);
	if ($code==0)
	{
		print "MD5 checksums of genome '{$build}.fa.md5' match.\n";
	}
	else
	{
		print "MD5 checksums of genome '{$build}.fa.md5' differ. Deleting old data!\n";
		exec2("rm -rf {$local_folder}/{$build}.*");
	}
}

//copy reference genome and index files
print "Copying genome files (if missing)...\n";
list($files) = exec2("ls {$genome_folder}/{$build}.*");
foreach($files as $file)
{
	$base = basename($file);
	if (!file_exists($local_folder.$base))
	{
		print "Copying genome file '$base'.\n";
		$parser->copyFile($genome_folder.$base, $local_folder.$base);
		if (!chmod($local_folder.$base, 0777))
		{
			trigger_error("Could not change privileges of local data folder '{$local_folder}{$base}'!", E_USER_ERROR);
		}
	}
}

######################### VEP annotation data #########################
if ($build=="GRCh37")
{
	$annotation_folder = get_path("vep_data")."/cache/";
	$local_annotation_folder = get_path("local_data")."/".basename(get_path("vep_data"))."/";

	print "\n";
	print "### VEP annotation data ###\n";
	print "from: {$annotation_folder}\n";
	print "to  : {$local_annotation_folder}\n";
	print "\n";

	//create local annotation folder if missing
	if (!file_exists($local_annotation_folder))
	{
		if (!mkdir($local_annotation_folder))
		{
			trigger_error("Could not create local data annotation folder '{$local_annotation_folder}'!", E_USER_ERROR);
		}
		if (!chmod($local_annotation_folder, 0777))
		{
			trigger_error("Could not change privileges of local annotation data folder '{$local_annotation_folder}'!", E_USER_ERROR);
		}
	}
	
	//remove outdated annotation data
	$update = true;
	$info = "/homo_sapiens/94_GRCh37/info.txt";
	if (file_exists("{$local_annotation_folder}/{$info}"))
	{
		exec("diff {$local_annotation_folder}/{$info} {$annotation_folder}/{$info}", $output, $code);
		if ($code==0)
		{
			print "Annotation infos in '$info' match.\n";
			$update = false;
		}
		else
		{
			print "Annotation infos in '$info' differ. Deleting old data!\n";
			exec2("rm -rf {$local_annotation_folder}/*");
		}
	}
	
	if ($update)
	{
		print "rsync-ing annotation data...\n";
		list($stdout) = exec2("rsync -a {$annotation_folder} {$local_annotation_folder}");
		foreach($stdout as $line)
		{
			print trim($line)."\n";
		}
	}
}

?>