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

//wait sync is done by another process 
$pid_file = sys_get_temp_dir()."/megSAP_data_setup_{$build}.txt";
if (file_exists($pid_file))
{
	$pid_old = trim(file_get_contents($pid_file));
	$iter = 0;
	while (posix_getpgid($pid_old)!==FALSE)
	{
		++$iter;
		if ($iter>30) break; //wait for 30 min max
		print "Process with PID {$pid_old} is already syncing. Waiting one minute...\n";
		sleep(60);
	}
}
if (!file_put_contents($pid_file, getmypid(), LOCK_EX))
{
	trigger_error("Could not create PID file {$pid_file}", E_USER_ERROR);
}
$pid_mod = substr(sprintf('%o', fileperms($pid_file)), -4);
if ($pid_mod!="0777")
{
	chmod($pid_file, 0777);
}

//copy reference genome and index files
print "Copying genome files...\n";
list($files) = exec2("ls {$genome_folder}/{$build}.*");
foreach($files as $file)
{
	$base = basename($file);
	print "  rsync-ing genome file '$base'.\n";
	$existed_before = file_exists($local_folder.$base);
	list($stdout) = exec2("rsync --archive --acls --no-perms --no-group --no-owner --chmod=ugo=rwX ".$genome_folder.$base." ".$local_folder.$base);
	foreach($stdout as $line)
	{
		$line = trim($line);
		if ($line=="") continue;
		print "    $line\n";
	}
	if (!$existed_before)
	{
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
	$info = "/homo_sapiens/95_GRCh37/info.txt";
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
			print "Annotation infos in '$info' differ. Deleting old data and performing update!\n";
			exec2("rm -rf {$local_annotation_folder}/*");
		}
	}
	
	//check that at least one rsync finished
	$finished = "/rsync_done.txt";
	if (!file_exists("{$local_annotation_folder}/{$finished}"))
	{
		print "No rsync finished. Performing update!\n";
		$update = true;
	}
		
	//perform update
	if ($update)
	{
		print "rsync-ing annotation data...\n";
		
		list($stdout) = exec2("rsync --archive --omit-dir-times --acls --no-perms --no-group --no-owner --chmod=ugo=rwX {$annotation_folder} {$local_annotation_folder}");
		foreach($stdout as $line)
		{
			print trim($line)."\n";
		}
		
		touch("{$local_annotation_folder}/{$finished}");
	}
}

?>