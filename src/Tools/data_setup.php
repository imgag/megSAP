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


//init
$data_folder = get_path("data_folder");
$local_data = get_path("local_data");
$rsync = "rsync --archive --omit-dir-times --acls --no-perms --no-group --no-owner --chmod=ugo=rwX";

######################### reference genome #########################
$genome_folder = "{$data_folder}/genomes/";
print "### Reference genome ###\n";
print "from: {$genome_folder}\n";
print "to  : {$local_data}\n";
print "\n";

//make local folder
if (!file_exists($local_data))
{
	if (!mkdir($local_data))
	{
		trigger_error("Could not create local data folder '{$local_data}'!", E_USER_ERROR);
	}
	if (!chmod($local_data, 0777))
	{
		trigger_error("Could not change privileges of local data folder '{$local_data}'!", E_USER_ERROR);
	}
}

//remove outdated genome FASTA files
if (file_exists("{$local_data}/{$build}.fa.md5"))
{
	exec("diff {$local_data}/{$build}.fa.md5 {$genome_folder}/{$build}.fa.md5", $output, $code);
	if ($code==0)
	{
		print "MD5 checksums of genome '{$build}.fa.md5' match.\n";
	}
	else
	{
		print "MD5 checksums of genome '{$build}.fa.md5' differ. Deleting old data!\n";
		exec2("rm -rf {$local_data}/{$build}.*");
	}
}

//wait sync is done by another processes
$pid_file = "{$local_data}/megSAP_data_setup_{$build}.txt";
if (file_exists($pid_file))
{
	$pid_old = trim(file_get_contents($pid_file));
	$iter = 0;
	while (posix_getpgid($pid_old)!==FALSE)
	{
		++$iter;
		if ($iter>180) break; //wait for 3h max
		print "Process with PID {$pid_old} is already syncing. Waiting one minute...\n";
		sleep(60);
		
		//update PID in case another process took over
		$pid_old = trim(file_get_contents($pid_file));
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
	$existed_before = file_exists($local_data.$base);
	list($stdout, $stderr) = exec2("{$rsync} {$genome_folder}{$base} {$local_data}{$base}");
	foreach(array_merge($stdout, $stderr) as $line)
	{
		$line = trim($line);
		if ($line=="") continue;
		print "    $line\n";
	}
	if (!$existed_before)
	{
		if (!chmod($local_data.$base, 0777))
		{
			trigger_error("Could not change privileges of local data folder '{$local_data}{$base}'!", E_USER_ERROR);
		}
	}
}

######################### VEP cache #########################
if ($build=="GRCh37")
{
	$annotation_folder = get_path("vep_data")."/cache/";
	$local_annotation_folder = "{$local_data}/".basename(get_path("vep_data"))."/";

	print "\n";
	print "### VEP cache ###\n";
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
		
		list($stdout, $stderr) = exec2("{$rsync} {$annotation_folder} {$local_annotation_folder}");
		foreach(array_merge($stdout, $stderr) as $line)
		{
			$line = trim($line);
			if ($line=="") continue;
			
			print "  {$line}\n";
		}
		
		touch("{$local_annotation_folder}/{$finished}");
	}
	else
	{
		print "No rsync necessary!\n";
	}
}

######################### VEP GeneSplicer data #########################
if ($build=="GRCh37")
{
	$annotation_folder = dirname(get_path("vep"))."/GeneSplicer/human/";
	$local_annotation_folder = "{$local_data}/GeneSplicer/";
	
	print "\n";
	print "### VEP GeneSplicer data ###\n";
	print "from: {$annotation_folder}\n";
	print "to  : {$local_annotation_folder}\n";
	print "\n";
		
	print "rsync-ing GeneSplicer data...\n";
	list($stdout, $stderr) = exec2("{$rsync} {$annotation_folder} {$local_annotation_folder}");
	foreach(array_merge($stdout, $stderr) as $line)
	{
			$line = trim($line);
			if ($line=="") continue;
			
			print "  {$line}\n";
	}
}


######################### VEP annotation databases #########################
if ($build=="GRCh37")
{
	$local_annotation_folder = "{$local_data}/ensembl-vep-dbs/";

	//create local folder if missing
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
	
	print "\n";
	print "### VEP annotation databases ###\n";
	print "to  : {$local_annotation_folder}\n";
	print "\n";
		
	print "rsync-ing annotation databases...\n";
	
	$db_files = array("/dbs/CADD/whole_genome_SNVs.tsv.gz", "/dbs/CADD/InDels.tsv.gz", "/dbs/REVEL/revel_all_chromosomes.tsv.gz", "/dbs/fathmm-MKL/fathmm-MKL_Current.tab.gz", "/dbs/dbscSNV/dbscSNV1.1_GRCh37.txt.gz", "/dbs/gnomAD/gnomAD_genome_r2.1.vcf.gz", "/dbs/RepeatMasker/RepeatMasker.bed.gz", "/dbs/ClinVar/clinvar_20190123_converted.vcf.gz", "/dbs/phyloP/hg19.100way.phyloP100way.bw", "/dbs/OMIM/omim.bed.gz", "/dbs/HGMD/HGMD_PRO_2018_4_fixed.vcf.gz");
	foreach($db_files as $db_file)
	{
		$source = $data_folder.$db_file;
		$target = $local_annotation_folder.basename($db_file);
		print "  rsync-ing database {$source}\n";
		list($stdout, $stderr) = exec2("{$rsync} {$source} {$target}");
		foreach(array_merge($stdout, $stderr) as $line)
		{
			$line = trim($line);
			if ($line=="") continue;
			
			print "    {$line}\n";
		}
		
		$source_index = $source.".tbi";
		if (file_exists($source_index))
		{
			print "  rsync-ing database index {$source_index}\n";
			
			$target_index = $target.".tbi";
			list($stdout, $stderr) = exec2("{$rsync} {$source_index} {$target_index}");
			foreach(array_merge($stdout, $stderr) as $line)
			{
				$line = trim($line);
				if ($line=="") continue;
				
				print "    {$line}\n";
			}
		}
	}
}

?>
