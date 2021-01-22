<?php
/** 
	@page export_raw_data 
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("export_raw_data", "Exports RAW data.");
$parser->addStringArray("samples", "Processed sample names (or file with one sample per line).", false);
$parser->addString("out", "Output folder and ZIP file name.", false);
$parser->addEnum("mode", "Export mode (FASTQ only, BAM only, whole analysis folder).", true, ["fastq", "bam", "folder"], "fastq");
$parser->addFlag("internal", "Use internal webserver and do not use password.");
extract($parser->parse($argv));

//init
$db = DB::getInstance("NGSD");
$ngsbits = get_path("ngs-bits");

//check folder
if (!file_exists($out) || !is_dir($out))
{
	print "Output folder '$out' is missing - creating it ...\n";
	mkdir($out, 0777, true);
}

//load samples from file
if (count($samples)==1 && file_exists($samples[0]))
{
	$file = file($samples[0]);
	$samples = array();
	foreach($file as $line)
	{
		$line = trim($line);
		if ($line=="" || $line[0]=='#') continue;
		$samples[] = $line;
	}
}

//check and create meta data file
print "Checking ".count($samples)." samples ...\n";
$meta = [];
$meta[] = "#sample\texternal_name\tgender\tdisease_group_and_status\tkit\tsequencer\n";
foreach($samples as $ps)
{
	$info = get_processed_sample_info($db, $ps);
	$quality = $info['ps_quality'];
	if ($quality=="bad") trigger_error("Sample '$ps' has 'bad' quality!", E_USER_ERROR);
	
	if ($mode=="folder")
	{
		$folder = $info['ps_folder'];
		if (!file_exists($folder))
		{
			trigger_error("Sample folder '$folder' is missing!", E_USER_ERROR);
		}
	}
	else if ($mode=="bam")
	{
		$bam = $info['ps_bam'];
		if (!file_exists($bam))
		{
			trigger_error("Sample '$ps': BAM file is missing!", E_USER_ERROR);
		}
	}
	else
	{
		$bam = $info['ps_bam'];
		$fastqs = glob($info['ps_folder']."/*.fastq.gz");
		if (count($fastqs)==0 && !file_exists($bam))
		{
			trigger_error("Sample '$ps': BAM file is missing and no FASTQs found!", E_USER_ERROR);
		}
	}
	
	$gender = $info['gender'];
	$external = $info['name_external'];
	$device = $info['device_type'];
	$system = $info['sys_name'];
	$d_group = $info['disease_group'];
	$d_status = $info['disease_status'];
	$meta[] = "{$ps}\t{$external}\t{$gender}\t{$d_group} ({$d_status})\t{$system}\t{$device}\n";
}
file_put_contents("$out/meta_data.tsv", $meta);

//export data
foreach($samples as $ps)
{
	print "$ps\n";
	$info = get_processed_sample_info($db, $ps);
	
	if ($mode=="folder")
	{
		$folder = $info['ps_folder'];
		exec2("ln -s {$folder} {$out}/Sample_{$ps}");
	}
	else if ($mode=="bam")
	{
		print "  Copying BAM file ...\n";
		$bam = $info['ps_bam'];
		exec2("ln -s {$bam} {$out}/".basename($bam));
	}
	else
	{		
		$bam = $info['ps_bam'];
		$fastqs = glob($info['ps_folder']."/*.fastq.gz");
		if (count($fastqs)>0)
		{
			print "  Copying FASTQ files ...\n";
			foreach($fastqs as $fastq)
			{
				exec2("ln -s {$fastq} {$out}/".basename($fastq));
			}
		}
		else
		{
			print "  Generating FASTQ files from BAM ...\n";
			exec2("{$ngsbits}/BamToFastq -in {$bam} -out1 {$out}/{$ps}_R1_001.fastq.gz -out2 {$out}/{$ps}_R2_001.fastq.gz");
		}
	}
}

//determine password
$password = "";
if(!$internal)
{
	$password = random_string(20);
}
 
//zip
print "Zipping output folder".($password=="" ? "" : " using the password '{$password}'")."...\n";
exec2("zip -1 ".($password=="" ? "" : "--password {$password} ")."-r {$out}.zip $out");

//move
print "You can move the file to the webserver using:\n";
$user = exec('whoami');
if($internal)
{
	print "  > mv {$out}.zip /mnt/users/bioinf/http_shareukt/{$user}/\n";
}
else
{
	print "  > scp {$out}.zip {$user}@imgag.de:/var/www/html/download/{$user}/\n";
}

print "The URL of the file is:\n";
if($internal)
{
	print "  http://srv018.img.med.uni-tuebingen.de:8080/{$user}/{$out}.zip\n";
}
else
{
	print "  https://download.imgag.de/{$user}/{$out}.zip\n";
}


?>