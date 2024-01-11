<?php

/** 
	@page mapping_minimap_test
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("mapping_minimap_test", "Maps nanopore reads to a reference genome using minimap2.");
$parser->addInfileArray("in",  "Input file(s) in FASTQ/BAM format.", false);
$parser->addOutfile("out",  "Output file in BAM format (sorted).", false);
//optional
$parser->addString("sample", "Sample name to use in BAM header. If unset the basename of the 'out' file is used.", true, "");
$parser->addInfile("system",  "Processing system INI file (automatically determined from NGSD if 'sample' or 'out' is a valid processed sample name).", true);
$parser->addInt("threads", "Maximum number of threads used.", true, 2);
extract($parser->parse($argv));

//check input
if (count($in) < 1) trigger_error("No input files provided!", E_USER_ERROR);

$file_types = array();
foreach ($in as $file) 
{
	if (ends_with($file, ".bam") || ends_with($file, ".cram"))
	{
		$file_types[] = "bam";
	}
	elseif (ends_with($file, ".fq.gz") || ends_with($file, ".fastq.gz"))
	{
		$file_types[] = "fastq";
	}
	else
	{
		trigger_error("Invalid file type '{$file}' provided! Valid file types are FastQ and BAM/CRAM.", E_USER_ERROR);
	}
}
//remove duplicates
$file_types = array_unique($file_types);
if (count($file_types) > 1)
{
	trigger_error("Multiple input file types ('bam', 'fastq') provided!", E_USER_ERROR);
}
//set mode
$mode = $file_types[0];

//init vars
if($sample == "") $sample = basename2($out);
$basename = dirname($out)."/".$sample;
$bam_current = $parser->tempFile(".bam", $sample);

//extract processing system information from DB
$sys = load_system($system, $sample);

//set read group information
$group_props = array();
$group_props[] = "ID:{$sample}";
$group_props[] = "SM:{$sample}";
$group_props[] = "LB:{$sample}";
$group_props[] = "CN:medical_genetics_tuebingen";
$group_props[] = "DT:".date("c");
$group_props[] = "PL:ONT";
if(db_is_enabled("NGSD"))
{
	$db_conn = DB::getInstance("NGSD");
	$psample_info = get_processed_sample_info($db_conn, $sample, false);
	if ($psample_info!=NULL)
	{
		$group_props[] = "PM:".$psample_info['device_type'];
		$group_props[] = "en:".$psample_info['sys_name'];
	}
}


$pipeline = array();

if ($mode == "fastq")
{
	//fastq mapping with minimap2
	$pipeline[] = array(get_path("minimap2"), "--MD -Yax map-ont --eqx -t {$threads} -R '@RG\\t".implode("\\t", $group_props)."' ".genome_fasta($sys['build'])." ".implode(" ", $in));
}
else
{
	//bam mode: convert bam with samtools fastq
	$contains_methylation = false;
	foreach ($in as $file) 
	{
		//add methylation tag to output
		if (contains_methylation($file)) $met_tag = " -TMM,ML ";
	} 
	
	if (count($in) > 1)
	{
		//concat multiple bam files
		$pipeline[] = array(get_path("samtools"), "cat ".implode(" ", $in));
		$pipeline[] = array(get_path("samtools"), "fastq ".$met_tag);
	}
	else
	{
		$pipeline[] = array(get_path("samtools"), "fastq ".$met_tag.implode(" ", $in));
	}
	//perform mapping from STDIN
	$pipeline[] = array(get_path("minimap2"), "--MD -Yyax map-ont --eqx -t {$threads} -R '@RG\\t".implode("\\t", $group_props)."' ".genome_fasta($sys['build'])." - ");
}

//convert sam to bam and sort with samtools
$tmp_for_sorting = $parser->tempFile();
$pipeline[] = array(get_path("samtools"), "sort -T $tmp_for_sorting -m 1G -@ ".min($threads, 4)." -o $bam_current -", true);

//execute 
$parser->execPipeline($pipeline, "mapping");

//create index
$parser->indexBam($bam_current, $threads);

//copy BAM to final output location
$parser->copyFile($bam_current, $out);
$parser->copyFile($bam_current.".bai", $out.".bai");

// check if copy was successful
if (!file_exists($out) || filesize($bam_current) != filesize($out))
{
	trigger_error("Error during coping BAM file! File sizes don't match!", E_USER_ERROR);
}

//run read/mapping QC
$read_qc_file = $basename."_stats_fastq.qcML";
$mapping_qc_file = $basename."_stats_map.qcML";
$params = array("-in $bam_current", "-long_read", "-out {$mapping_qc_file}", "-read_qc {$read_qc_file}", "-ref ".genome_fasta($sys['build']), "-build ".ngsbits_build($sys['build']));
if ($sys['target_file']=="" || $sys['type']=="lrGS")
{
	$params[] = "-wgs";
}
else
{
	$params[] = "-roi ".$sys['target_file'];
}
if ($sys['build']!="GRCh38")
{
	$params[] = "-no_cont";
}

$parser->exec(get_path("ngs-bits")."MappingQC", implode(" ", $params), true);


?>
