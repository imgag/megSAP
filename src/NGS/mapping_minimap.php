<?php

/** 
	@page mapping_minimap
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("mapping_minimap", "Maps nanopore reads to a reference genome using minimap2.");
$parser->addInfileArray("in",  "Input file(s) in FASTQ format.", false);
$parser->addOutfile("out",  "Output file in BAM format (sorted).", false);
//optional
$parser->addString("sample", "Sample name to use in BAM header. If unset the basename of the 'out' file is used.", true, "");
$parser->addInfile("system",  "Processing system INI file (automatically determined from NGSD if 'sample' or 'out' is a valid processed sample name).", true);
$parser->addInt("threads", "Maximum number of threads used.", true, 2);
extract($parser->parse($argv));


//init vars
if($sample == "") $sample = basename($out, ".bam");
$basename = dirname($out)."/".$sample;
print $basename;
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

//debug
print get_path("minimap2")." --MD -ax map-ont --eqx -t {$threads} -R '@RG\\t".implode("\\t", $group_props)."' ".genome_fasta($sys['build'])." ".implode(" ", $in);

//mapping with minimap2
$pipeline[] = array(get_path("minimap2"), " --MD -ax map-ont --eqx -t {$threads} -R '@RG\\t".implode("\\t", $group_props)."' ".genome_fasta($sys['build'])." ".implode(" ", $in));

//convert sam to bam with samtools
$tmp_unsorted = $parser->tempFile("_unsorted.bam");

//sort BAM by coordinates
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

//run mapping QC
$stafile2 = $basename."_stats_map.qcML";
$params = array("-in $bam_current", "-out $stafile2", "-ref ".genome_fasta($sys['build']), "-build ".ngsbits_build($sys['build']));
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
