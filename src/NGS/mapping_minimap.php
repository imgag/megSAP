<?php

/** 
	@page mapping_minimap
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("mapping_minimap", "Maps nanopore reads to a reference genome using minimap2.");

$parser->addOutfile("out",  "Output file in BAM format (sorted).", false);
//optional
$parser->addInfile("in_bam",  "Input BAM file.", true);
$parser->addInfileArray("in_fastq",  "Input file(s) in FASTQ format.", false);
$parser->addString("sample", "Sample name to use in BAM header. If unset the basename of the 'out' file is used.", true, "");
$parser->addInfile("system",  "Processing system INI file (automatically determined from NGSD if 'sample' or 'out' is a valid processed sample name).", true);
$parser->addInt("threads", "Maximum number of threads used.", true, 2);
extract($parser->parse($argv));


//init vars
if($sample == "") $sample = basename2($out);
$basename = dirname($out)."/".$sample;
print $basename;
$bam_current = $parser->tempFile(".bam", $sample);

//extract processing system information from DB
$sys = load_system($system, $sample);

//check input
if(is_null($in_fastq) && is_null($in_bam)) trigger_error("Either BAM or FastQ files have to be provided!", E_USER_ERROR);
if(!is_null($in_fastq) && !is_null($in_bam)) trigger_error("You have to provide either BAM or FastQ files, not both!", E_USER_ERROR);

if(!is_null($in_fastq))
{
	$in = $in_fastq;
	$bam_modus = false;
}
else //BAM provided
{
	$tmp1 = $parser->tempFile(".fastq.gz");
	$parser->exec("{$ngsbits}BamToFastq", "-mode single-end -in $bam_file -out1 $tmp1", true);
	$in = array($tmp1);
	$bam_modus = true;
}


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
print get_path("minimap2")." --MD -Yax map-ont --eqx -t {$threads} -R '@RG\\t".implode("\\t", $group_props)."' ".genome_fasta($sys['build'])." ".implode(" ", $in);

//mapping with minimap2
$pipeline[] = array(get_path("minimap2"), " --MD -Yax map-ont --eqx -t {$threads} -R '@RG\\t".implode("\\t", $group_props)."' ".genome_fasta($sys['build'])." ".implode(" ", $in));

//annotate methylation
if($bam_modus)
{
	$pipeline[] = array(get_path("fgbio"), "--compression 0 ZipperBams -u {$in_bam} -r ".genome_fasta($sys['build']));
}

//sort BAM by coordinates
$tmp_for_sorting = $parser->tempFile();
$pipeline[] = array(get_path("samtools"), "sort -T $tmp_for_sorting -m 1G -@ ".min($threads, 4)." -o $bam_current -", true);

//execute 
$parser->execPipeline($pipeline, "mapping");

//create index
$parser->indexBam($bam_current, $threads);

//TODO: test
//check readcounts
// if($bam_modus)
// {
// 	list($input_readcount) = $parser->exec(get_path("samtools"), "view -F 256 -c -@ {$threads} {$in_bam}");
// 	$input_readcount = (int) $input_readcount;
// }
// else
// {
// 	$pipeline = array();
// 	$pipeline[] = array("zcat", implode(" ", $in_fastq));
// 	$pipeline[] = array("wc", "-l");
// 	list($stdout, $stderr) = $parser->execPipeline($pipeline, "FastQ read counts", true);
// 	$input_readcount = ((int) $stdout[0])/4;
// }
// list($stdout, $stderr) = $parser->exec(get_path("samtools"), "view -F 256 -c -@ {$threads} {$bam_current}");
// $output_readcount = (int) $stdout[0];
// if($input_readcount != $output_readcount) trigger_error("Read count of input and output differs! \n Input:\t{$input_readcount}\nOutput:\t{$output_readcount}", E_USER_ERROR);


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
