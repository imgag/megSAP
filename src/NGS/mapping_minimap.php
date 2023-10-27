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
$parser->addString("sample", "Sample name to use in BAM header. If unset the basename of the 'out' file is used.", true, "");
$parser->addInfileArray("in_fastq",  "Input file(s) in FASTQ format.", true);
$parser->addInfileArray("in_bam", "Input BAM file(s), with modified bases information (MM ML tags).", true);
$parser->addInfile("system",  "Processing system INI file (automatically determined from NGSD if 'sample' or 'out' is a valid processed sample name).", true);
$parser->addInt("threads", "Maximum number of threads used.", true, 2);
extract($parser->parse($argv));

if ((is_null($in_fastq) && is_null($in_bam)) || (count($in_fastq) > 0 && count($in_bam) > 0))
{
	trigger_error("Please specify either 'in_fastq' or 'in_bam'!", E_USER_ERROR);
}

//init vars
if($sample == "") $sample = basename2($out);
$basename = dirname($out)."/".$sample;
$bam_current = $parser->tempFile(".bam", $sample);

//extract processing system information from DB
$sys = load_system($system, $sample);

//set read group information
$group_props = [];
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

// alignment pipeline:
// BAM input available:
// samtools cat <input bams> | samtools fastq | minimap | zipperbams | samtools sort
// FASTQ input:
// zcat <input fastqs>                        | minimap              | samtools sort

$bam_input = count($in_bam) > 0;

$minimap_options = [
	"-a",
	"--MD",
	"-x map-ont",
	"--eqx",
	"-t {$threads}",
	"-R '@RG\t" . implode("\\t", $group_props) . "'",
	genome_fasta($sys['build'])
];

$pipeline = [];

//start from BAM input if available
if ($bam_input)
{
	$pipeline[] = [get_path("samtools"), "cat " . implode(" ", $in_bam)];
	$pipeline[] = [get_path("samtools"), "fastq -0 - -"];
}
else
{
	$pipeline[] = ["zcat", implode(" ", $in_fastq)];
}

//mapping with minimap2
$pipeline[] = [get_path("minimap2"), implode(" ", $minimap_options)];

//add tags from unmapped modified bases BAM
if ($bam_input)
{
	//TODO replace fgbio ZipperBams with a different tool
	//TODO fgbio ZipperBams does not support multiple BAM files, use concatenated version
	$pipeline[] = ["/mnt/storage2/users/ahadmaj1/.snakemake-conda/56f717cad63f6f01dc29688088922398_/bin/fgbio", "--compression 0 ZipperBams --unmapped {$in_bam[0]} --ref ".genome_fasta($sys['build'])];
}

//sort BAM by coordinates
$tmp_for_sorting = $parser->tempFile();
$pipeline[] = [get_path("samtools"), "sort -T {$tmp_for_sorting} -m 1G -@ ".min($threads, 4)." -o {$bam_current} -", true];

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
$qcml = $basename."_stats_map.qcML";
$params = ["-in $bam_current", "-out $qcml", "-ref ".genome_fasta($sys['build']), "-build ".ngsbits_build($sys['build'])];
if ($sys["target_file"]=="" || $sys["type"]=="lrGS")
{
	$params[] = "-wgs";
}
else
{
	$params[] = "-roi ".$sys["target_file"];
}
if ($sys['build']!="GRCh38")
{
	$params[] = "-no_cont";
}

$parser->exec(get_path("ngs-bits")."MappingQC", implode(" ", $params), true);

?>