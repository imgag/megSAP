<?php

/** 
	@page mapping_minimap
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("mapping_minimap", "Maps reads to a reference genome using minimap2.");
$parser->addOutfile("out",  "Output file in BAM format (sorted).", false);
//optional
$parser->addString("sample", "Sample name to use in BAM header. If unset the basename of the 'out' file is used.", true, "");
$parser->addInfileArray("in_fastq",  "Input file(s) in FASTQ format.", true);
$parser->addInfileArray("in_bam", "Input BAM file(s), with modified bases information (MM ML tags).", true);
$parser->addInfile("system",  "Processing system INI file (automatically determined from NGSD if 'sample' or 'out' is a valid processed sample name).", true);
$parser->addString("qc_fastq", "Output qcML file with read statistics.", true, "");
$parser->addString("qc_map", "Output qcML file with mapping statistics.", true, "");
$parser->addOutfile("local_bam", "Filename the local BAM file is written to. It can be used to speed up variant calling etc. Not created if unset.", true);
$parser->addFlag("bam_output", "Output is BAM instead of CRAM.");
$parser->addInt("threads", "Maximum number of threads used.", true, 2);
$parser->addFlag("softclip_supplements", "Use soft clipping for supplementary alignments.");
extract($parser->parse($argv));

if ((is_null($in_fastq) && is_null($in_bam)) || (!is_null($in_fastq) && (count($in_fastq) > 0) && !is_null($in_bam) && (count($in_bam) > 0)))
{
	trigger_error("Please specify either 'in_fastq' or 'in_bam'!", E_USER_ERROR);
}

//init vars
if($sample == "") $sample = basename2($out);
$basename = dirname($out)."/".$sample;
$bam_current = $parser->tempFile(".bam", $sample);
$qcml_reads = $basename."_stats_fastq.qcML";
$qcml_map = $basename."_stats_map.qcML";


//extract processing system information from DB
$sys = load_system($system, $sample);
$genome_fasta = genome_fasta($sys['build']);

//set read group information
$group_props = [];
$group_props[] = "ID:{$sample}";
$group_props[] = "SM:{$sample}";
$group_props[] = "LB:{$sample}";
$group_props[] = "DT:".date("c");
//default PL to ONT, change to PacBio if necessary
if ($sys['name_short'] == "LR-PB-SPRQ")
{
	$group_props[] = "PL:PACBIO";
}
else
{
	$group_props[] = "PL:ONT";
}
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

//set preset, default to ONT
if ($sys['name_short'] == "LR-PB-SPRQ")
{
	$preset = "map-hifi";
}
else
{
	$preset = "map-ont";
}


// alignment pipeline:
// BAM input available:
// samtools cat <input bams> | samtools fastq | minimap | samtools sort
// FASTQ input:
// zcat <input fastqs>                        | minimap | samtools sort

$bam_input = !is_null($in_bam) && (count($in_bam) > 0);

if ($bam_input)
{
	//get read group description (base calling) from input bams
	$rg_description = array();
	foreach ($in_bam as $bam_file) 
	{
		$tmp = get_read_group_description($bam_file);
		$rg_description = array_merge($rg_description, $tmp);	
	}
	//remove duplicates
	$rg_description = array_unique($rg_description);

	//add to '@RG' entry
	foreach ($rg_description as $description) 
	{
		$group_props[] = $description;
	}
}

$minimap_options = [
	"-a",
	"--MD",
	"-x {$preset}",
	"--eqx",
	"-t {$threads}",
	"-R '@RG\\t" . implode("\\t", $group_props) . "'",
	$genome_fasta
];
//include methylation
if ($bam_input) $minimap_options[] = "-y";

//soft-clip supplement alignments
if ($softclip_supplements) $minimap_options[] = "-Y";

$pipeline = [];

//start from BAM input if available
$tmp_fastq = $parser->tempFile(".fastq.gz");
if ($bam_input)
{
	//bam mode: convert bam with samtools fastq
	$met_tag = "";
	foreach ($in_bam as $file) 
	{
		//add methylation tag to output
		if (contains_methylation($file, 100, $sys['build']))
		{
			$met_tag = " -TMM,ML ";
			break;
		} 
	} 
	// make separate calls of samtools fastq in a subshell
	// samtools cat <aligned.bam> <unaligned.bam> creates invalid BAMs due to different headers
	$fastq_cmds = [];
	foreach ($in_bam as $file)
	{
		$command = $parser->execApptainer("samtools", "samtools", "fastq --reference {$genome_fasta} -o /dev/null {$met_tag} {$file}", [$file, $genome_fasta], [], true);
		$fastq_cmds[] = $command;
	}
	$fastq_cmds_str = implode("; ", $fastq_cmds);
	$pipeline[] =  ["", "({$fastq_cmds_str})"];	//perform mapping from STDIN
	$command = $parser->execApptainer("minimap2", "minimap2", implode(" ", $minimap_options)." - ", [genome_fasta($sys['build'])], [], true);
	$pipeline[] = ["", $command];
}
else //fastq_mode
{	
	$in_files = array();
	$in_files = array_merge($in_files, $in_fastq);
	$in_files[] = genome_fasta($sys['build']);

	//FastQ mapping
	$command = $parser->execApptainer("minimap2", "minimap2", implode(" ", $minimap_options)." ".implode(" ", $in_fastq), $in_files, [], true);
	$pipeline[] = ["", $command];
}

//sort BAM by coordinates
$tmp_for_sorting = $parser->tempFile();
$command = $parser->execApptainer("samtools", "samtools", "sort --reference {$genome_fasta} -T {$tmp_for_sorting} -m 1G -@ ".min($threads, 4)." -o {$bam_current} -", [$genome_fasta], [], true);
$pipeline[] = ["", $command];//execute 
$parser->execPipeline($pipeline, "mapping");

//create index
$parser->indexBam($bam_current, $threads);

//check for missing reads
if ($bam_input)
{
	//get read count of input BAMs
	$rc_input_bam = 0;
	foreach ($in_bam as $bam_file) 
	{
		$rc_input_bam += get_read_count($bam_file, max(8, $threads), array("-F", "2304"), $sys['build']);
	}
	//get read count of output BAM
	$rc_output_bam = get_read_count($bam_current, max(8, $threads), array("-F", "2304"), $sys['build']);

	// get rel. difference
	$diff = floatval(abs($rc_input_bam - $rc_output_bam));
	$rel_diff = floatval($diff /(($rc_input_bam + $rc_output_bam)/2.0));

	//allow 0.01% difference between input and output 
	if ($rel_diff > 0.0001) 
	{
		trigger_error("Mapped read counts doesn't match input read counts! (Input: {$rc_input_bam} reads vs. Output: {$rc_output_bam} reads (rel. difference: ".($rel_diff*100.0)."%))", E_USER_ERROR);
	}
	else
	{
		trigger_error("Input BAMs: {$rc_input_bam} reads vs. Output BAMs: {$rc_output_bam} reads (rel. difference: ".($rel_diff*100.0)."%)", E_USER_NOTICE);
	}
}

//run mapping QC
if ($qc_map !== "")
{
	$in_files = array();
	$in_files[] = genome_fasta($sys["build"]);

	$params = [
		"-in $bam_current",
		"-out $qcml_map",
		"-read_qc $qcml_reads",
		"-ref ".genome_fasta($sys["build"]),
		"-build ".ngsbits_build($sys["build"]),
		"-long_read"
	];

	if ($sys['target_file']=="" || $sys['type']=="lrGS")
	{
		$params[] = "-wgs";
	}
	else
	{
		$params[] = "-roi ".realpath($sys['target_file']);
		$in_files[] = $sys['target_file'];
	}
	if ($sys['build']!="GRCh38")
	{
		$params[] = "-no_cont";
	}

	$parser->execApptainer("ngs-bits", "MappingQC", implode(" ", $params), $in_files, [dirname($out)]);
}

//create CRAM/BAM in output folder
if ($bam_output)
{
	//copy BAM to final output location
	$parser->copyFile($bam_current, $out);
	$parser->copyFile($bam_current.".bai", $out.".bai");

	// check if copy was successful
	if (!file_exists($out) || filesize($bam_current) != filesize($out))
	{
		trigger_error("Error during coping BAM file! File sizes don't match!", E_USER_ERROR);
	}
}
else
{
	$out_cram = $basename.".cram";
	$parser->execTool("Tools/bam_to_cram.php", "-bam {$bam_current} -cram {$out} -build ".$sys['build']." -threads {$threads}");

	//check if conversion was successful
	compare_bam_read_count($bam_current, $out_cram, max(8, $threads), true, true, 0.0, array(), $sys['build']);
	
	//in case we re-map an old analysis with BAM output, we need to delete the BAM file
	$bam = basename2($out).".bam";
	if (file_exists($bam)) unlink($bam);
	$bai = basename2($out).".bam.bai";
	if (file_exists($bai)) unlink($bai);
}

// rename tmp BAM to allow using it for variant calling etc
if (!empty($local_bam))
{
	$parser->moveFile($bam_current, $local_bam);
	$parser->moveFile($bam_current.".bai", $local_bam.".bai");
}


?>