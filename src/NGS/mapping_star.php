<?php

/**
	@page mapping_star
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("mapping_star", "Alignment of RNA-seq FASTQ files to a reference genome.");
//mandatory arguments
$parser->addInfile("in1", "Input forward reads in FASTQ(.GZ) format.", false);
$parser->addOutfile("out", "Output BAM file name.", false);

//optional arguments
$parser->addInfile("in2",  "Input reverse reads in FASTQ(.GZ) format for paired-end alignment.", true);
$parser->addString("genome", "STAR reference genome index.", true, get_path("data_folder")."genomes/STAR/GRCh37");

$parser->addFlag("uncompressed", "FASTQ input files are uncompressed.");

$parser->addFlag("long_reads", "Use STAR version suitable for very long reads > 500nt.");
$parser->addInt("threads", "Number of parallel threads.", true, 4);

$parser->addFlag("unstranded_xs", "For unstranded data, add the XS strand attribute for spliced alignments for cufflinks compatibility. Note: Alignments with undefined strandedness will be removed!");
$parser->addFlag("skip_dedup", "Skip alignment duplication marking.");
$parser->addFlag("no_splicing", "Prevent reads from getting spliced");

$parser->addFlag("all_junctions", "Disable filtering of reported junctions (affects splicing output only, not alignment).");

$parser->addInt("sj_overhang", "Minimum overhang for non-annotated splice junctions.", true, 8);
$parser->addInt("sjdb_overhang", "Minimum overhang for annotated splice junctions.", true, 1);
$parser->addInt("Lmax", "Max length of read fraction for seed search", true, 50);

extract($parser->parse($argv));

$outdir = realpath(dirname($out))."/";
$prefix = $outdir.basename($out, ".bam");

// Temporary prefix where STAR stores all files.
$STAR_tmp_folder = $parser->tempFolder();

//build command
//add a read group line to the final BAM file
$group_props = array();
$basename = basename($out, ".bam");
list($gs, $ps) = explode("_", $basename."_99");
$group_props[] = "ID:{$basename}";
$group_props[] = "SM:{$gs}";
$group_props[] = "LB:{$gs}_{$ps}";
$group_props[] = "CN:medical_genetics_tuebingen";
$group_props[] = "DT:".date("c");
$group_props[] = "PL:ILLUMINA";
$psample_info = get_processed_sample_info($basename,false);
if (!is_null($psample_info))
{
	$group_props[] = "PM:".$psample_info['device_type'];
	$group_props[] = "en:".$psample_info['sys_name'];
}

$arguments = array(
	"--readFilesIn", isset($in2) ? "{$in1} {$in2}" : $in1,
	"--genomeDir {$genome}",
	"--outFileNamePrefix {$STAR_tmp_folder}/",
	"--outStd BAM_Unsorted",
	"--outSAMtype BAM Unsorted",
	"--outSAMunmapped Within",
	"--runThreadN {$threads}",
	"--outSAMattributes All",
	"--chimOutType WithinBAM SoftClip",
	"--chimSegmentMin 12",
	"--chimJunctionOverhangMin 12",
	"--chimSegmentReadGapMax 3",
	"--seedSearchStartLmax $Lmax",
	"--alignMatesGapMax 1000000",
	"--alignSJoverhangMin $sj_overhang",
	"--alignSJDBoverhangMin $sjdb_overhang",
	"--alignSJstitchMismatchNmax 5 -1 5 5",
	"--outSAMattrRGline", "'".implode(" ", $group_props)."'"
);
	
if ($unstranded_xs) $arguments[] = "--outSAMstrandField intronMotif";
if (!$uncompressed) $arguments[] = "--readFilesCommand zcat";
if ($all_junctions) $arguments[] = "--outSJfilterDistToOtherSJmin 0 0 0 0 --outSJfilterOverhangMin 1 1 1 1 --outSJfilterCountUniqueMin 1 1 1 1 --outSJfilterCountTotalMin 1 1 1 1";

if ($no_splicing)
{
	//min has to be larger than max to prevent spliced alignments
	$arguments[] = "--alignIntronMin 2";
	$arguments[] = "--alignIntronMax 1";
	$arguments[] = "--twopassMode None";
}
else
{
	$arguments[] = "--alignIntronMax 1000000 --alignIntronMin 20";
	$arguments[] = "--twopassMode Basic";
}

//STAR or STARlong program
$star = $long_reads ? get_path("STAR")."long" : get_path("STAR");

//mapping with STAR
$pipeline = array();
$pipeline[] = array($star, implode(" ", $arguments));

$pipeline[] = array(get_path("samtools"), "view -h");

//duplicate flagging with samblaster
if (!$skip_dedup) $pipeline[] = array(get_path("samblaster"), "");

//sort BAM by coordinates
$tmp_for_sorting = $parser->tempFile();
$pipeline[] = array(get_path("samtools"), "sort -T $tmp_for_sorting -m 1G -@ ".min($threads, 4)." -o $out -", true);

//execute (STAR -> samblaster -> samtools SAM to BAM -> samtools sort)
$parser->execPipeline($pipeline, "mapping");

//create BAM index file
$parser->indexBam($out, $threads);

//downstream analysis files
$outfile_splicing= "{$prefix}_splicing.tsv";
$splicing_header = array(
	"#chr",
	"start",
	"end",
	"strand",
	"intron_motif",
	"annotated",
	"n_unique_reads",
	"n_multi_reads",
	"max_overhang"
);
file_put_contents($outfile_splicing, implode("\t", $splicing_header)."\n");
file_put_contents($outfile_splicing, file_get_contents("{$STAR_tmp_folder}/SJ.out.tab"), FILE_APPEND);

$outfile_chimeric = "{$prefix}_chimeric.tsv";
$chimeric_header = array(
	"#donor_chr",
	"donor_start",
	"donor_strand",
	"acceptor_chr",
	"acceptor_start",
	"acceptor_strand",
	"junction_type",
	"left_repeat_len",
	"right_repeat_len",
	"read_name",
	"seg1_first_base",
	"seg1_cigar",
	"seg2_first_base",
	"seg2_cigar"
);
file_put_contents($outfile_chimeric, implode("\t", $chimeric_header)."\n");
file_put_contents($outfile_chimeric, file_get_contents("{$STAR_tmp_folder}/Chimeric.out.junction"), FILE_APPEND);

//write the final log file into the tool log
$final_log = "{$STAR_tmp_folder}/Log.final.out";
$parser->log("STAR Log.final.out", file($final_log));

?>