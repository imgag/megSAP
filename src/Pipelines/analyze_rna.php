<?php

/**
	@page analyze_rna
	
	@todo test if STAR-Fusion and STAR can use the genome reference, perhaps even the same genome as for DNA
	@todo check when duplicate removal is really necessary for RNA => try samblaster instead of picard (we use it for DNA)
	@todo check if/when indel realignment is really necessary for RNA => ABRA as indel realigner instead of GATK (license problems!)
	@todo if we continue using GATK for indel realignment: check if GATK/hg19.fa can be used as GATKReference instead of GATK/hg19/hg19_GATK.fa
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("analyze_rna", "RNA Mapping pipeline using STAR.");
$parser->addInfileArray("in_for", "Forward reads in fastq.gz file(s).", false);
$parser->addString("out_folder", "Output folder.", false);
$parser->addString("out_name", "Output base file name, typically the processed sample ID (e.g. 'GS120001_01').", false, NULL);
//optional
$parser->addInfileArray("in_rev", "Reverse reads fastq.gz file(s) for paired end alignment.", true);
$parser->addInfile("system", "Processing system INI file (Determined from 'out_name' by default).", true);
$steps_all = array("ma", "rc", "an", "fu", "db");
$parser->addString("steps", "Comma-separated list of processing steps to perform.", true, implode(",", $steps_all));
$parser->addString("reference", "FASTA File that contains the reference genome used for mapping.", true, get_path("data_folder")."genomes/STAR/hg19");
$parser->addString("fusionDetectionReference", "Reference genome used for STAR-Fusion", true, get_path("data_folder")."/genomes/STAR-Fusion/hg19");
$parser->addString("gtfFile", "GTF File containing feature annotations used for read counting.", true, get_path("data_folder")."genomes/gtf/ucsc_refseq_hg19.gtf");
$parser->addString("GATKReference", "GATK reference fasta file", true, get_path("data_folder")."genomes/GATK/hg19/hg19_GATK.fa");
$parser->addString("featureType", "Feature type used for mapping reads to features.", true, "exon");
$parser->addString("gtfAttribute", "GTF attribute used as feature ID.", true, "gene_id");
$parser->addInt("threads", "The maximum number of threads used.", true, 4);
$parser->addFlag("indelRealign", "Perform indel realignment. By default is is skipped.");
$parser->addFlag("stranded", "Specify whether a stranded protocol was used during library preparation. Default is non-stranded.");
$parser->addFlag("keepUnmapped", "Save unmapped reads as fasta files.");
$parser->addFlag("rmdup", "Remove PCR duplicates after mapping.");
$parser->addFlag("sharedMemory", "Use shared memory for running STAR alignment jobs.");

extract($parser->parse($argv));

//check steps
$steps = explode(",", $steps);
foreach($steps as $step)
{
	if (!in_array($step, $steps_all))
	{
		trigger_error("Unknown processing step '$step'!", E_USER_ERROR);
	}
}

//init
$prefix = $out_folder."/".$out_name;
$sys = load_system($system, $out_name);
$paired = isset($in_rev);

//mapping and QC
$final_bam = $prefix.".bam";
$qc_fastq = $prefix."_stats_fastq.qcML";
$qc_map = $prefix."_stats_map.qcML";
if(in_array("ma", $steps))
{	
	//check FASTQ quality encoding
	$files = $paired ? array_merge($in_for, $in_rev) : $in_for;
	foreach($files as $file)
	{
		list($stdout, $stderr) = $parser->exec(get_path("ngs-bits")."FastqFormat", "-in $file", true);
		if (!contains($stdout[2], "Sanger"))
		{
			trigger_error("Input file '$in_for' is not in Sanger/Illumina 1.8 format!", E_USER_ERROR);
		}
	}

	//check that adapters are specified
	if ($sys["adapter1_p5"]=="" || $sys["adapter2_p7"]=="")
	{
		trigger_error("No forward and/or reverse adapter sequence given!\nForward: ".$sys["adapter1_p5"]."\nReverse: ".$sys["adapter2_p7"], E_USER_ERROR);
	}
	
	//adapter trimming + QC (SeqPurge for paired-end, Skewer/ReadQC for single-end)
	if($paired)
	{
		$fastq_trimmed1 = $parser->tempFile("_trimmed.fastq.gz");
		$fastq_trimmed2 = $parser->tempFile("_trimmed.fastq.gz");
		$parser->exec(get_path("ngs-bits")."SeqPurge", "-in1 ".implode(" ", $in_for)." -in2 ".implode(" ", $in_rev)." -out1 $fastq_trimmed1 -out2 $fastq_trimmed2 -a1 ".$sys["adapter1_p5"]." -a2 ".$sys["adapter2_p7"]." -qc $qc_fastq -threads ".$threads, true);
	}
	else
	{
		$parser->exec(get_path("ngs-bits")."ReadQC", "-in1 ".implode(" ", $in_for)." -out $qc_fastq", true);
		
		$fastq_trimmed1 = $parser->tempFile("_trimmed.fastq.gz");
		$fastq_trimmed_basename = substr($fastq_trimmed1, 0, -9);
		$parser->exec("zcat", implode(" ", $in_for)." | ".get_path("skewer")." -x ".$sys["adapter1_p5"]." -y ".$sys["adapter2_p7"]." -m any -threads $threads -z -o $fastq_trimmed_basename --quiet -", true);
		rename($fastq_trimmed_basename."-trimmed.fastq.gz", $fastq_trimmed1);
		
		//remove the log file of skewer
		unlink($fastq_trimmed_basename."-trimmed.log");
	}

	//mapping
	$args = array("-out $final_bam", "-genome $reference", "-p $threads", "-in1 $fastq_trimmed1");
	if($paired) $args[] = "-in2 $fastq_trimmed2";
	if($keepUnmapped) $args[] = " -keepUnmapped";
	if($rmdup) $args[] = "-rmdup";
	if($sharedMemory)
	{
		if (in_array("fu", $steps)) //for STAR 2.5.2b this does not work, but it might change in newer STAR releases
		{
			$parser->log("Using shared memory and detecting fusion proteins is not possible at the same time. Disabling shared memory.");
		}
		else
		{
			$args[] = "-useSharedMemory";
		}
	}
	$parser->execTool("NGS/mapping_star_htseq.php", implode(" ", $args));
	
	//indel realignment
	if($indelRealign)
	{
		//index
		$parser->exec(get_path("samtools"), "index $final_bam", true);
		
		//split and trim reads with n in the cigar string -> spliced reads
		$bam_split = $parser->tempFile(".sorted.split.bam");
		$parser->exec(get_path("GATK"), "-T SplitNCigarReads -R $GATKReference -I $bam_mapped -o $bam_split -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS", true);
		
		//perform indel realignment on split reads
		$parser->execTool("NGS/indel_realign.php", "-in $bam_split -out $final_bam");
	}
}

//read counting
$counts_raw = $prefix."_counts_raw.tsv";
$counts_fpkm = $prefix."_counts_fpkm.tsv";
if(in_array("rc", $steps))
{
	$args = array();
	if($paired) $args[] = "-paired";
	if($stranded) $args[] = "-stranded";
	$parser->execTool("NGS/read_counting_featureCounts.php", "-in $final_bam -out $counts_raw -threads $threads -gtfFile $gtfFile -featureType $featureType -gtfAttribute $gtfAttribute ".implode(" ", $args));

	//normalize expression values
	$counts_tmp = $parser->tempFile("_counts_tmp.tsv"); //remove the header line and just take the gene ID and counts column
	$parser->exec("tail", "-n +2 $counts_raw | cut -f 1,7 > $counts_tmp", false);
	//normalize read counts
	$parser->execTool("NGS/normalize_read_counts.php", "-in $counts_tmp -out $counts_fpkm -gtf $gtfFile -method rpkm -feature $featureType -idattr $gtfAttribute -header");
}

//annotate
if(in_array("an", $steps))
{
	$parser->execTool("NGS/annotate_count_file.php", "-in $counts_fpkm -out $counts_fpkm");
}

//detect fusions
if(in_array("fu",$steps))
{
	$parser->exec(get_path("STAR-Fusion"), "--genome_lib_dir $fusionDetectionReference -J {$prefix}Chimeric.out.junction --output_dir $out_folder", true);
	
	//cleanup
	$parser->exec("rm -r", "$out_folder/star-fusion.filter.intermediates_dir", false);
	$parser->exec("rm -r", "$out_folder/star-fusion.predict.intermediates_dir", false);
	unlink("$out_folder/star-fusion.fusion_candidates.final");
	unlink("$out_folder/star-fusion.fusion_candidates.preliminary");
	unlink("$out_folder/star-fusion.fusion_candidates.preliminary.wSpliceInfo");
	unlink("$out_folder/star-fusion.fusion_candidates.preliminary.wSpliceInfo.ok");
	unlink("$out_folder/star-fusion.STAR-Fusion.filter.ok");
	unlink("$out_folder/star-fusion.STAR-Fusion.predict.ok");
	rename("$out_folder/star-fusion.fusion_candidates.final.abridged", $prefix."_var_fusions.tsv");
}

//import to database
$log_db  = $prefix."_log_db.log";
if (in_array("db", $steps))
{
	if(file_exists($log_db)) unlink($log_db);
	$parser->execTool("NGS/db_check_gender.php", "-in $final_bam -pid $out_name");
		
	//update last_analysis column of processed sample in NGSD
	updateLastAnalysisDate($out_name, $final_bam);

	//import QC data
	$parser->execTool("NGS/db_import_qc.php", "-id $out_name -files $qc_fastq $qc_map -force");
}
?>
