<?php

/**
 * @page analyze_rna
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("analyze_rna", "RNA mapping pipeline using STAR.");

//mandatory
$parser->addString("folder", "Analysis data folder.", false);
$parser->addString("name", "Base file name, typically the processed sample ID (e.g. 'GS120001_01').", false);

//optional
$parser->addInfile("system", "Processing system INI file (determined from NGSD via the 'name' by default).", true);
$steps_all = array("ma", "rc", "an", "fu", "db");
$parser->addString("steps", "Comma-separated list of steps to perform:\nma=mapping, rc=read counting, an=annotation, fu=fusion detection, db=import into NGSD", true, implode(",", $steps_all));

$parser->addEnum("library_type", "Specify the library type, i.e. the strand R1 originates from (dUTP libraries correspond to reverse).", true, array("unstranded", "reverse", "forward"), "reverse");
$parser->addFlag("no_splicing", "Disable spliced read alignment.");
$parser->addFlag("abra", "Enable indel realignment with ABRA.");

$parser->addString("out_folder", "Folder where analysis results should be stored. Default is same as in '-folder' (e.g. Sample_xyz/).", true, "default");
$parser->addInt("threads", "The maximum number of threads to use.", true, 4);

extract($parser->parse($argv));

//resolve out_folder
if($out_folder=="default")
{
	$out_folder = $folder;
}

// create logfile in output folder if no filepath is provided:
	if ($parser->getLogFile() == "") $parser->setLogFile($out_folder."/analyze_rna_".date("YmdHis").".log");

//check steps
$steps = explode(",", $steps);
foreach($steps as $step)
{
	if (!in_array($step, $steps_all))
	{
		trigger_error("Unknown processing step '$step'!", E_USER_ERROR);
	}
}

//log server name
list($server) = exec2("hostname -f");
$user = exec('whoami');
$parser->log("Executed on server: ".implode(" ", $server)." as ".$user);

//init
$prefix = $out_folder."/".$name;
$sys = load_system($system, $name);
$build = $sys['build'];
$target_file = $sys['target_file'];

//determine genome from build
$genome = get_path("data_folder")."/genomes/STAR/{$build}/";

//determine gtf from build
$gtfFile = get_path("data_folder")."/dbs/gene_annotations/{$build}.gtf";

//find FASTQ files
$in_for = glob($folder."/*_R1_001.fastq.gz");
$in_rev = glob($folder."/*_R2_001.fastq.gz");

if ((count($in_for) == 0) && (count($in_rev) == 0))
{
	trigger_error("No FASTQ files found!", E_USER_ERROR);
}

$paired = (count($in_rev) != 0);

//mapping and QC
$final_bam = $prefix.".bam";
$qc_fastq = $prefix."_stats_fastq.qcML";
$qc_map = $prefix."_stats_map.qcML";
if (in_array("ma", $steps) || in_array("fu", $steps))
{
	//check FASTQ quality encoding
	$files = array_merge($in_for, $in_rev);
	foreach ($files as $file)
	{
		list($stdout, $stderr) = $parser->exec(get_path("ngs-bits")."FastqFormat", "-in {$file}", true);
		if (!contains($stdout[2], "Sanger"))
		{
			trigger_error("Input file '{$file}' is not in Sanger/Illumina 1.8 format!", E_USER_ERROR);
		}
	}

	//check that adapters are specified
	if ($sys["adapter1_p5"]=="" || $sys["adapter2_p7"]=="")
	{
		trigger_error("No forward and/or reverse adapter sequence given!\nForward: ".$sys["adapter1_p5"]."\nReverse: ".$sys["adapter2_p7"], E_USER_ERROR);
	}

	//adapter trimming + QC (SeqPurge for paired-end, ReadQC+skewer for single-end)
	if ($paired)
	{
		$fastq_trimmed1 = $parser->tempFile("_trimmed.fastq.gz");
		$fastq_trimmed2 = $parser->tempFile("_trimmed.fastq.gz");
		$seqpurge_params = array(
			"-in1", implode(" ", $in_for),
			"-in2", implode(" ", $in_rev),
			"-out1 {$fastq_trimmed1}",
			"-out2 {$fastq_trimmed2}",
			"-a1", $sys["adapter1_p5"],
			"-a2", $sys["adapter2_p7"],
			"-qc", $qc_fastq,
			"-threads", bound($threads, 1, 6),
			"-qcut 0"
			);
		$parser->exec(get_path("ngs-bits")."SeqPurge", implode(" ", $seqpurge_params), true);
	}
	else
	{
		$parser->exec(get_path("ngs-bits")."ReadQC", "-in1 ".implode(" ", $in_for)." -out $qc_fastq", true);

		$fastq_trimmed_tmpdir = $parser->tempFolder();
		$fastq_trimmed1 = "{$fastq_trimmed_tmpdir}/{$name}-trimmed.fastq.gz";
		
		$skewer_params = [
			"-x", $sys["adapter1_p5"],
			"--mode tail",
			"--threads", bound($threads, 1, 6),
			"--compress",
			"--output", "{$fastq_trimmed_tmpdir}/{$name}",
			"-"
		];
		
		$pipeline = [
			["zcat", implode(" ", $in_for)],
			[get_path("skewer"), implode(" ", $skewer_params)]
		];
		$parser->execPipeline($pipeline, "skewer");
	}
}
if (in_array("ma", $steps))
{
	//mapping
	$args = array(
		"-out", $final_bam,
		"-threads", $threads,
		"-in1", $fastq_trimmed1,
		"-genome", $genome,
		"--log", $parser->getLogFile()
	);

	if ($paired) $args[] = "-in2 $fastq_trimmed2";
	if ($no_splicing) $args[] = "-no_splicing";
	
	$parser->execTool("NGS/mapping_star.php", implode(" ", $args));

	//indel realignment
	if ($abra)
	{
		$junction_file = "{$prefix}_splicing.tsv";
		if (!file_exists($junction_file))
		{
			trigger_error("Could not open junction file '$junction_file' needed for indel realignment. Please re-run mapping step.", E_USER_ERROR);
		}
		
		$abra_out = $parser->tempFile("_abra_realigned.bam");
		
		$abra_params = array(
			"-in", $final_bam,
			"-out", $abra_out,
			"-threads", $threads,
			"-build", $build,
			"-gtf", $gtfFile,
			"-junctions", $junction_file
		);
		if (!$paired) $abra_params[] = "-se";
		if (isset($target_file) && $target_file != "") $abra_params[] = "-roi {$target_file}";
		
		$parser->execTool("NGS/indel_realign_abra.php", implode(" ", $abra_params));
		$parser->moveFile($abra_out, $final_bam);
		$parser->indexBam($final_bam, $threads);
	}

	//mapping QC
	$mappingqc_params = array(
		"-in", $final_bam,
		"-out", $qc_map
	);

	$mappingqc_params[] = (isset($target_file) && $target_file != "") ? "-roi {$target_file}" : "-rna";
	if ($build!="GRCh37") $mappingqc_params[] = "-no_cont";

	$parser->exec(get_path("ngs-bits")."MappingQC", implode(" ", $mappingqc_params), true);
}

//read counting
$counts_raw = $prefix."_counts_raw.tsv";
$counts_exon_raw = $prefix."_counts_exon_raw.tsv";
$counts_normalized = $prefix."_counts.tsv";
$counts_qc = $prefix."_stats_rc.tsv";
$repair_bam = $final_bam;
if (in_array("rc", $steps))
{
	if ($paired)
	{
		$tmpdir = $parser->tempFolder();
		$repair_bam = "{$tmpdir}/{$name}.bam";
		$parser->exec(dirname(get_path("feature_counts")) . "/utilities/repair",
			"-i {$final_bam} -o {$repair_bam} -T {$threads}", true);
	}
	
	$args_common = array(
		"-in", $repair_bam,
		"-library_type", $library_type,
		"-gtf_file", $gtfFile,
		"-threads", $threads
	);
	
	if (!$paired)
	{
		$args_common[] = "-single_end";
	}

	$args = array_merge($args_common, [
		"-out", $counts_raw,
		"-qc_file", $counts_qc
	]);

	$parser->execTool("NGS/rc_featurecounts.php", implode(" ", $args));

	// exon-level counting
	$args_exon = array_merge($args_common, [
		"-exon_level",
		"-out", $counts_exon_raw
	]);
	$parser->execTool("NGS/rc_featurecounts.php", implode(" ", $args_exon));

	// read count normalization
	$parser->execTool("NGS/rc_normalize.php", "-in $counts_raw -out $counts_normalized");

	// re-run read counting without duplicate alignments
	$counts_nodup = $prefix."_counts_nodup_raw.tsv";
	$args_dup = array_merge($args_common, [
		"-ignore_dup",
		"-out", $counts_nodup
	]);
	$parser->execTool("NGS/rc_featurecounts.php", implode(" ", $args_dup));
}

//annotate
if (in_array("an", $steps))
{
	$parser->execTool("NGS/rc_annotate.php", "-in $counts_normalized -out $counts_normalized -gtfFile $gtfFile -annotationIds gene_name,gene_biotype");
}

//detect fusions
if (in_array("fu",$steps))
{
	//path of STAR-Fusion index
	$fusion_index = get_path("data_folder")."/genomes/STAR-Fusion/{$build}";
	
	if (is_dir($fusion_index)) {
	
		//add samtools to path
		putenv("PATH=" . implode(":", [
			dirname(get_path("STAR-Fusion_samtools")),
			dirname(get_path("STAR")),
			getenv("PATH")
		]));

		$fusion_tmp_folder = $parser->tempFolder();
		$chimeric_file = "{$prefix}_chimeric.tsv";
		if (!file_exists($chimeric_file))
		{
			trigger_error("Could not open chimeric file '$chimeric_file' needed for STAR-Fusion. Please re-run mapping step.", E_USER_ERROR);
		}

		//remove header from chimeric file for STAR-Fusion
		$chimeric_file_tmp = $parser->tempFile("_STAR_chimeric.tsv");
		$chimeric_file_f = file($chimeric_file);
		array_shift($chimeric_file_f);
		file_put_contents($chimeric_file_tmp, $chimeric_file_f);

		$starfusion_params = [
			"--genome_lib_dir", $fusion_index,
			"--chimeric_junction", $chimeric_file_tmp,
			"--output_dir", $fusion_tmp_folder,
			"--examine_coding_effect",
			"--CPU", $threads,
			"--min_FFPM", "0"
		];

		$input_reads_available = isset($fastq_trimmed1) &&
			is_file($fastq_trimmed1) &&
			(!$paired || (isset($fastq_trimmed2) && is_file($fastq_trimmed2)));
		if ($input_reads_available)
		{
			$starfusion_params[] = "--FusionInspector inspect";
			$starfusion_params[] = "--left_fq ". $fastq_trimmed1;
			if ($paired)
			{
				$starfusion_params[] = "--right_fq " . $fastq_trimmed2;
			}
		}
		//specify neccessary pythonpath
		putenv("PYTHONPATH=" . get_path("STAR-Fusion_pythonpath"));
		$parser->exec(get_path("STAR-Fusion"), implode(" ", $starfusion_params), true);

		$output_files = [ "{$fusion_tmp_folder}/star-fusion.fusion_predictions.abridged.coding_effect.tsv" => "{$prefix}_var_fusions.tsv" ];
		if ($input_reads_available)
		{
			$output_files = array_filter(array_merge([
				"{$fusion_tmp_folder}/FusionInspector-inspect/finspector.consolidated.cSorted.bam" => "{$prefix}_var_fusions.bam",
				"{$fusion_tmp_folder}/FusionInspector-inspect/finspector.consolidated.cSorted.bam.bai" => "{$prefix}_var_fusions.bam.bai",
				"{$fusion_tmp_folder}/FusionInspector-inspect/finspector.fa" => "{$prefix}_var_fusions.fa",
				"{$fusion_tmp_folder}/FusionInspector-inspect/finspector.fa.fai" => "{$prefix}_var_fusions.fa.fai",
				"{$fusion_tmp_folder}/FusionInspector-inspect/finspector.gtf" => "{$prefix}_var_fusions.gtf",
				"{$fusion_tmp_folder}/FusionInspector-inspect/finspector.fusion_inspector_web.html" => "{$prefix}_fusions.html"
			], $output_files), "file_exists", ARRAY_FILTER_USE_KEY);
		}
		foreach ($output_files as $src => $dest)
		{
			$parser->moveFile($src, $dest);
		}
		if ($input_reads_available)
		{
			$igv_session_file = "{$prefix}_var_fusions.xml";
			$igv_tracks = array_filter([
					"{$prefix}_var_fusions.bam",
					"{$prefix}_var_fusions.gtf"
				], "file_exists");
			if (count($igv_tracks) > 0)
			{
				$igv_tracks_arg = implode(" ", $igv_tracks);
				$genome_rel = relative_path(dirname($igv_session_file), "{$prefix}_var_fusions.fa");
				$parser->execTool("NGS/igv_session.php", "-genome {$genome_rel} -out {$igv_session_file} -in {$igv_tracks_arg} -relative");
			}

		}
	}
	else
	{
		trigger_error("STAR-Fusion index not present. Skipping fusion detection step.", E_USER_WARNING);
	}
}

//import to database
if (in_array("db", $steps))
{
	$parser->execTool("NGS/db_check_gender.php", "-in $final_bam -pid $name --log ".$parser->getLogFile());

	//import QC data
	$parser->execTool("NGS/db_import_qc.php", "-id $name -files $qc_fastq $qc_map -force --log ".$parser->getLogFile());
}

?>
