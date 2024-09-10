<?php

/*
	@page indel_realign_abra
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("indel_realign_abra", "Perform indel realignment using ABRA2.");
$parser->addInfileArray("in",  "Input BAM file(s).", false);
$parser->addOutfileArray("out",  "Output BAM file(s). Must be different from input BAM file(s). No index BAI files are created!", false);
//optional
$parser->addString("build", "The genome build to use.", true, "GRCh38");
$parser->addInt("threads", "Maximum number of threads used.", true, 1);
$parser->addInfile("roi", "Target region for realignment.", true, "");
$parser->addFloat("mer",  "ABRA2 minimum edge pruning ratio parameter. Default value is for germline - use 0.02 for somatic data.", true, 0.1);
$parser->addFloat("mad",  "ABRA2 downsampling depth parameter. Default value is for germline - use 5000 for somatic data.", true, 250);
$parser->addFlag("skip_kmer", "Skips k-mer precalcuation (for benchmarking only).");
$parser->addFlag("se", "RNA: Single-end input");
$parser->addInfile("gtf",  "RNA: GTF annotation file.", true, "");
$parser->addInfile("junctions",  "RNA: Junctions output file from STAR mapping.", true, "");
extract($parser->parse($argv));

//init
$ref_genome = genome_fasta($build);
if (count($in)!=count($out))
{
	trigger_error("The number of input and output files must match!", E_USER_ERROR);
}

//create k-mer folder
if (!$skip_kmer && isset($roi))
{
	$kmer_folder = get_path("data_folder")."/dbs/ABRA/";
	if (!file_exists($kmer_folder) && !mkdir($kmer_folder, 0777, true))
	{
		trigger_error("Could not create ABRA2 database folder '$kmer_folder'!", E_USER_ERROR);
	}

	//pre-calcalculate target region file with k-mer size (if not present)
	$kmer_file = $kmer_folder.basename($roi, ".bed")."_".md5_file($roi).".bed";
	if (!file_exists($kmer_file))
	{
		$read_length = 100;
		$kmer_tmp = $parser->tempFile();
		$parser->execSingularity("abra2", get_path("container_abra2"), "java -Xmx16G -cp /opt/abra2-2.23/abra2.jar abra.KmerSizeEvaluator", "$read_length $ref_genome $kmer_tmp $threads $roi", [$ref_genome, $roi]);
		$parser->moveFile($kmer_tmp, $kmer_file);
		chmod($kmer_file, 0777);
	}
}

//indel realignment with ABRA
$params = array();
$in_files = array();
$out_files = array();
$params[] = "--in ".implode(",", $in);
$in_files = array_merge($in_files, $in);
$params[] = "--out ".implode(",", $out);
$out_files = array_merge($out_files, $out);
//manta fix
$params[] = "--no-edge-ci";
if (!$skip_kmer && isset($roi))
{
	$params[] = "--target-kmers ".$kmer_file;
	$in_files[] = $kmer_file;
}
$params[] = "--threads ".$threads;
$params[] = "--tmpdir ".$parser->tempFolder("abra");
$params[] = "--ref {$ref_genome}";
$in_files[] = $ref_genome;
$params[] = "--mer ".$mer;
$params[] = "--mad ".$mad;
if ($se) {
	$params[] = "--single";
}
if (isset($gtf)) {
	$params[] = "--gtf ".$gtf;
	$in_files[] = $gtf;
}
if (isset($junctions)) {
	$params[] = "--junctions ".$junctions;
	$in_files[] = $junctions;
}
$parser->execSingularity("abra2", get_path("container_abra2"), "java -Xmx16G -jar /opt/abra2-2.23/abra2.jar", implode(" ", $params), $in_files, $out_files);

?>