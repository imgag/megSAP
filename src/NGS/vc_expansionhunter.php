<?php

/**
  @page vc_expansionhunter
  
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_expansionhunter", "Call repeat expansions with Expansion Hunter. Creates an VCF file.");
$parser->addInfile("in", "Input BAM file.", false);
$parser->addOutfile("out", "Output VCF file.", false);
//optional
$parser->addInt("threads", "The maximum number of threads used.", true, 2);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
$parser->addString("pid", "Processed sample name (e.g. 'GS120001_01'). If unset BAM file name will be used.", true);
$parser->addString("variant_catalog", "Variant catalog to use.", true, repository_basedir()."/data/repeat_expansions/ExpansionHunter_variant_catalog_grch38.json");
$parser->addFlag("no_images", "Skip creation of SVG images");
extract($parser->parse($argv));

// use BAM file name as fallback if no processed sample name is provided
if(!isset($pid)) $pid = basename($in);

//init
$out_prefix = dirname($out)."/".basename($out, ".vcf");
$expansion_hunter_binary = get_path("expansion_hunter");


//### perform RE calling ###

// prepare command
$args = [];
$args[] = "--threads $threads";
$args[] = "--reads $in";
$args[] = "--reference ".genome_fasta($build);
$args[] = "--output-prefix {$out_prefix}";
$args[] = "--variant-catalog $variant_catalog";
		
// add gender info
if (db_is_enabled("NGSD"))
{
	$db = DB::getInstance("NGSD", false);
	$info = get_processed_sample_info($db, $pid, false);
	if (!is_null($info))
	{
		if ($info['gender'] != 'n/a') $args[] = "--sex ".$info['gender'];
	}
}

// run Expansion Hunter
$parser->exec($expansion_hunter_binary, implode(" ", $args));


// add ExpansionHunter version string to VCF
list($stdout, $stderr, $ec) = $parser->exec($expansion_hunter_binary, "-v");

// parse STDERR
list($date, $tool_version) = explode(",", $stdout[0]);
$tool_version = strtr($tool_version, array("[Starting " => "", "]" => ""));

// read VCF file into memory
$vcf_file_content = file($out);

// write back to disk and insert additional comments
$fh = fopen2($out, 'w');
$header_written = false;
foreach($vcf_file_content as $line)
{
	fwrite($fh, $line);
	// add info after the first line (##fileformat=...)
	if (!$header_written && starts_with($line, "##"))
	{
		// add file/call info to file
		fwrite($fh, "##filedate=$date\n");
		fwrite($fh, "##source={$tool_version}\n");
		fwrite($fh, "##reference=".genome_fasta($build)."\n");
		$header_written = true;
	}
}

//### use REViewer to create SVG for each repeat expansion ###
if (!$no_images)
{
	// sort and index ExpansionHunter BAM
	$parser->sortBam($out_prefix."_realigned.bam", $out_prefix."_realigned.bam", 1);
	$parser->indexBam($out_prefix."_realigned.bam", 1);

	// create output folder
	$svg_folder = dirname($out)."/repeat_expansions/";
	if (!file_exists($svg_folder))
	{
		mkdir($svg_folder);
	}
	else
	{
		// delete previous SVGs
		foreach (glob("{$svg_folder}*.svg") as $file) unlink($file);
	}

	// remove loci which are not covered by sufficient reads (no GT called or no flanking reads)
	$re_usable = array();
	foreach (file($out_prefix.".vcf") as $line) 
	{
		// skip comment/header lines
		if (starts_with($line, "#")) continue;

		// split VCF line
		$vcf_line = explode("\t", $line);

		// extract locus
		$info = explode(";", $vcf_line[7]);
		$rep_id = "";
		foreach ($info as $value) 
		{
			if (starts_with(trim($value), "REPID="))
			{
				$rep_id = trim(explode("=", $value)[1]);
				break;
			}
		}
		if ($rep_id == "") trigger_error("ERROR: VCF line without repeat id!", E_USER_ERROR);

		// extract genotype
		$format = explode(":", $vcf_line[8]);
		$sample = explode(":", $vcf_line[9]);
		
		$re_usable[$rep_id] = ($sample[array_search("GT", $format, true)] != "./.") && ($sample[array_search("ADFL", $format, true)] != "0/0") 
								&& ($sample[array_search("GT", $format, true)] != ".") && ($sample[array_search("ADFL", $format, true)] != "0/0");
	}

	// extract RE names
	$loci = array();
	$variant_catalog_content = json_decode(file_get_contents($variant_catalog), true);
	foreach ($variant_catalog_content as $repeat) 
	{
		if (isset($repeat["LocusId"]))
		{
			$locus = $repeat["LocusId"];
			if (isset($re_usable[$locus]) && $re_usable[$locus]) 
			{
				$loci[] = $locus;
			}
		}
	}

	// create SVG for each RE
	$reviewer_binary = get_path("REViewer");
	$args = [];
	$args[] = "--reads {$out_prefix}_realigned.bam";
	$args[] = "--vcf {$out_prefix}.vcf";
	$args[] = "--reference ".genome_fasta($build);
	$args[] = "--catalog {$variant_catalog}";
	foreach ($loci as $locus) 
	{
		$prefix = $svg_folder.basename($out, ".vcf");
		$parser->exec($reviewer_binary, implode(" ", $args)." --output-prefix $prefix --locus $locus", true, false, true);
		
		//rename to keep the naming consistent with v0.1.1
		$source_file = $prefix.".".$locus.".svg";
		if (file_exists($source_file))
		{
			$parser->moveFile($source_file, $prefix."_".$locus.".svg");
		}
		else
		{
			trigger_error("SVG file for RE locus {$locus} not found: {$source_file}", E_USER_NOTICE);
		}
	}
}

?>