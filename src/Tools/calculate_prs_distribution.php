<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("calculate_prs_distribution", "Calculates the PRS distribution for all diagnostic genome samples.");
$parser->addInfileArray("in", "List of PRS VCFs", false);
$parser->addOutfile("out",  "Output TSV file containing the percentiles for each PRS.", false);
$parser->addOutfile("out2",  "Optional TSV file containing the scores for each Sample.", true, "");
$parser->addString("exclude_disease_group", "Name of a disease group which should be excluded from calculation", true, "");
$parser->addString("processing_system", "Processing system short name to which the cohort should be limited", true, "");
$parser->addInfile("custom_sample_table", "Custom table which should be used for cohort. Must contain columns 'name' and 'path', does no filtering of the table.", true, "");
extract($parser->parse($argv));




// init
$write_sample_output = isset($out2) && ($out2 != "");
$ngs_bits_path = get_path("ngs-bits");

//check excluded disease group
$valid_disease_groups = array("n/a", "Neoplasms", "Diseases of the blood or blood-forming organs", "Diseases of the immune system", "Endocrine, nutritional or metabolic diseases", "Mental, behavioural or neurodevelopmental disorders", 
							"Sleep-wake disorders", "Diseases of the nervous system", "Diseases of the visual system", "Diseases of the ear or mastoid process", "Diseases of the circulatory system", "Diseases of the respiratory system", 
							"Diseases of the digestive system", "Diseases of the skin", "Diseases of the musculoskeletal system or connective tissue", "Diseases of the genitourinary system", "Developmental anomalies", "Other diseases");
if ($exclude_disease_group != "" && !in_array($exclude_disease_group, $valid_disease_groups)) trigger_error("Invalid disease group '".$exclude_disease_group."' given!", E_USER_ERROR);

if ($custom_sample_table == "")
{
	// get all (matching) diagnostic WGS samples
	$export_table = $parser->tempFile("_diag_wgs.tsv");
	
	$args = array();
	$args[] = "-no_bad_samples";
	$args[] = "-no_tumor";
	$args[] = "-no_ffpe";
	$args[] = "-run_finished";
	$args[] = "-no_bad_runs";
	$args[] = "-system_type WGS";
	$args[] = "-ancestry EUR";
	$args[] = "-project_type diagnostic";
	$args[] = "-add_path SAMPLE_FOLDER";
	if ($processing_system != "") $args[] = "-system {$processing_system}";
	$args[] = "-out {$export_table}";
	$parser->exec($ngs_bits_path."NGSDExportSamples", implode(" ", $args));
	//exclude specific disease group

	if ($exclude_disease_group != "") 
	{
		$tmp_table = $parser->tempFile("_diag_wgs.tsv");
		$pipeline[] = $parser->exec($ngs_bits_path."TsvFilter", "-v -filter 'disease_group is {$exclude_disease_group}' -in {$export_table} -out {$tmp_table}");
		$export_table = $tmp_table;
	}
	
	$sample_sheet = Matrix::fromTSV($export_table);
}
else
{
	//use custom table
	$sample_sheet = Matrix::fromTSV($custom_sample_table);
}



$name_idx = $sample_sheet->getColumnIndex("name");
$path_idx = $sample_sheet->getColumnIndex("path");

$temp_out = $parser->tempFile("_prs.tsv");

$results = array();

// iterate over all samples and compute PRS
for ($row_idx=0; $row_idx < $sample_sheet->rows(); $row_idx++) 
{
	$name = $sample_sheet->get($row_idx, $name_idx);
	$source_path = $sample_sheet->get($row_idx, $path_idx);

	print "Calculating PRS for Sample $name (".($row_idx +1)."/".$sample_sheet->rows().")...\n";

	// run VcfCalculatePRS
	$sample_vcf = "$source_path/{$name}_var.vcf.gz";
	if (!file_exists($sample_vcf)) 
	{
		print "\tmissing VCF, skipping sample!\n";
		continue;
	}
	$sample_bam = "$source_path/{$name}.bam";
	if (!file_exists($sample_bam)) 
	{
		$sample_bam = "$source_path/{$name}.cram";
		if (!file_exists($sample_bam)) 
		{
			print "\tmissing BAM/CRAM, skipping sample!\n";
			continue;
		}
	}

	list($stdout, $stderr, $exitcode) = $parser->exec($ngs_bits_path."VcfCalculatePRS", "-in $sample_vcf -bam $sample_bam -out $temp_out -prs ".implode(" ", $in));

	foreach($stdout as $line)
	{
		$pgs_id = trim(explode(":", $line)[0]);
		$prs_score = explode("=", explode(" ", explode(":", $line)[1])[2])[1];
		$results[$pgs_id][$name] = $prs_score;
	}

}

// write results to file
$out_h = fopen2($out, "w");

// write header
$range = range(0.01, 1.0, 0.01);
fwrite($out_h, "#pgs_id\tsample_count\t".implode("\t", $range)."\n");

if ($write_sample_output)
{
	$out2_h = fopen2($out2, "w");
	fwrite($out2_h, "#pgs_id\tsample\tscore\n");
}

foreach ($results as $pgs_id => $prs_values) 
{
	// write scores for each sample to file
	if ($write_sample_output)
	{
		foreach ($prs_values as $sample => $score) 
		{
			fwrite($out2_h, "$pgs_id\t$sample\t$score\n");
		}
	}
	
	$prs_scores = array_values($prs_values);
	if (!sort($prs_scores, SORT_NUMERIC))
	{
		trigger_error("Error sorting PRS scores for '$pgs_id'!", E_USER_ERROR);
	}

	$percentiles = array();
	$sample_count = count($prs_scores);
	print "Number of samples: ".count($prs_scores)."\n";
	foreach ($range as $percentile) 
	{ 
		$float_idx = $percentile*$sample_count;
		$lower_idx = floor($float_idx);
		$upper_idx = ceil($float_idx);
		$ratio = $float_idx - $lower_idx;
		$percentiles[] = ($prs_scores[$lower_idx-1] * $ratio) + ($prs_scores[$upper_idx-1] * (1 - $ratio));
	}

	fwrite($out_h, "$pgs_id\t$sample_count\t".implode("\t", $percentiles)."\n");
}

fclose($out_h);
if ($write_sample_output) fclose($out2_h);

?>