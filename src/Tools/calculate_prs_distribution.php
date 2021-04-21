<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("calculate_prs_distribution", "Calculates the PRS distribution for all diagnostic genome samples.");
$parser->addInfileArray("in", "List of PRS VCFs", false);
$parser->addOutfile("out",  "Output TSV file containing the percentiles for each PRS.", false);
$parser->addOutfile("out2",  "Optional TSV file containing the scores for each Sample.", true, "");
extract($parser->parse($argv));

// init
$write_sample_output = isset($out2) && ($out2 != "");

// get all diagnostic WGS samples
$export_table = $parser->tempFile("_diag_wgs.tsv");
$ngs_bits_path = get_path("ngs-bits");
$pipeline = array();
$pipeline[] = array($ngs_bits_path."NGSDExportSamples", "-no_bad_samples -no_tumor -no_ffpe -run_finished -no_bad_runs -add_path SAMPLE_FOLDER");
$pipeline[] = array($ngs_bits_path."TsvFilter", "-filter 'project_type is diagnostic'");
$pipeline[] = array($ngs_bits_path."TsvFilter", "-filter 'system_type is WGS' -out $export_table");
$parser->execPipeline($pipeline, "Export Samples");

$sample_sheet = Matrix::fromTSV($export_table);


$name_idx = $sample_sheet->getColumnIndex("name");
$path_idx = $sample_sheet->getColumnIndex("path");

$temp_out = $parser->tempFile("_prs.tsv");

$results = array();
foreach ($in as $vcf) 
{
	$results[basename($vcf, ".vcf")] = array();
}

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

	list($stdout, $stderr, $exitcode) = $parser->exec($ngs_bits_path."VcfCalculatePRS", "-in $sample_vcf -out $temp_out -prs ".implode(" ", $in));

	foreach($stdout as $line)
	{
		$split_line = explode(":\t", $line);
		$results[trim($split_line[0])][$name] = trim($split_line[1]);
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