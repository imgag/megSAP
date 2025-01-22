<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("an_scarHRD", "Run scarHRD on the CNVs of a given sample. Prefilters the CNVs.");
$parser->addInfile("cnvs", "ClinCNV result file with the called cnvs.", false);
$parser->addString("tumor", "Name of the tumor sample.", false);
$parser->addString("normal", "Name of the normal sample.", false);
$parser->addString("out_folder", "Folder in which the hrd score file will be saved.", false);
$parser->addFloat("min_clone", "Minimum clonality of a CNV to be included in scarHRD analysis. (Prefilter)", true, 0.2);
$parser->addFloat("min_sll", "Minimum scaled loglikelihood of a CNV to be included in scarHRD analysis. (Prefilter)", true, 0.75);
$parser->addFlag("skip_prefilter", "Skip prefiltering of the CNVs", true);
$parser->addFlag("report_filtered", "Calculate the HRD score based on the filtered CNVs. Removes the CNVs flagged as artifacts in the somatic report config in NGSD instead of prefiltering.", true);
$parser->addString("db", "Database to be used for the filtering of CNVs (only used when -filtered is given)", true, "NGSD");
extract($parser->parse($argv));

$out_folder = realpath($out_folder); //R cannot deal with relative paths...

function prepare_clincnv($clincnv, $sample)
{
	$tmp = temp_file(".tsv", "scarhrd_in_");
	$lines = array();
	$lines[] = "SampleID\tChromosome\tStart_position\tEnd_position\ttotal_cn\tA_cn\tB_cn\tploidy";
	
	$sample_id = basename($clincnv, "_clincnv.tsv");
	
	$ploidy = -1;
	$count = 0;
	foreach(file($clincnv) as $line)
	{
		$line = trim($line);
		if ($line[0] == '#')
		{
			if (starts_with($line, "##ploidy:"))
			{
				$ploidy = floatval(explode(" ", $line)[1]);
			}
			continue;
		}
		$parts = explode("\t", $line);
		
		$lines[] = "{$sample}\t{$parts[0]}\t{$parts[1]}\t{$parts[2]}\t{$parts[5]}\t{$parts[7]}\t{$parts[8]}\t{$ploidy}";
		$count++;
	}

	file_put_contents($tmp, implode("\n", $lines)."\n");
	return [$count, $tmp];
}

function generate_report_filtered_cnvs($tumor, $clincnvs, $db)
{
	global $parser;
	
	$db = DB::getInstance($db);
	
	//get somatic report config for the cnvs
	$db_cnv_result = $db->executeQuery("SELECT CONCAT(s.name, \"_\", LPAD(ps.process_id, 2, \"0\")) as ps_tumor, src_cnv.somatic_cnv_id, src_cnv.exclude_artefact, cnv.chr, cnv.start, cnv.end" .
	" FROM sample s" .
		" LEFT JOIN processed_sample ps on s.id = ps.sample_id" .
		" LEFT JOIN somatic_report_configuration src on ps.id = src.ps_tumor_id" .
		" LEFT JOIN somatic_report_configuration_cnv src_cnv on src.id = src_cnv.somatic_report_configuration_id" .
		" LEFT JOIN somatic_cnv cnv on cnv.id = src_cnv.somatic_cnv_id" .
	" WHERE CONCAT(s.name, \"_\", LPAD(ps.process_id, 2, \"0\")) = \"{$tumor}\" and src.id is not null");

	if (count($db_cnv_result) == 0)
	{
		return $clincnvs;
	}
	
	$cnv_exclude_key = array();
	foreach($db_cnv_result as $res_cnv)
	{
		$chr = $res_cnv["chr"];
		$start = $res_cnv["start"];
		$end = $res_cnv["end"];
		$cnv_exclude_key["{$chr}:{$start}-{$end}"] = boolval($res_cnv["exclude_artefact"]);
	}
	
	//filter cnvs that were marked as artifacts out of the cnvs called by clincnv
	$filtered_clincnv = array();
	$cnvs_passed = 0;
	$cnvs_total = 0;
	foreach(file($clincnvs) as $line)
	{
		$line = trim($line);
		
		if ($line[0] == "#")
		{
			$filtered_clincnv[] = $line;
			continue;
		}
		$cnvs_total++;
		$parts = explode("\t", $line);
		
		$chr = $parts[0];
		$start = $parts[1];
		$end = $parts[2];
		
		if (array_key_exists("{$chr}:{$start}-{$end}", $cnv_exclude_key) and $cnv_exclude_key["{$chr}:{$start}-{$end}"])
		{
			continue;	
		}
		
		$cnvs_passed++;
		$filtered_clincnv[] = $line;
	}
	
	// echo $sample->to_str(). "\t has $cnvs_passed\t/$cnvs_total passed cnvs.\n";
	$tmp = temp_file(".tsv", "clincnv_filtered_");
	file_put_contents($tmp, implode("\n", $filtered_clincnv));
	
	return $tmp;
}

function generate_prefiltered_cnvs($tumor, $clincnvs)
{
	global $parser, $min_clone, $min_sll;

	$headers = array();
	$header = "";
	$prefiltered_cnvs = array();
	$cnvs_passed = 0;
	$cnvs_total = 0;
	foreach(file($clincnvs) as $line)
	{
		if (starts_with($line, "##"))
		{
			$headers[] = $line;
			continue;
		};
		if (starts_with($line, "#chr")) 
		{
			$header = $line;
			continue;
		}
		$cnvs_total++;
		
		$header_parts = explode("\t", $header);
		$header_parts[count($header_parts)-1] = trim($header_parts[count($header_parts)-1]);
		$cnv_parts = explode("\t", $line);
		$cnv_parts[count($cnv_parts)-1] = trim($cnv_parts[count($cnv_parts)-1]);
		
		$tumor_clonality = $cnv_parts[array_search("tumor_clonality", $header_parts)];
		$tumor_CN_change = $cnv_parts[array_search("tumor_CN_change", $header_parts)];
		$number_of_regions = $cnv_parts[array_search("number_of_regions", $header_parts)];
		$loglikelihood = $cnv_parts[array_search("loglikelihood", $header_parts)];
		
		$parser->log("CNV: $tumor_clonality, ". floatval($loglikelihood) / floatval($number_of_regions));
		
		if (floatval($tumor_clonality) < $min_clone && floatval($loglikelihood) / floatval($number_of_regions) < $min_sll) continue;
		$parser->log("passed");
		$cnvs_passed++;
		$prefiltered_cnvs[] = $line;
	}

	$parser->log("Total: $cnvs_total || passed: $cnvs_passed");

	$headers_text = implode("", $headers).$header;
	$cnvs_prefiltered_text = implode("", $prefiltered_cnvs);
	
	$tmp = temp_file(".tsv", "clincnv_filtered_");
	file_put_contents($tmp, $headers_text.$cnvs_prefiltered_text);
	
	return $tmp;
}

function write_empty_result($out, $prefix)
{
	$no_cnvs_file = array();
	$no_cnvs_file[] = "\"\"\t\"HRD\"\t\"Telomeric AI\"\t\"LST\"\t\"HRD-sum\"";
	$no_cnvs_file[] = "\"{$prefix}\"\t0\t0\t0\t0";
	
	file_put_contents($out, implode("\n", $no_cnvs_file)."\n");	
}

//run scarHRD if there are cnvs else write an "empty" result file with scores=0. (scarHRD throws error if input is empty).
function run_scarHRD($parser, $cnvs, $count, $prefix, $out_folder)
{
	if ($count != 0)
	{	
		$wd = "/opt/scarHRD/R";
		$parser->execApptainer("scarHRD", "cli_scarHRD.R", "-s $cnvs -o $out_folder -w $wd", [$cnvs], [$out_folder]);
	}
	else
	{
		$out = $out_folder."/{$prefix}_HRDresults.txt";
		write_empty_result($out, $prefix);
	}
}


// MAIN: 
$prefix = "{$tumor}-{$normal}";
list($count, $seg) = prepare_clincnv($cnvs, $prefix);

if ($report_filtered)
{
	$filtered_clincnv = generate_report_filtered_cnvs($tumor, $cnvs, $db);
	list($count, $seg_filtered) = prepare_clincnv($filtered_clincnv, $prefix);
	run_scarHRD($parser, $seg_filtered, $count, $prefix, $out_folder);
}
else if ($skip_prefilter)
{
	run_scarHRD($parser, $seg, $count, $prefix, $out_folder);
} 
else
{
	$filtered_clincnv = generate_prefiltered_cnvs($tumor, $cnvs);
	list($count, $seg_filtered) = prepare_clincnv($filtered_clincnv, $prefix);
	run_scarHRD($parser, $seg_filtered, $count, $prefix, $out_folder);
}
	


?>