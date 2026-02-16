<?php

/** 
	@page extract_signatures
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("extract_signatures", "Extracting genome mutational signatures from a vcf file.");
$parser->addInfile("in",  "VCF file or folder with multiple VCF files for SNV Signatures or somatic ClinCNV segmentation file for CNV signatures.", false);
$parser->addString("out",  "Output folder.", false);
//optional
$parser->addInt("minSig", "Minimum number of Signatures tested.", true, 1);
$parser->addInt("maxSig", "Maximum number of Signatures tested.", true, 10);
$parser->addInt("nmfRep", "Number of nmf replicates.", true, 50);
$parser->addString("reference", "Reference genome used. Options are: GRCh37, GRCh38, mm9, mm10.", true, "GRCh38");
$parser->addString("mode", "The kind of signatures to be calculated. Options are: snv, cnv", true, "snv");
$parser->addFlag("fullOutput", "copies the full result folder to the given out.");
$parser->addInt("threads", "Maximum number of threads used.", true, 4);
$parser->addString("seeds", "VCF file or folder with multiple VCF files for SNV Signatures or somatic ClinCNV segmentation file for CNV signatures.", true, "random");

extract($parser->parse($argv));

function csv2tsv($csv)
{
	if (!file_exists($csv)) return;
	
	//parse CSV
	$file = file($csv);
	$headers = explode(", ", trim($file[0]));
	if (count($headers)!=7) trigger_error("CNV line does not have 7 comma-separted parts: ".file[0], E_USER_ERROR);
	$parts = explode(", ", trim($file[1]));
	if (count($parts)!=7) trigger_error("CNV line does not have 7 comma-separted parts: ".file[1], E_USER_ERROR);
	
	//write TSV
	$tsv = substr($csv, 0, -4).".tsv";
	$h = fopen2($tsv, "w");
	for ($i=2; $i<count($headers); ++$i)
	{
		fputs($h, "##".$headers[$i].": ".$parts[$i]."\n");
	}
	fputs($h, "#type\tsignature\tvalue\n");
	for ($i=0; $i<2; ++$i)
	{
		$parts2 = explode("&", $parts[$i]);
		foreach($parts2 as $entry)
		{
			$entry = trim($entry);
			if ($entry=="") continue;
			
			$signature = $entry;
			$value = "";
			if (contains($entry, "(") && contains($entry, "%)"))
			{
				$pos = strrpos($entry, "(");
				$signature = trim(substr($entry, 0, $pos));
				$value = trim(substr($entry, $pos+1, -2));
			}
			fputs($h, $headers[$i]."\t{$signature}\t{$value}\n");
		}
	}
	fclose($h);
	
	//delete csv
	unlink($csv);
}

function calculate_cnv_signatures($in, $minSig, $maxSig, $nmfRep, $seeds)
{
	global $parser;
	global $reference;
	global $threads;
	global $data;

	$result_folder = temp_folder();
	$tmp = cnv_prepare_clincnv($in);

	if ($tmp == "")
	{
		trigger_error("Warning: No CNVs called! - Skipping CNV signature calculation.", E_USER_NOTICE);
		return;
	}
	
	$in_files = [];
	if (is_file(($seeds)) || is_dir($seeds))
	{
		$in_files[] = $seeds;
	}
	
	$parser->execApptainer("SigProfilerExtractor", "python3 -c", "'from SigProfilerExtractor import sigpro as sig; sig.sigProfilerExtractor(\"seg:FACETS\", \"$result_folder\", \"$tmp\", reference_genome=\"{$reference}\", minimum_signatures={$minSig}, maximum_signatures={$maxSig}, nmf_replicates={$nmfRep}, cpu={$threads}, seeds=\"{$seeds}\", volume=\"{$data}\")'", $in_files, [$data]);

	copy_cnv_result_files($result_folder);
}

function cnv_prepare_clincnv($cnvFile)
{
	$tmp = temp_file(".seg", "cnv_input");
	$lines = array();
	foreach(file($cnvFile) as $line)
	{
		$line = trim($line);
		if ($line === "" || (isset($line[0]) && $line[0] == '#')) continue;

		$parts = explode("\t", $line);
		$lines[] = "{$parts[0]}\t{$parts[1]}\t{$parts[2]}\t{$parts[5]}\t{$parts[8]}\t".basename($cnvFile, "_clincnv.tsv")."\t0";
	}

	if (empty($lines)) return "";

	$header = "chr\tstart\tend\ttcn.em\tlcn.em\tsample\tid";
	array_unshift($lines, $header);
	file_put_contents($tmp, implode("\n", $lines)."\n");
	return $tmp;
}

function copy_cnv_result_files($result_folder)
{
	global $fullOutput;
	global $out;

	if ($fullOutput)
	{
		exec("cp -r {$result_folder}/* {$out}");
	}
	else
	{
		$files_to_copy = [
			"/JOB_METADATA.txt",
			"/Seeds.txt",
			"/CNV48/Suggested_Solution/COSMIC_CNV48_Decomposed_Solution/CNV_Decomposition_Plots.pdf",
			"/CNV48/Suggested_Solution/COSMIC_CNV48_Decomposed_Solution/De_Novo_map_to_COSMIC_CNV48.csv"
		];

		foreach($files_to_copy as $file)
		{
			$result_file = "{$result_folder}{$file}";
			if (is_file($result_file))
			{
				exec("cp -r {$result_file} {$out}");
			}
		}
	}
	
	csv2tsv($out."/De_Novo_map_to_COSMIC_CNV48.csv");
}

function calculate_snv_signatures($in, $minSig, $maxSig, $nmfRep, $seeds)
{
	global $parser;
	global $reference;
	global $threads;
	global $data;

	$result_folder = temp_folder();
	$in_dir = prepare_input_files($in);

	if (is_file(($seeds)) || is_dir($seeds))
	{
		$in_files = [$seeds];
	}

	$in_files[] = $in_dir;
	$in_files[] = $data;
	
	//run container
	$parser->execApptainer("SigProfilerExtractor", "python3 -c", "'from SigProfilerExtractor import sigpro as sig; sig.sigProfilerExtractor(\"vcf\", \"$result_folder\", \"{$in_dir}\", reference_genome=\"{$reference}\", minimum_signatures={$minSig}, maximum_signatures={$maxSig}, nmf_replicates={$nmfRep}, cpu={$threads}, seeds=\"{$seeds}\", volume=\"{$data}\")'", $in_files, [$data]);

	copy_snv_result_files($result_folder);
}

function prepare_input_files($in)
{
	global $parser;

	if (is_dir($in)) return $in;
	else if (is_file($in))
	{
		$in_dir = temp_folder();
		$parser->copyFile($in, $in_dir."/".basename($in));
		return $in_dir;
	}
}

function copy_snv_result_files($result_folder)
{
	global $fullOutput;
	global $out;

	if ($fullOutput)
	{
		exec("cp -r {$result_folder}/* {$out}");
	}
	else
	{
		$files_to_copy = [
			"/JOB_METADATA.txt",
			"/Seeds.txt",
			"/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/SBS96_Decomposition_Plots.pdf",
			"/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/De_Novo_map_to_COSMIC_SBS96.csv",
			"/ID83/Suggested_Solution/COSMIC_ID83_Decomposed_Solution/ID83_Decomposition_Plots.pdf",
			"/ID83/Suggested_Solution/COSMIC_ID83_Decomposed_Solution/De_Novo_map_to_COSMIC_ID83.csv",
			"/DBS78/Suggested_Solution/COSMIC_DBS78_Decomposed_Solution/De_Novo_map_to_COSMIC_DBS78.csv",
			"/DBS78/Suggested_Solution/COSMIC_DBS78_Decomposed_Solution/DBS78_Decomposition_Plots.pdf"
		];

		foreach($files_to_copy as $file)
		{
			$result_file = "{$result_folder}{$file}";
			if (is_file($result_file))
			{
				exec("cp -r {$result_file} {$out}");
			}
		}
	}
	
	csv2tsv($out."/De_Novo_map_to_COSMIC_SBS96.csv");
	csv2tsv($out."/De_Novo_map_to_COSMIC_ID83.csv");
	csv2tsv($out."/De_Novo_map_to_COSMIC_DBS78.csv");
}

function test_input($in, $mode, $out_folder, $minSig, $maxSig, $nmfRep, $seeds)
{
	if (!file_exists($in)) 
	{
		trigger_error("Input: '$in' is neither a file nor a directory.", E_USER_ERROR);
    }

    if (is_dir($in) && $mode === "cnv") 
	{
		trigger_error("The input cannot be a directory for CNV signature calculation.", E_USER_ERROR);
    }

    if (!is_dir($out_folder)) 
	{
        mkdir($out_folder, 0777, true);
    }

	if ($minSig < 1 || $maxSig < 1)
	{
		trigger_error("Neither the minimum nor the maximum number of signatures is allowed to be smaller than one.", E_USER_ERROR);
	}

	if ($maxSig < $minSig)
	{
		trigger_error("The minimum number of signatures has to be less or equal to the maximum number of signatures.", E_USER_ERROR);
	}

	if ($nmfRep < 1)
	{
		trigger_error("Number of nmf replicates cannot be lower than 1", E_USER_ERROR);
	}

	if ($seeds != "random" && !is_file($seeds))
	{
		trigger_error("The given file for argument --seeds does not exist.", E_USER_ERROR);
	}
}

function reference_installed($data, $reference)
{
	
}

function install_reference($parser, $data, $reference)
{
	global $full_env_prefix;
	
	if (! is_dir($data)) mkdir($data);
	$parser->log("Installing reference '{$reference}' for MatrixGenerator in folder '{$data}'");
	$cmd = $parser->execApptainer("SigProfilerExtractor", "python3 -c", "'from SigProfilerMatrixGenerator.scripts import reference_genome_manager; rgm = reference_genome_manager.ReferenceGenomeManager(\"${data}\"); rgm.download_genome(\"GRCh38\")'", [], [$data], true);
	$parser->exec("{$full_env_prefix} {$cmd}", "");
}


//MAIN
$data = get_path("sigprofilerextractor_data").get_path("container_SigProfilerExtractor")."/";


$env_var1 = "SIGPROFILERMATRIXGENERATOR_VOLUME={$data}";
$env_var2 = "SIGPROFILERASSIGNMENT_VOLUME={$data}";
$prefix = container_platform()=='apptainer' ? "APPTAINERENV" : "SINGULARITYENV";
$full_env_prefix = "{$prefix}_{$env_var1} {$prefix}_{$env_var2}";

if (true || ! reference_installed($data, $reference))
{
	install_reference($parser, $data, $reference);
}



test_input($in, $mode, $out, $minSig, $maxSig, $nmfRep, $seeds);
if ($mode === "snv")
{
	calculate_snv_signatures($in, $minSig, $maxSig, $nmfRep, $seeds);
}
elseif ($mode === "cnv")
{
	calculate_cnv_signatures($in, $minSig, $maxSig, $nmfRep, $seeds);
}

?>