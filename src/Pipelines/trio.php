<?php

/**
	@page trio
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

function determine_index($name, $parts)
{
	$index = array_search($name, $parts);
	
	if ($index===FALSE || $index==-1)
	{
		trigger_error("Could not determine index of column '$name' in header line: ".implode(" ", $parts), E_USER_ERROR);
	}
	
	return $index;
}

//parse command line arguments
$parser = new ToolBase("trio", "Trio analysis pipeline.");
$parser->addInfile("f", "BAM file of father.", false, true);
$parser->addInfile("m", "BAM file of mother.", false, true);
$parser->addInfile("c", "BAM file of child (index).", false, true);
$parser->addString("out_folder", "Output folder name.", false);
//optional
$parser->addInfile("system",  "Processing system INI file used for all samples (created from NGSD via processed sample 'c' by default).", true);
$steps_all = array("vc", "an", "cn");
$parser->addString("steps", "Comma-separated list of steps to perform:\nvc=variant calling, an=annotation, cn=copy-number analysis.", true, implode(",", $steps_all));
$parser->addFlag("no_check", "Skip gender check of parents and parent-child correlation check (otherwise done before variant calling)");

extract($parser->parse($argv));

//init
$sample_c = basename($c, ".bam");
$sample_f = basename($f, ".bam");
$sample_m = basename($m, ".bam");

//check steps
$steps = explode(",", $steps);
foreach($steps as $step)
{
	if (!in_array($step, $steps_all)) trigger_error("Unknown processing step '$step'!", E_USER_ERROR);
}

//extract processing system information from DB
$sys = load_system($system, basename($c, ".bam"));
$target_file = $sys['type']!="WGS" ? $sys['target_file'] : get_path("data_folder")."/enrichment/ssHAEv6_2017_01_05.bed";
if ($target_file=="")
{
	trigger_error("Cannot perform trio analysis without target region (processing systems of {$c} is '".$sys["name_short"]."')!", E_USER_ERROR);
}

//prepare multi-sample paramters
$args_multisample = [
	"-bams $f $m $c",
	"-status control control affected",
	"-out_folder $out_folder",
	"-system $system",
	"-prefix trio",
	];

//variant calling
if (in_array("vc", $steps))
{
	//check parent-child correlation
	if (!$no_check)
	{
		$min_corr = 0.45;
		$output = $parser->exec(get_path("ngs-bits")."SampleSimilarity", "-in $f $c -mode bam -max_snps 4000", true);
		$correlation = explode("\t", $output[0][1])[3];
		if ($correlation<$min_corr)
		{
			trigger_error("The genotype correlation of father and child is {$correlation}; it should be above {$min_corr}!", E_USER_ERROR);
		}
		$output = $parser->exec(get_path("ngs-bits")."SampleSimilarity", "-in $m $c -mode bam -max_snps 4000", true);
		$correlation = explode("\t", $output[0][1])[3];
		if ($correlation<$min_corr)
		{
			trigger_error("The genotype correlation of mother and child is {$correlation}; it should be above {$min_corr}!", E_USER_ERROR);
		}
	}
	
	//check gender of parents if NGSD is available and check is not disabled
	if (!$no_check && db_is_enabled("NGSD"))
	{
        //check gender of father
        $parser->execTool("NGS/db_check_gender.php", " -in $f -pid $sample_f -gender male");
        $parser->execTool("NGS/db_check_gender.php", " -in $m -pid $sample_m -gender female");
	}

	//variant calling with multi-sample pipeline
	$parser->execTool("Pipelines/multisample.php", implode(" ", $args_multisample)." -steps vc", true); 
	
	//UPD detection
	$upd_file = "$out_folder/trio_upd.tsv";
	$args_upd = [
		"-in {$out_folder}/all.vcf.gz",
		"-c {$sample_c}",
		"-f {$sample_f}",
		"-m {$sample_m}",
		"-out {$upd_file}",
		];
	$c_cnvs = substr($c, 0, -4)."_cnvs.tsv";
	if (file_exists($c_cnvs))
	{
		$args_upd[] = "-exclude {$c_cnvs}";
	}
	else
	{
		trigger_error("Child CNV file not found for UPD detection: {$c_cnvs}", E_USER_WARNING);
	}
	$parser->exec(get_path("ngs-bits")."UpdHunter", implode(" ", $args_upd), true);
}

//annotation
if (in_array("an", $steps))
{
	//annotation with multi-sample pipeline
	$parser->execTool("Pipelines/multisample.php", implode(" ", $args_multisample)." -steps an", true);	
	
	//determine mendelian error rate
	$vars_all = 0;
	$vars_high_depth = 0;
	$vars_mendelian_error = 0;
	$gsvar = "$out_folder/trio.GSvar";
	$h = fopen($gsvar, "r");
	while(!feof($h))
	{
		$line = trim(fgets($h));
		if ($line=="") continue;
		
		//skip comments
		if (starts_with($line, "##")) continue;
		
		//header columns (determine column indices)		
		if (starts_with($line, "#"))
		{
			if ($sample_c=="" || $sample_f=="" || $sample_m=="")
			{
				trigger_error("At least one of the sample names could not be determined! C={$sample_c} F={$sample_f} M={$sample_m}", E_USER_ERROR);
			}
			
			$parts = explode("\t", $line);
			
			$i_c = determine_index($sample_c, $parts);
			$i_f = determine_index($sample_f, $parts);
			$i_m = determine_index($sample_m, $parts);
			$i_quality = determine_index("quality", $parts);
			
			continue;
		}
		
		++$vars_all;
		
		//parse content lines
		$parts = explode("\t", $line);
		$chr = $parts[0];
		$geno_c = $parts[$i_c];
		$geno_f = $parts[$i_f];
		$geno_m = $parts[$i_m];
		
		//count mendelian errors
		if ($chr!="chrX" && $chr!="chrY" && $chr!="chrMT") //only for autosomes
		{
			//only for variants with DP>=20
			$high_depth = true;
			$q_parts = explode(";", $parts[$i_quality]);
			foreach($q_parts as $q_part)
			{
				if (starts_with($q_part, "DP="))
				{
					$dps = explode(",", substr($q_part, 3));
					foreach($dps as $dp)
					{
						if ($dp<20) $high_depth = false;
					}
					break;
				}
			}
			
			if ($high_depth)
			{
				++$vars_high_depth;
				
				//hom, hom => het/wt
				if ($geno_f=="hom" && $geno_m=="hom" && $geno_c!="hom") $vars_mendelian_error += 1;
				//hom, het => wt
				else if ($geno_f=="hom" && $geno_m=="het" && $geno_c=="wt") $vars_mendelian_error += 1;
				else if ($geno_f=="het" && $geno_m=="hom" && $geno_c=="wt") $vars_mendelian_error += 1;
				//het, wt => hom
				else if ($geno_f=="het" && $geno_m=="wt" && $geno_c=="hom") $vars_mendelian_error += 1;
				else if ($geno_f=="wt" && $geno_m=="het" && $geno_c=="hom") $vars_mendelian_error += 1;
				//wt, wt  => het/hom
				else if ($geno_f=="wt" && $geno_m=="wt" && $geno_c!="wt") $vars_mendelian_error += 1;
			}
		}
	}
	print "Overall variants: {$vars_all}\n";
	print "Medelian errors: ".number_format(100.0*$vars_mendelian_error/$vars_high_depth, 2)."% (of {$vars_high_depth} high-depth, autosomal variants)\n";
	
	//determine gender of child
	list($stdout, $stderr) = $parser->exec(get_path("ngs-bits")."SampleGender", "-method hetx -in $c", true);
	$gender_data = explode("\t", $stdout[1])[1];
	if ($gender_data!="male" && $gender_data!="female") $gender_data = "n/a";
	print "Gender of child (from data): {$gender_data}\n";
	
	//write output file
	$tmp = $parser->tempFile(".GSvar");
	$h2 = fopen($tmp, "w");
	rewind($h);
	while(!feof($h))
	{
		$line = trim(fgets($h));
		if ($line=="") continue;
		
		//compare child gender from data with gender from header
		if (starts_with($line, "##SAMPLE=<ID={$sample_c},"))
		{
			//determine gender from header
			$parts = explode(",", substr($line, 10, -1));
			foreach($parts as $part)
			{
				if (starts_with($part, "Gender="))
				{
					list(, $gender_header) = explode("=", $part);
				}
			}
			print "Gender of child (from header): {$gender_data}\n";
			
			//deviating => error
			if ($gender_data!=$gender_header && $gender_header!="n/a" && $gender_data!="n/a")
			{
				trigger_error("Gender of child from sample header '{$gender_header}' deviates from gender from data '{$gender_data}'", E_USER_ERROR);
			}
				
			//replace gender in header by gender from data
			if ($gender_data!=$gender_header)
			{
				$line = strtr($line, array("Gender={$gender_header}"=>"Gender={$gender_data}"));
			}
		}
		
		//update analysis type
		if (starts_with($line, "##ANALYSISTYPE="))
		{
			$line = "##ANALYSISTYPE=GERMLINE_TRIO";
		}
		
		fwrite($h2, "$line\n");
	}
	fclose($h);
	fclose($h2);
	$parser->moveFile($tmp, $gsvar);
}

//copy-number
if (in_array("cn", $steps))
{
	$parser->execTool("Pipelines/multisample.php", implode(" ", $args_multisample)." -steps cn", true);	
}

?>
