<?php

/**
	@page multisample
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

function extract_info($format, $data)
{
	if ($data==".")
	{
		return array('NONE', 0, 0);
	}
	
	$data = explode(":", $data);
	$data = array_combine($format, $data);
	$depth = array_sum(explode(",",$data['DP'])); 
	$ao = array_sum(explode(",",$data['AO']));
	$genotype = vcfgeno2human($data['GT']);
	return array($genotype, $depth, number_format($ao/$depth,2));
}


//parse command line arguments
$parser = new ToolBase("multisample", "Multi-sample analysis pipeline.");
$parser->addInfileArray("bams", "Input BAM files.", false);
$parser->addStringArray("status", "List of affected status of the input samples (BAMs) - can be 'affected' or 'control'.", false);
$parser->addString("out_folder", "Output folder name.", false);
//optional
$parser->addInfile("system",  "Processing system INI file used for all samples (created from NGSD via first argument of 'bams' by default).", true);
$steps = array("vc", "an");
$parser->addEnum("start", "Start processing at step.", true, $steps, "vc");
extract($parser->parse($argv));

//check input counts
if (count($bams)!=count($status))
{
	trigger_error("Different number of arguments for 'bams' and 'status' parameters! BAMs: ".count($bams)." Status: ".count($status)."!", E_USER_ERROR);
}

//check input status
foreach($status as $stat)
{
	$valid = array("affected", "control");
	if (!in_array($stat, $valid))
	{
		trigger_error("Invalid status '$stat' given for '$bam'. Valid are '".implode("', '", $valid)."'", E_USER_ERROR);
	}
}
$status = array_combine($bams, $status);

//check input sample names
$names = array();
$tmp = array();
foreach($bams as $bam)
{
	$name = basename($bam, ".bam");
	if (isset($tmp[$name]))
	{
		trigger_error("Sample file name '$name' occurs twice in input file list. Each sample must be uniquely indentifyable by name!", E_USER_ERROR);
	}
	$names[$bam] = $name;
	$tmp[$name] = true;
}

//extract processing system information from DB
$sys = load_system($system, $names[$bams[0]]);
$target_file = $sys['type']!="WGS" ? $sys['target_file'] : get_path("data_folder")."/enrichment/ssHAEv6_2017_01_05.bed";
if ($target_file=="")
{
	trigger_error("Cannot perform multi-sample analysis without target region (processing systems of ".$bams[0]." is '".$sys["name_short"]."')!", E_USER_ERROR);
}

//create output folder if missing
$out_folder .= "/";
if (!file_exists($out_folder)) mkdir($out_folder);

//(1) variant calling of all samples together (with very conservative parameters)
$vcf_all = $out_folder."all.vcf.gz";
if ($start=="vc")
{
	$parser->execTool("NGS/vc_freebayes.php", "-bam ".implode(" ", $bams)." -out $vcf_all -target $target_file -min_mq 20 -min_af 0.1 -build ".$sys['build'], true);	
}

//(2) Convert VCF to single-sample format
$indices = array();
$h1 = gzopen($vcf_all, "r");
if ($h1===FALSE) trigger_error("Could not open file '" + $vcf_all + "'.", E_USER_ERROR);
$vcf = $parser->tempFile("_unzipped.vcf");
$h2 = fopen($vcf, "w");
if ($h2===FALSE) trigger_error("Could not open file '" + $vcf + "'.", E_USER_ERROR);
while(!gzeof($h1))
{
	$line = trim(gzgets($h1));
	if (strlen($line)==0) continue;
	
	if ($line[0]=="#" && $line[1]=="#") //comments
	{
		fwrite($h2, $line."\n");
	}
	else if ($line[0]=="#") //header
	{
		//add multi-sample comments
		fwrite($h2, "##FORMAT=<ID=MULTI,Number=.,Type=String,Description=\"Multi-sample genotype information (genotype, depth).\">\n");		
		fwrite($h2, "##ANALYSISTYPE=GERMLINE_MULTISAMPLE\n");
		foreach($bams as $bam)
		{
			fwrite($h2, gsvar_sample_header($names[$bam], array("Status"=>$status[$bam])));
		}
		
		//determine indices for each sample	
		$parts = explode("\t", $line);
		foreach($bams as $bam)
		{
			$indices[$bam] = vcf_column_index($names[$bam], $parts); 
		}
		
		//write main header
		fwrite($h2, implode("\t", array_slice($parts, 0, 9))."\tmulti\n");
	}
	else //content
	{
		$parts = explode("\t", $line);
		$format = explode(":", $parts[8]);
		$muti_info = array();
		foreach($bams as $bam)
		{
			
			$index = $indices[$bam];
			list($gt, $dp) = extract_info($format, $parts[$index]);
			$muti_info[] = $names[$bam]."=$gt|$dp";
		}
		
		//update format field and remove mother/father
		$parts[8] = "MULTI";
		fwrite($h2, implode("\t", array_slice($parts, 0, 9))."\t".implode(",", $muti_info)."\n");
	}
}
fclose($h1);
fclose($h2);

//(3) zip variant list
$vcf_zipped = $out_folder."multi_var.vcf.gz";
$parser->exec("bgzip", "-c $vcf > $vcf_zipped", false); //no output logging, because Toolbase::extractVersion() does not return
$parser->exec("tabix", "-p vcf $vcf_zipped", false); //no output logging, because Toolbase::extractVersion() does not return

//(4) basic annotation
$parser->execTool("Pipelines/annotate.php", "-out_name multi -out_folder $out_folder -system $system -multi");

?>
