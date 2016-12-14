<?php

/**
	@page trio
	
	@todo AO should not be summed up, but the correct allele should be counted only - but wait for fix by vcflib team.
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);


function convert_genotype($gt)
{
	$gt = strtr($gt, "/.", "|0");
	if ($gt=="0|0")
	{
		return "WT";
	}
	else if ($gt=="0|1" || $gt=="1|0")
	{
		return "HET";
	}
	else if ($gt=="1|1")
	{
		return "HOM";
	}
	else
	{
		trigger_error("Invalid genotype '$gt'!", E_USER_ERROR);
	}
}

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
	$genotype = convert_genotype($data['GT']);
	return array($genotype, $depth, number_format($ao/$depth,2));
}


//parse command line arguments
$parser = new ToolBase("trio", "Trio analysis pipeline.");
$parser->addInfile("f", "BAM file of father.", false, true);
$parser->addInfile("m", "BAM file of mother.", false, true);
$parser->addInfile("c", "BAM file of child (index).", false, true);
$parser->addString("out_folder", "Output folder name.", false);
//optional
$parser->addInfile("system",  "Processing system INI file used for all samples (created from NGSD via processed sample 'c' by default).", true);
$steps = array("check", "realign", "vc", "an");
$parser->addEnum("start", "Start processing at step.", true, $steps, "check");
extract($parser->parse($argv));

//extract processing system information from DB
$sys = load_system($system, basename($c, ".bam"));
if ($sys['target_file']=="")
{
	trigger_error("Cannot perform trio analysis without target region (processing systems of child is '".$sys["name_short"]."')!", E_USER_ERROR);
}

//create output folder if missing
$out_folder .= "/";
if (!file_exists($out_folder)) mkdir($out_folder);

//(1) check genders of parents
if ($start=="check")
{
	list($stdout) = $parser->exec(get_path("ngs-bits")."SampleGender", "-method hetx -in $f", true);
	$gender = trim(substr($stdout[count($stdout)-1], 7));
	if (starts_with($gender,"unknown"))
	{
		$parser->log("Could not check gender of father. It could not be determine from BAM file.");
	}
	else if (!starts_with($gender,"male"))
	{
		trigger_error("Gender of father is not male: '$gender'!", E_USER_ERROR);
	}
	list($stdout) = $parser->exec(get_path("ngs-bits")."SampleGender", "-method hetx -in $m", true);
	$gender = trim(substr($stdout[count($stdout)-1], 7));
	if (starts_with($gender,"unknown"))
	{
		$parser->log("Could not check gender of mother. It could not be determine from BAM file.");
	}
	else if (!starts_with($gender,"female"))
	{
		trigger_error("Gender of mother is not female: '$gender'!", E_USER_ERROR);
	}
}

//(2) check parent-child correlation
if ($start=="check")
{
	$output = $parser->exec(get_path("ngs-bits")."SampleCorrelation", "-in1 $f -in2 $c -bam -max_snps 4000", true);
	$parts = explode(":", $output[0][1]);
	if ($parts[1]<0.5)
	{
		trigger_error("The genotype correlation of father and child is ".$parts[1]."; it should be above 0.5!", E_USER_ERROR);
	}
	$output = $parser->exec(get_path("ngs-bits")."SampleCorrelation", "-in1 $m -in2 $c -bam -max_snps 4000", true);
	$parts = explode(":", $output[0][1]);
	if ($parts[1]<0.5)
	{
		trigger_error("The genotype correlation of mother and child is ".$parts[1]."; it should be above 0.5!", E_USER_ERROR);
	}
}

//(4) variant calling of all three samples together (with very conservative parameters)
$vcf_all = $out_folder.basename($c, ".bam")."_all.vcf.gz";
if ($start=="check" || $start=="realign" || $start=="vc")
{
	$parser->execTool("NGS/vc_freebayes.php", "-bam $c $m $f -out $vcf_all -target ".$sys['target_file']." -min_mq 20 -min_af 0.1 -build ".$sys['build'], true);	
}

//(5) convert VCF to file with only child column, but with TRIO annotation
$h1 = gzopen($vcf_all, "r");
if ($h1===FALSE) trigger_error("Could not open file '" + $vcf_all + "'.", E_USER_ERROR);
$vcf = $parser->tempFile("_unzipped.vcf");
$h2 = fopen($vcf, "w");
if ($h2===FALSE) trigger_error("Could not open file '" + $vcf + "'.", E_USER_ERROR);
while(!gzeof($h1))
{
	$line = trim(gzgets($h1));
	if (strlen($line)==0) continue;
	
	if ($line[0]=="#" && $line[1]=="#") //comment
	{
		if (starts_with($line, "##FORMAT=<ID=QA")) //add trio header
		{
			fwrite($h2, $line."\n");
			fwrite($h2, "##FORMAT=<ID=TRIO,Number=.,Type=String,Description=\"Trio information: GT/DP/AF of child, GT/DP/AF of mother, GT/DP/AF of father\">\n");
		}
		else //write out other headers
		{
			fwrite($h2, $line."\n");
		}
	}
	else if ($line[0]=="#") //header
	{
		$parts = explode("\t", $line);
		$parts = array_slice($parts, 0, -2);
		fwrite($h2, implode("\t", $parts)."\n");
	}
	else //content
	{
		$parts = explode("\t", $line);
		$format = explode(":", $parts[8]);
		if ($parts[9]==".") continue;
		$c_info = extract_info($format, $parts[9]);
		$m_info = extract_info($format, $parts[10]);
		$f_info = extract_info($format, $parts[11]);
		
		//update format field and remove mother/father
		$parts[8] = $parts[8].":TRIO";
		$parts = array_slice($parts, 0, -2);
		
		fwrite($h2, implode("\t", $parts).":".implode(",", $c_info).",".implode(",", $m_info).",".implode(",", $f_info)."\n");
	}
}
fclose($h1);
fclose($h2);

//(7) sort variants by genomic position and zip
$tmp = $parser->tempFile(".vcf");
$parser->exec(get_path("ngs-bits")."VcfStreamSort","-in $vcf -out $tmp", true);
$vcf_zipped = $out_folder.basename($c, ".bam")."_var.vcf.gz";
$parser->exec("bgzip", "-c $tmp > $vcf_zipped", false); //no output logging, because Toolbase::extractVersion() does not return
$parser->exec("tabix", "-p vcf $vcf_zipped", false); //no output logging, because Toolbase::extractVersion() does not return

//(8) basic annotation
$parser->execTool("Pipelines/annotate.php", "-out_name ".basename($c, ".bam")." -out_folder $out_folder -system $system");

//(9) add trio annotation column
$gsvar = $out_folder.basename($c, ".bam").".GSvar";
$parser->exec(get_path("ngs-bits")."TrioAnnotation", "-in $gsvar -out $gsvar", true);

?>
