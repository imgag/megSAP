<?php

/**
	@page trio
	
	@todo use three separate columns for genotypes and handle trios similar to multi-sample analyses.
	@todo AO should not be summed up, but the correct allele should be counted only - but wait for fix by vcflib team.
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
	$genotype = vcfgeno2human($data['GT'], true);
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
$steps = array("check", "vc", "an");
$parser->addEnum("start", "Start processing at step.", true, $steps, "check");
extract($parser->parse($argv));

//extract processing system information from DB
$sys = load_system($system, basename($c, ".bam"));
$target_file = $sys['type']!="WGS" ? $sys['target_file'] : get_path("data_folder")."/enrichment/ssHAEv6_2017_01_05.bed";
if ($target_file=="")
{
	trigger_error("Cannot perform trio analysis without target region (processing systems of child is '".$sys["name_short"]."')!", E_USER_ERROR);
}

//create output folder if missing
$out_folder .= "/";
if (!file_exists($out_folder)) mkdir($out_folder);

//(1.1) check genders of parents
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


//(1.2) determine gender of index
list($stdout) = $parser->exec(get_path("ngs-bits")."SampleGender", "-method hetx -in $c", true);
$gender = trim(substr($stdout[count($stdout)-1], 7));
if (starts_with($gender,"unknown"))
{
	trigger_error("Gender of index could not be determined!", E_USER_ERROR);
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

//(3) variant calling of all three samples together (with very conservative parameters)
$vcf_all = $out_folder.basename($c, ".bam")."_all.vcf.gz";
if ($start=="check" || $start=="vc")
{
	$parser->execTool("NGS/vc_freebayes.php", "-bam $c $m $f -out $vcf_all -target $target_file -min_mq 20 -min_af 0.1 -build ".$sys['build'], true);	
}

//(4) convert VCF to file with only child column, but with TRIO annotation
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
		//write sample headers			
		fwrite($h2, gsvar_sample_header(basename($c, ".bam"), array("ID"=>"genotype", "DiseaseStatus"=>"affected", "SampleName"=>basename($c, ".bam"))));
		fwrite($h2, gsvar_sample_header(basename($f, ".bam"), array("DiseaseStatus"=>"unaffected", "Gender"=>"male")));
		fwrite($h2, gsvar_sample_header(basename($m, ".bam"), array("DiseaseStatus"=>"unaffected", "Gender"=>"female")));
		fwrite($h2, "##ANALYSISTYPE=GERMLINE_TRIO\n");

		//determine indices for each sample	
		$parts = explode("\t", $line);
		$c_idx = vcf_column_index(basename($c, ".bam"), $parts); 
		$f_idx = vcf_column_index(basename($f, ".bam"), $parts); 
		$m_idx = vcf_column_index(basename($m, ".bam"), $parts); 
		fwrite($h2, implode("\t", array_slice($parts, 0, 9))."\t".$parts[$c_idx]."\n");
	}
	else //content
	{
		$parts = explode("\t", $line);
		$format = explode(":", $parts[8]);
		if ($parts[$c_idx]==".") continue;
		$c_info = extract_info($format, $parts[$c_idx]);
		$m_info = extract_info($format, $parts[$m_idx]);
		$f_info = extract_info($format, $parts[$f_idx]);
		
		//update format field and remove mother/father
		$parts[8] = $parts[8].":TRIO";
		fwrite($h2, implode("\t", array_slice($parts, 0, 9))."\t".$parts[$c_idx].":".implode(",", $c_info).",".implode(",", $m_info).",".implode(",", $f_info)."\n");
	}
}
fclose($h1);
fclose($h2);

//(5) zip variant list
$vcf_zipped = $out_folder.basename($c, ".bam")."_var.vcf.gz";
$parser->exec("bgzip", "-c $vcf > $vcf_zipped", false); //no output logging, because Toolbase::extractVersion() does not return
$parser->exec("tabix", "-p vcf $vcf_zipped", false); //no output logging, because Toolbase::extractVersion() does not return

//(6) basic annotation
$parser->execTool("Pipelines/annotate.php", "-out_name ".basename($c, ".bam")." -out_folder $out_folder -system $system");

//(7) add trio annotation column
$gsvar = $out_folder.basename($c, ".bam").".GSvar";
$parser->exec(get_path("ngs-bits")."TrioAnnotation", "-in $gsvar -out $gsvar -gender $gender", true);

?>
