<?php
/** 
	@page contamination_detection
	
	Test data: DX180347_01 is contaminated with DX180330_01
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("contamination_detection", "Detects which sample caused the contamination of a second sample.");
$parser->addString("ps", "Processed sample name of contaminated sample.", false);
$parser->addInfile("roi",  "Target region for variant calling. If unset, determined from 'ps'.", true);
$parser->addInt("min_var_qual", "Minimum quality for variant calls.", true,  30);
$parser->addInt("min_var_depth", "Minimum depth for variant calls.", true, 50);
$parser->addInt("min_var_mapq", "Minimum mapping quality for variant calls.", true, 55);
$parser->addFloat("max_af", "Maximum gnomAD allele frequency of variant used for comparison.", true, 0.02);
$parser->addFloat("max_ngsd", "Maximum occurences in NGSD of variants used for comparison.", true, 50);
$parser->addFloat("min_hits", "Minimum number of variant hits for reporting a sample.", true, 3);
$parser->addInfile("ann_vcf", "Annotated VCF with variants which should be used for comparison (by a previous run with 'gsvar_output' option).", true);
$parser->addOutfile("gsvar_output", "Output GSvar file containing the (filtered) detected variants used for comparison.", true);
$parser->addInt("threads", " The maximum number of threads used for variant calling.", true, 5);

extract($parser->parse($argv));

if (is_null($gsvar_output) || $gsvar_output=="")
{
	$gsvar_output = $parser->tempFile(".GSvar");
}

function vcf_lines($filename)
{
	$lines = 0;
	
	$h = fopen2($filename, 'r');
	while(!feof($h))
	{
		$line = trim(fgets($h));
		if ($line=="" || $line[0]=='#') continue;
		
		++$lines;
	}
	fclose($h);
	
	return $lines;
}

//init
$db = DB::getInstance("NGSD");
$ngsbits = get_path("ngs-bits");
$info = get_processed_sample_info($db, $ps);

//roi
if ($roi=="")
{
	$dummy = "";
	$sys = load_system($dummy, $ps);
	$roi = trim($sys['target_file']);
}
if ($roi=="")
{
	trigger_error("No ROI set for variant calling!", E_USER_ERROR);
}	

//get required files
$ps_bam = $info['ps_bam'];
$ps_vcf = $info['ps_folder']."/".$ps."_var.vcf.gz";

//perform additional VC and filtering if no precalculated file is provided
if (is_null($ann_vcf) || $ann_vcf=="")
{
	//additional VC
	$vc_vcf = $parser->tempFile("_vc_noFilter.vcf.gz");
	$args = array();
	$args[] = "-bam $ps_bam";
	$args[] = "-out $vc_vcf";
	$args[] = "-threads $threads";
	$args[] = "-min_af 0.01";
	$args[] = "-min_mq 50";
	$args[] = "-min_bq 25";
	$args[] = "-no_ploidy";
	$args[] = "-target ".$roi;
	$args[] = "-target_extend 50";
	$parser->execTool("NGS/vc_freebayes.php", implode(" ", $args));

	//parse original VCF and create a array with all variants
	$decomp_ps_vcf = $parser->tempFile("_decomp.vcf");
	$parser->exec("zcat", "$ps_vcf > $decomp_ps_vcf");
	$vars = array();
	$fh = fopen2($decomp_ps_vcf, "r");
	$n_var = 0;
	while(!feof($fh))
	{
		$line = trim(fgets($fh));
		// skip comments and header
		if ($line=="" || $line[0]=="#") continue;

		// extract variant info
		list($chr, $pos, , $ref, $obs) = explode("\t", $line);

		//store info:
		$vars["$chr:$pos $ref>$obs"] = true;
		$n_var++;
	}
	fclose($fh);
	print "##Germline VCF file contains {$n_var} variants.\n";

	// parse new vcf and remove all variants which are already found in first run
	$decomp_vc_vcf = $parser->tempFile("_vc_decomp.vcf");
	$parser->exec("zcat", "$vc_vcf > $decomp_vc_vcf");
	$input_fh = fopen2($decomp_vc_vcf, "r");
	$unique_vcf = $parser->tempFile("_unique.vcf");
	$output_fh = fopen2($unique_vcf, "w");
	$n_var = 0;
	$n_new_var = 0;
	while(!feof($input_fh))
	{
		$line = trim(fgets($input_fh));
		if ($line=="") continue;
		
		// write comments and header
		if ($line[0]=="#") 
		{
			if (starts_with($line, "##"))
			{
				fwrite($output_fh, $line."\n");
			}
			else
			{
				fwrite($output_fh, "##SAMPLE=<ID={$ps},DiseaseStatus=affected>\n");
				fwrite($output_fh, $line."\n");
			}
			continue;
		}
		// extract variant info
		list($chr, $pos, , $ref, $obs) = explode("\t", $line);
		$n_var++;

		// skip all variants which are already found in original VC
		if (isset($vars["$chr:$pos $ref>$obs"])) continue;

		// write unique variants
		fwrite($output_fh, $line."\n");
		$n_new_var++;
	}
	fclose($input_fh);
	fclose($output_fh);
	print "##Low-AF variant calling found $n_var variants.\n";

	//store annotate file in output directory for additional parameter tuning
	$ann_vcf = strtr($gsvar_output, [".GSvar"=>"_ann.vcf"]);
	
	//annotate VCF
	$args = [];
	$args[] = "-in ".$unique_vcf;
	$args[] = "-out ".$ann_vcf;
	$args[] = "-threads ".$threads;
	$args[] = "-ps_name ".$ps;
	$parser->execTool("NGS/an_vep.php", implode(" ", $args));
}
else
{
	$n_new_var = vcf_lines($ann_vcf);
}
print "##{$n_new_var} variant found that are not in the germline VCF.\n";

//filter VCF for variant quality
$filtered_vcf = $parser->tempFile("_filtered.vcf");
$parser->exec("{$ngsbits}VcfFilter", "-qual $min_var_qual -sample 'DP >= $min_var_depth' -info 'MQM >= $min_var_mapq' -in $ann_vcf -out $filtered_vcf");

//convert to GSvar
$gsvar = $parser->tempFile(".GSvar");
$parser->execTool("NGS/vcf2gsvar.php", "-in $filtered_vcf -out $gsvar");
print "##".vcf_lines($gsvar)." variant remain after quality filter (QUAL={$min_var_qual}, DP={$min_var_depth}, MQM={$min_var_mapq})\n";

//filter GSvar file for AF and NGSD count
$filters = [
	"Allele frequency	max_af=".(100.0*$max_af),
	"Allele frequency (sub-populations)	max_af=".(100.0*$max_af),
	"Count NGSD	max_count={$max_ngsd}	ignore_genotype=false",
	"Genotype affected	genotypes=het"
	];
$filter_file = $parser->tempFile("_filters.txt");
file_put_contents($filter_file, implode("\n", $filters));

$parser->exec("{$ngsbits}VariantFilterAnnotations", "-filters {$filter_file} -in {$gsvar} -out {$gsvar_output}");
print "##".vcf_lines($gsvar_output)." variant remain after frequency filter (AF={$max_af}, NGSD={$max_ngsd})\n";

//find variants in NGSD
$vars_refound = [];
$c_vars_ngsd = 0;
$gsvar_filtered = file($gsvar_output);
foreach($gsvar_filtered as $line)
{
	$line = trim($line);
	if ($line=="" || $line[0]=="#") continue;
	
	//extract variant info
	list($chr, $start, $end, $ref, $obs) = explode("\t", $line);
	
	$res = $db->executeQuery("SELECT id, gnomad FROM variant WHERE chr='{$chr}' AND start='{$start}' AND end='{$end}' AND ref='{$ref}' AND obs='{$obs}'");
	if (count($res)==0) continue;
	++$c_vars_ngsd;
	
	$var_id = $res[0]['id'];
	// $af_1000g = $res[0]['1000g'];
	$af_gnomad = $res[0]['gnomad'];
	
	//determine samples that contain the variant
	$res2 = $db->executeQuery("SELECT processed_sample_id, genotype FROM detected_variant WHERE variant_id={$var_id}");
	foreach($res2 as $row)
	{
		$ps_id = $row['processed_sample_id'];
		$gt = $row['genotype'];
		
		$vars_refound[$ps_id][] = "{$chr}:{$start} {$ref}>{$obs} GT={$gt} gnomAD_af={$af_gnomad} NGSD_count=".count($res2);
	}
}
print "##{$c_vars_ngsd} variant found in NGSD.\n";

//print output
print "#sample\toverlap_perc\toverlap_count\toverlap_variants\n";
foreach($vars_refound as $ps_id => $vars)
{
	$ol_count = count($vars);
	if ($ol_count<$min_hits) continue;
	
	$overlap_perc = number_format(100.0*$ol_count/$c_vars_ngsd, 2);
	$ps2 = $db->getValue("SELECT CONCAT(s.name,'_',LPAD(ps.process_id,2,'0')) FROM processed_sample as ps, sample as s WHERE ps.sample_id = s.id AND ps.id={$ps_id}");
	
	//skip same samples
	if (explode("_", $ps)[0]==explode("_", $ps2)[0]) continue;
	
	print "{$ps2}\t{$overlap_perc}\t{$ol_count}\t".implode(" // ", $vars)."\n";
}

?>
