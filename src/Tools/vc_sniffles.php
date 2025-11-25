<?php

/**
  @page vc_sniffles
  
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_sniffles", "Call of structural variants with sniffles. Creates an VCF file.");
$parser->addInfileArray("bam", "Indexed and sorted BAM file(s).", false);
$parser->addOutfile("out", "Output VCF file (gzipped and tabix indexed).", false);
//optional
$parser->addStringArray("sample_ids", "Optional sample id(s)/name(s) for VCF output.", true);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
$parser->addFlag("somatic", "Use somatic mode for SV calling (only single sample).");
$parser->addFlag("use_tandem_repeat_file", "Use tandem repeat BED file for calling (only supported for GRCh38).");
$parser->addFlag("add_unfiltered_output", "Also output unfiltered VCF (only single sample).");
$parser->addFlag("include_mosaic", "Also do mosaic calling for germline samples (since sniffles 2.7).");
$parser->addInfile("target",  "Optional target region to limit SV calls to certain areas", true, true);
$parser->addInt("threads", "Number of threads used.", true, 4);

extract($parser->parse($argv));

//get sample names
if(!is_null($sample_ids) && (count($sample_ids) != 0))
{
	//check provided sample names
	if(count($sample_ids) != count($sample_ids)) trigger_error("Number of provided BAM files and sample names does not match!", E_USER_ERROR);
}
else
{
	$sample_ids = array();
	foreach($bam as $single_sample_bam)
	{
		$sample_ids[] = basename2($single_sample_bam);
	}
}

if(count($bam) == 0)
{
	trigger_error("No BAM file(s) provided!", E_USER_ERROR);
}
else if(count($bam) == 1)
{
  	//single sample
	$single_sample_bam = $bam[0];

	// get sample name
	$name = $sample_ids[0];

	//add tandem repeat file
	if ($use_tandem_repeat_file)
	{
		if ($build != "GRCh38") trigger_error("Tadem repeat file only supported for GRCh38!", E_USER_ERROR);
		$args[] = "--tandem-repeats ".get_path("data_folder")."/dbs/tandem-repeats/human_GRCh38_no_alt_analysis_set.trf.bed";
	}

	$tmp_vcf = $parser->tempFile(".vcf", "sniffles");
    $args = array();
	$in_files = array();
    $args[] = "--input ".$single_sample_bam;
    $args[] = "--vcf ".$tmp_vcf;
    $args[] = "--threads ".$threads;
    $args[] = "--reference ".genome_fasta($build);
    $args[] = "--allow-overwrite";
    $args[] = "--sample-id ".$name;
	$args[] = "--output-rnames";
    if($somatic) $args[] = "--non-germline";
    if($include_mosaic) $args[] = "--mosaic-include-germline";
	if(isset($target))
	{
		$args[] = "--regions ".$target;
		$in_files[] = $target;
	}

	//set bind path for sniffles container
	$in_files[] = $single_sample_bam;
	$in_files[] = genome_fasta($build);

	//execute sniffles container
	$parser->execApptainer("sniffles", "sniffles", implode(" ", $args), $in_files);


	if ($add_unfiltered_output)
	{
		//also create second call file without filtering
		$args = array();
		$in_files = array();
		$tmp_vcf_nofilter = $parser->tempFile(".vcf", "sniffles");
		$args[] = "--input ".$single_sample_bam;
		$args[] = "--vcf ".$tmp_vcf_nofilter;
		$args[] = "--threads ".$threads;
		$args[] = "--reference ".genome_fasta($build);
		$args[] = "--no-qc";
		$args[] = "--sample-id ".$name;
		$args[] = "--output-rnames";
		if($somatic) $args[] = "--non-germline";
		if($include_mosaic) $args[] = "--mosaic-include-germline";
		if(isset($target))
		{
			$args[] = "--regions ".$target;
			$in_files[] = $target;
		} 
		//set bind path for sniffles container
		$in_files[] = $single_sample_bam;
		$in_files[] = genome_fasta($build);
		//execute sniffles container
		$parser->execApptainer("sniffles", "sniffles", implode(" ", $args), $in_files);
	}
}
else
{
	//multisample
	if($somatic) trigger_error("Somatic calling not supported for multisample!", E_USER_ERROR);

	//create single sample snf files
	$snfs = array();
	foreach($bam as $single_sample_bam)
	{
		// get sample name
		$name = $sample_ids[count($snfs)];
		$tmp_snf = $parser->tempFile(".snf", "sniffles");
		$args = array();
		$in_files = array();
    	$args[] = "--input ".$single_sample_bam;
    	$args[] = "--snf ".$tmp_snf;
    	$args[] = "--threads ".$threads;
    	$args[] = "--reference ".genome_fasta($build);
    	$args[] = "--allow-overwrite";
		$args[] = "--sample-id ".$name;
		$args[] = "--output-rnames";
		if($include_mosaic) $args[] = "--mosaic-include-germline";
		if(isset($target))
		{
			$args[] = "--regions ".$target;
			$in_files[] = $target;
		} 

		//set bind path for sniffles container
		$in_files[] = $single_sample_bam;
		$in_files[] = genome_fasta($build);

		//execute sniffles container
		$parser->execApptainer("sniffles", "sniffles", implode(" ", $args), $in_files);

		$snfs[] = $tmp_snf;
	}
	
	//combine snfs to multisample VCF
	$tmp_vcf = $parser->tempFile(".vcf", "sniffles");

	$args = array();
	$in_files = array();
	$args[] = "--input ".implode(" ", $snfs);
	$args[] = "--vcf ".$tmp_vcf;
	$args[] = "--threads ".$threads;
	$args[] = "--reference ".genome_fasta($build);
	$args[] = "--allow-overwrite";
	if(isset($target))
	{
		$args[] = "--regions ".$target;
		$in_files[] = $target;
	} 
	$in_files[] = genome_fasta($build);

	//execute sniffles container
	$parser->execApptainer("sniffles", "sniffles", implode(" ", $args), $in_files);
}

//mask low evidence wildtype calls (single sample only)
if(count($bam) == 1)
{
	$masked_tmp_vcf = $parser->tempFile(".vcf", "sniffles_masked");
	$parser->execApptainer("ngs-bits", "SnifflesVcfFix", "-in {$tmp_vcf} -out {$masked_tmp_vcf}");
	$tmp_vcf = $masked_tmp_vcf;
}

if ($add_unfiltered_output)
{
	//copy unfiltered to output
	$out_nofilter = dirname($out)."/".basename($out, ".vcf.gz")."_no_filter.vcf.gz";
	$parser->execApptainer("htslib", "bgzip", "-c {$tmp_vcf_nofilter} > {$out_nofilter}", [], [dirname($out)]);
	$parser->execApptainer("htslib", "tabix", "-f -p vcf {$out_nofilter}", [], [dirname($out)]);
}

//rarely single SVs are in the wrong order:
$sorted_tmp_vcf = $parser->tempFile(".vcf", "sniffles_sorted");
$parser->execApptainer("ngs-bits", "VcfSort", "-in {$tmp_vcf} -out {$sorted_tmp_vcf}");
$tmp_vcf = $sorted_tmp_vcf;

//add name/pipeline info to VCF header
$vcf = Matrix::fromTSV($tmp_vcf);
$comments = $vcf->getComments();
$comments[] = "#reference=".genome_fasta($build)."\n";
$comments[] = "#PIPELINE=".repository_revision(true)."\n";
if(count($bam) == 1)
{
	$comments[] = gsvar_sample_header($name, array("DiseaseStatus"=>"affected"), "#", "");
}
else
{
	foreach($sample_ids as $id)
	{
		$comments[] = gsvar_sample_header($id, array(), "#", "");
	}
}
$vcf->unique(); //remove duplicate variant calls 
$vcf->setComments($comments);
$vcf->toTSV($tmp_vcf);

$parser->execApptainer("htslib", "bgzip", "-c $tmp_vcf > $out", [], [dirname($out)]);

//index output file
$parser->execApptainer("htslib", "tabix", "-f -p vcf $out", [], [dirname($out)]);

?>