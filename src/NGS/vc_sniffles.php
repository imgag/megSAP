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

	//TODO: add tandem repeat file

    $tmp_vcf = $parser->tempFile(".vcf", "sniffles");

    $args = array();
    $args[] = "--input ".$single_sample_bam;
    $args[] = "--vcf ".$tmp_vcf;
    $args[] = "--threads ".$threads;
    $args[] = "--reference ".genome_fasta($build);
    $args[] = "--allow-overwrite";
    $args[] = "--sample-id ".$name;
    if($somatic) $args[] = "--non-germline";

    $parser->exec(get_path("sniffles"), implode(" ", $args), true);

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
    	$args[] = "--input ".$single_sample_bam;
    	$args[] = "--snf ".$tmp_snf;
    	$args[] = "--threads ".$threads;
    	$args[] = "--reference ".genome_fasta($build);
    	$args[] = "--allow-overwrite";
		$args[] = "--sample-id ".$name;

		$parser->exec(get_path("sniffles"), implode(" ", $args), true);

		$snfs[] = $tmp_snf;
	}
	

	//combine snfs to multisample VCF
	$tmp_vcf = $parser->tempFile(".vcf", "sniffles");

	$args = array();
	$args[] = "--input ".implode(" ", $snfs);
	$args[] = "--vcf ".$tmp_vcf;
	$args[] = "--threads ".$threads;
	$args[] = "--reference ".genome_fasta($build);
	$args[] = "--allow-overwrite";

    $parser->exec(get_path("sniffles"), implode(" ", $args), true);

}

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

$vcf->setComments(sort_vcf_comments($comments));
$vcf->toTSV($tmp_vcf);

$parser->exec("bgzip", "-c $tmp_vcf > $out", false);

//index output file
$parser->exec("tabix", "-f -p vcf $out", false); //no output logging, because Toolbase::extractVersion() does not return

?>