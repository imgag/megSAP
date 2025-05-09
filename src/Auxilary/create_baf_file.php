<?php
/** 
	@page create_baf_file
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("create_baf_file", "Creates BAF file from VCF file.");
$parser->addOutfile("out_file", "Target to store BAF file.", false);
$parser->addInfile("vcf", "Input VCF file.", false);
$parser->addString("s_col", "Column name of sample.", false);
//optional
$parser->addString("build", "The genome build to use.", true, "GRCh38");
extract($parser->parse($argv));

//Abort if out_file exists to prevent interference with other jobs
if(file_exists($out_file)) return;
$genome = genome_fasta($build);

$in_handle  = gzopen2($vcf,"r");
$out_handle = fopen2($out_file,"w");

//process VCF header (extract variant caller and sample index)
$var_caller = false;
$source_lines = 0;
$sample_idx = false;
while(!feof($in_handle))
{
	$line = nl_trim(fgets($in_handle));
	if($line=="") continue;
	
	//variant caller
	if(starts_with($line, "##source="))	
	{
		#++$source_lines;
		$line = strtolower($line);
		if(contains($line, "freebayes")) $var_caller = "freebayes";
		if(contains($line, "strelka")) $var_caller = "strelka";
		if(contains($line, "umivar2") ) $var_caller = "umivar2";
		if(contains($line, "dragen") ) $var_caller = "dragen";
		if(contains($line, "deepvariant") ) $var_caller = "deepvariant";

	}
	
	//header line
	if(starts_with($line, "#CHROM")) 
	{
		$parts = explode("\t", $line);
		$sample_idx = array_search($s_col, $parts);
		break;
	}
}
#if($source_lines!=1) trigger_error("Found ".($source_lines==0 ? "no" : "several")." 'source' entries in VCF header (needed to identify the variant caller).", E_USER_ERROR);
if($var_caller===false) trigger_error("Unknown variant caller from 'source' entry in VCF header.", E_USER_ERROR);
if($sample_idx===false) trigger_error("Could not identify sample column '$s_col' in VCF header.", E_USER_ERROR);

//process variants
while(!feof($in_handle))
{
    $line = trim(fgets($in_handle));
    if(starts_with($line,"#")) continue;
    if(empty($line)) continue;
    
	$cols = explode("\t", $line);
	if (count($cols)<9) trigger_error("VCF file line contains less than 10 columns: $line\n", E_USER_ERROR);
	list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format) = $cols;

    if(strlen($ref) != 1 || strlen($alt) != 1 ) continue; //Skip INDELs
	if($chr=="chrMT" || $filter == "mosaic" || $filter == "low_mappability") continue; //Skip mito, mosaic and low_mappability vars

    //calculate DP/AF
	if($var_caller=="freebayes")
	{
		list($dp, $af) = vcf_freebayes($format, $cols[$sample_idx]);
	}
	else if($var_caller=="strelka")
	{
    	list($dp, $af) = vcf_strelka_snv($format, $cols[$sample_idx], $alt);
	}
	else if ($var_caller == "dragen")
	{
		list($dp, $af) = vcf_dragen_var($format, $cols[$sample_idx]);
	}
	else if($var_caller == "umivar2")
	{
		list($dp, $af, $p_value, $alt_count, $strand, $seq_context, $homopolymer, $m_af, $m_ref, $m_alt) = vcf_umivar2($filter, $info, $format, $cols[$sample_idx]);
	}
	else if ($var_caller =="deepvariant")
	{
		list($dp, $af) = vcf_deepvariant($format, $cols[$sample_idx]);
	}

    fputs($out_handle,"{$chr}\t{$pos}\t{$pos}\t{$chr}_{$pos}\t{$af}\t{$dp}\n");
}
gzclose($in_handle);
fclose($out_handle);

?>