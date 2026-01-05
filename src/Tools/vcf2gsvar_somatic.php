<?php 
/** 
	@page vcf2gsvar_somatic
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("vcf2gsvar_somatic", "Converts an annotated VCF file to a GSvar file.");
$parser->addInfile("in", "Input file in VCF format.", false);
$parser->addOutfile("out", "Output file in GSvar format.", false);
$parser->addString("t_col", "Column name of tumor sample.", false);
//optional
$parser->addString("n_col", "Column name of normal sample.", true, "na");
$parser->addFlag("updown", "Don't discard up- or downstream annotations (5000 bases around genes).");
$parser->addFlag("cfdna", "Set the analysis type to cfDNA");
$parser->addEnum("db", "Database to connect to", true, db_names(), "NGSD");
extract($parser->parse($argv));

//init
$tumor_only = empty($n_col) || $n_col == "na";

//convert base data with vcf2gsvar
$tmp = $parser->tempFile(".GSvar");
$args = ["-in $in", "-out $tmp", "-genotype_mode skip", "-custom custom_columns"];
if ($updown) $args[] = "-updown";
$parser->execTool("Tools/vcf2gsvar.php", implode(" ", $args));
$gsvar = Matrix::fromTSV($tmp);

//add meta data to header
$comments = [];
if ($cfdna)
{
	$comments[] = "#ANALYSISTYPE=CFDNA";
}
else
{
	$comments[] = "#ANALYSISTYPE=SOMATIC_".($tumor_only ? "SINGLESAMPLE" : "PAIR");
}
$comments[] = "#PIPELINE=".repository_revision(true);
$gsvar->setComments(array_merge($comments, $gsvar->getComments()));

//add tumor/normal AF/DP columns (before filter column)
$index = $gsvar->getColumnIndex("filter");
$empty_col = array_fill(0, $gsvar->rows(), "");
$gsvar->insertCol($index, $empty_col, "tumor_af", "Mutant allele frequency in tumor (Sample ".$t_col.").");
$i_tumor_af = $index++;
$gsvar->insertCol($index, $empty_col, "tumor_dp", "Tumor Depth (Sample ".$t_col.").");
$i_tumor_dp = $index++;
if(!$tumor_only)
{
	$gsvar->insertCol($index, $empty_col, "normal_af", "Mutant allele frequency in normal (Sample ".$n_col.").");
	$i_normal_af = $index++;
	$gsvar->insertCol($index, $empty_col, "normal_dp", "Normal depth (Sample ".$n_col.").");
	$i_normal_dp = $index++;
}

//process VCF header (extract variant caller and tumor/normal index)
$handle = fopen2($in, "r");
$var_caller = false;
$source_lines = 0;
$tumor_idx = false;
$normal_idx = false;
while(!feof($handle))
{
	$line = nl_trim(fgets($handle));
	if($line=="") continue;
	
	//variant caller
	if(starts_with($line, "##source="))	
	{
		++$source_lines;
		$line = strtolower($line);
		if(contains($line, "freebayes")) $var_caller = "freebayes";
		if(contains($line, "strelka")) $var_caller = "strelka";
		if(contains($line, "varscan2") ) $var_caller = "varscan2";
		if(contains($line, "umivar2") ) $var_caller = "umivar2";
		if(contains($line, "dragen") ) $var_caller = "dragen";
		if(contains($line, "deepsomatic")) $var_caller = "deepsomatic";
		echo $line."\n";

		//add umivar specific columns
		if ($var_caller == "umivar2")
		{
			$gsvar->insertCol($index, $empty_col, "p-value", "Uncorrected p-value");
			$i_p_value = $index++;
			$gsvar->insertCol($index, $empty_col, "Alt counts", "");
			$i_alt_counts = $index++;
			
			//add multi UMI columns
			$gsvar->insertCol($index, $empty_col, "m_AF", "Multi-UMI alternative allele fraction");
			$i_m_af = $index++;
			$gsvar->insertCol($index, $empty_col, "m_REF", "Multi-UMI reference-like read counts");
			$i_m_ref = $index++;
			$gsvar->insertCol($index, $empty_col, "m_ALT", "Multi-UMI alternative read counts");
			$i_m_alt = $index++;

			$gsvar->insertCol($index, $empty_col, "Strand", "Alleles in strands: Alternative forward, Alternative reverse, Reference forward, Reference reverse");
			$i_strand = $index++;
			$gsvar->insertCol($index, $empty_col, "Sequence context", "Sequence +/- 5bp around the variant.");
			$i_seq_context = $index++;
			$gsvar->insertCol($index, $empty_col, "Homopolymer", "Indicates if variant is in a homopolymer region.");
			$i_homopolymer = $index++;
		}
	}
	
	//header line
	if(starts_with($line, "#CHROM")) 
	{
		$parts = explode("\t", $line);
		$tumor_idx = array_search($t_col, $parts);
		$normal_idx = array_search($n_col, $parts);
		break;
	}
}
if($source_lines!=1) trigger_error("Found ".($source_lines==0 ? "no" : "several")." 'source' entries in VCF header (needed to identify the variant caller).", E_USER_ERROR);
if($var_caller===false) trigger_error("Unknown variant caller from 'source' entry in VCF header.", E_USER_ERROR);
if($tumor_idx===false) trigger_error("Could not identify tumor column '$t_col' in VCF header.", E_USER_ERROR);
if(!$tumor_only && $normal_idx===false && $var_caller !== "deepsomatic") trigger_error("Could not identify normal column '$n_col' in VCF header.", E_USER_ERROR);

//process variants
$r = -1;
$i_quality = $gsvar->getColumnIndex("quality");
while(!feof($handle))
{
	$line = nl_trim(fgets($handle));
	if ($line=="" || $line[0]=="#") continue;	
	
	$cols = explode("\t", $line);
	if (count($cols)<9) trigger_error("VCF file line contains less than 10 columns: $line\n", E_USER_ERROR);
	list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format) = $cols;
	
	//skip bad chromosomes
	if(chr_check($chr, 22, false) === FALSE) continue; 
	++$r;
	
	//abort if multi-allelic
	if (contains($alt, ",")) trigger_error("Multi-allelic variants not supported: $chr:$pos $ref>$alt", E_USER_ERROR);
		
	//variant type
	$is_indel = strlen($ref) > 1 || strlen($alt) > 1;

	//check that we annotate the right variant
	$start = $pos;
	$end = $start;
	if($is_indel) //correct indels
	{
		list($start, $end, $ref, $alt) = correct_indel($start, $ref, $alt);
	}
	if (array($chr, $start, $end, $ref, $alt) != array_slice($gsvar->getRow($r), 0, 5))
	{
		trigger_error("Internal error - this should not happen!\nFound variant: ".implode(" ", array($chr, $start, $end, $ref, $alt))."\tExpected variant: ".implode(" ", array_slice($gsvar->getRow($r), 0, 5)), E_USER_ERROR);
	}

	//set quality
	$snp_q = $qual;
	if ($var_caller=="strelka")
	{
		$field_name = $is_indel ? "QSI=" : "QSS=";
		foreach(explode(";", $info) as $entry)
		{
			if (starts_with($entry, $field_name))
			{
				$snp_q = substr($entry, strlen($field_name));
			}
		}
	}
	else if ($var_caller=="umivar2")
	{
		$f_parts = explode(":", $format);
		$i_p_value = array_search("Pval", $f_parts);
		if($i_p_value!==false)
		{
			$parts = explode(":", $cols[$tumor_idx]);
			$snp_q = -10 * log10($parts[$i_p_value]); //calculate phred score from pvalue
			$snp_q = min(floor($snp_q), 255);
		}
	}
	else if ($var_caller=="dragen")
	{
		$f_parts = explode(":", $format);
		$i_sq_value = array_search("SQ", $f_parts);
		
		if($i_sq_value!==false)
		{
			$parts = explode(":", $cols[$tumor_idx]);
			$snp_q = $parts[$i_sq_value];
		}
	}
	$gsvar->set($r, $i_quality, "QUAL={$snp_q}");

	//calculate DP/AF tumor
	if($var_caller=="freebayes")
	{
		list($tumor_dp, $tumor_af) = vcf_freebayes($format, $cols[$tumor_idx]);
	}
	else if($var_caller=="strelka")
	{
		if($is_indel)
		{
			list($tumor_dp, $tumor_af) = vcf_strelka_indel($format, $cols[$tumor_idx]);
		}
		else
		{
			list($tumor_dp, $tumor_af) = vcf_strelka_snv($format, $cols[$tumor_idx], $alt);
		}
	}
	else if ($var_caller == "dragen")
	{
		list($tumor_dp, $tumor_af) = vcf_dragen_var($format, $cols[$tumor_idx]);
	}
	else if($var_caller == "varscan2")
	{
		list($tumor_dp, $tumor_af) = vcf_varscan2($format, $cols[$tumor_idx]);
	}
	else if($var_caller == "umivar2")
	{
		list($tumor_dp, $tumor_af, $p_value, $alt_count, $strand, $seq_context, $homopolymer, $m_af, $m_ref, $m_alt) = vcf_umivar2($filter, $info, $format, $cols[$tumor_idx]);

		// add umiVar2 specific values
		if($var_caller == "umivar2")
		{
			$gsvar->set($r, $i_p_value, $p_value);
			$gsvar->set($r, $i_alt_counts, $alt_count);
			
			$gsvar->set($r, $i_m_af, $m_af);
			$gsvar->set($r, $i_m_ref, $m_ref);
			$gsvar->set($r, $i_m_alt, $m_alt);

			$gsvar->set($r, $i_strand, $strand);
			$gsvar->set($r, $i_seq_context, $seq_context);
			$gsvar->set($r, $i_homopolymer, $homopolymer);
		}
	}
	else if ($var_caller == "deepsomatic")
	{
		list($tumor_dp, $tumor_af) = vcf_deepvariant($format, $cols[$tumor_idx]);
	}
	
	$gsvar->set($r, $i_tumor_dp, $tumor_dp);
	$gsvar->set($r, $i_tumor_af, number_format($tumor_af, 4));

	//calculate DP/AF normal
	if (!$tumor_only)
	{
		if($var_caller=="freebayes")
		{
			list($normal_dp, $normal_af) = vcf_freebayes($format, $cols[$normal_idx]);
		}
		else if($var_caller=="strelka")
		{
			if($is_indel)
			{
				list($normal_dp, $normal_af) = vcf_strelka_indel($format, $cols[$normal_idx]);
			}
			else
			{
				list($normal_dp, $normal_af) = vcf_strelka_snv($format, $cols[$normal_idx], $alt);
			}
		}
		else if ($var_caller == "dragen")
		{
			list($normal_dp, $normal_af) = vcf_dragen_var($format, $cols[$normal_idx]);
		}
		else if ($var_caller == "deepsomatic")
		{
			list($normal_dp, $normal_af) = vcf_deepsomatic($info, $n_col);
		}
		$gsvar->set($r, $i_normal_dp, $normal_dp);
		$gsvar->set($r, $i_normal_af, number_format($normal_af, 4));
	}
}
gzclose($handle);

//store output
$gsvar->toTSV($out);

?>
