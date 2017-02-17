<?php
/**
	@page filter_vcf
	@todo revisit somatic and somatic_ds (freebayes!) filter
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("filter_vcf", "Filter VCF-files according to different filter criteria.");
$parser->addInfile("in", "Input variant file in VCF format containing all necessary columns (s. below for each filter).", false);
$parser->addOutfile("out", "Output variant file in VCF format.", false);
$filter = array('somatic', 'somatic_ds', 'coding', 'non_synonymous', 'somatic_diag_capa', 'iVac','all');
$parser->addString("type", "Filter set to use, can be comma delimited. Valid are: ".implode(",",$filter).".",false);
//optional
$parser->addFlag("keep", "Keep all variants. Otherwise only variants passing all filters will be kept");
$parser->addString("roi", "Target region BED file (for off-target filter).", true, "");
$parser->addFloat("contamination", "Estimated fraction of tumor cells in normal sample.",true,0.00);
$parser->addFloat("min_af", "Minimum variant allele frequency in tumor.",true,0.05);
extract($parser->parse($argv));

//zipped
$zipped = false;
if(ends_with($in, ".gz"))	$zipped = true;	

//check selected filter types
$types = explode(",",$type);
foreach($types as $t)	if(!in_array($t,$filter))	trigger_error("Unknown filter '".$t."' given.",E_USER_ERROR);
$in_file = Matrix::fromTSV($in);

//reset filter column
$filter_column = "FILTER";
$f = $in_file->getColumnIndex($filter_column, false, false);
if($f!=6)	trigger_error("Wrong column index for filter column! Is ".(emtpy($f)?"empty":$f)." should be 6.",E_USER_ERROR);
$prefix = "fv.";
$comments = $in_file->getComments();
$var_caller = NULL;
if($f!==FALSE)
{
	//clean up filter descriptions
	$comments = array();
	foreach($in_file->getComments() as $c)
	{
		if(strpos($c,"#FILTER=<ID=".$prefix)===0)	continue;	//skip filter_vcf description
		$comments[] = $c;
	}
	$in_file->setComments($comments);
	
	//clean up filter column
	for($i=0;$i<$in_file->rows();++$i)
	{
		$filter = $in_file->get($i,$f);
		$tmp1 = explode(";",$filter);
		$tmp2 = array();
		for($j=0;$j<count($tmp1);++$j)
		{
			if(strpos($tmp1[$j],$prefix)===0)	continue;
			$tmp2[] = $tmp1[$j];
		}
		$in_file->set($i,$f,implode(";",$tmp2));
	}
}

//extract tumor and normal column from pedigree tag
//e.g. ##PEDIGREE=<Tumor_DNA=GS150893,Normal_DNA=GS150892>
$tumor_id = NULL;
$tumor_col = NULL;
$normal_id = NULL;
$normal_col = NULL;
$found = false;
foreach($in_file->getComments() as $comment)
{
	if(strpos($comment,"#PEDIGREE=<")===0)
	{
		if($found==true)	trigger_error("Multiple PEDIGREE comments given!",E_USER_ERROR);
		
		list(,$value) = explode("=",$comment,2);
		$samples = explode(",",trim($value,'<>'));
		foreach($samples as $s)
		{
			list($type,$id) = explode("=",$s);
			if(strpos($type,"Tumor")===0)	
			{
				$tumor_id = $id;
				$tumor_col = $in_file->getColumnIndex($id);
			}
			if(strpos($type,"Normal")===0)
			{
				$normal_id = $id;
				$normal_col = $in_file->getColumnIndex($id);
			}
		}
	}

	if(strpos($comment,"#source")===0)	//information about variant caller
	{
		$tmp_vcaller = "";
		if(stripos($comment,"freebayes")!==FALSE)	$tmp_vcaller = "freebayes";
		if(stripos($comment,"strelka")!==FALSE)	$tmp_vcaller = "strelka";
		if(!is_null($var_caller) && $var_caller!=$tmp_vcaller)	trigger_error("Two different source informations identified ($var_caller != $tmp_vcaller!)",E_USER_ERROR);
		if(!is_null($var_caller) && $var_caller==$tmp_vcaller)	trigger_error("Source information used twice!",E_USER_WARNING);
		$var_caller = $tmp_vcaller;
	}
}
if(is_null($tumor_col))	trigger_error("No tumor column given!",E_USER_ERROR);
if(is_null($var_caller))	trigger_error("Variant caller not identified.",E_USER_ERROR);

//extract relevant MISO terms for filtering
$obo = MISO::fromOBO();
$miso_terms_coding = array();
$ids = array("SO:0001580","SO:0001568");
foreach($ids as $id)
{
	$miso_terms_coding[] = $obo->getTermNameByID($id);
	$tmp = $obo->getTermChildrenNames($id,true);
	$miso_terms_coding = array_merge($miso_terms_coding,$tmp);
}
$miso_terms_coding = array_unique($miso_terms_coding);
//prepare variant_type filter (valid effects for coding & splicing)
//children of synonymous_variant (SO:0001819)
$obo = MISO::fromOBO();
$miso_terms_synonymous = array();
$ids = array("SO:0001819");
foreach($ids as $id)
{
	$miso_terms_synonymous[] = $obo->getTermNameByID($id);
	$tmp = $obo->getTermChildrenNames($id,true);
	$miso_terms_synonymous = array_merge($miso_terms_synonymous,$tmp);
}
$miso_terms_synonymous = array_unique($miso_terms_synonymous);

//compare to bed-file
$targets = array();
if($roi!="") //skip for enrichments without target
{
	$target_bed = Matrix::fromTSV($roi);
	for($i=0;$i<$target_bed->rows();++$i)
	{
		$row = $target_bed->getRow($i);
		if(!isset($targets[$row[0]])) $targets[$row[0]] = array();
		$targets[$row[0]][] = array($row[1]+1, $row[2]);
	}
}

//set FILTER column variants
//example for filter description (vcf format 4.1):
//FILTER - filter status: PASS if this position has passed all filters, i.e. a call is made at this position. Otherwise,
//if the site has not passed all filters, a semicolon-separated list of codes for filters that fail. e.g. "q10;s50" might
//indicate that at this site the quality is below 10 and the number of samples with data is below 50% of the total
//number of samples. '0' is reserved and should not be used as a filter String. If filters have not been applied,
//then this field should be set to the missing value (= '.'). (String, no white-space or semi-colons permitted)
////##FILTER=<ID=ID,Description="description">
$filter = array();
$all_filter = array();
for($i=0;$i<$in_file->rows();++$i)
{
	$filter = array();
	$row = $in_file->getRow($i);

	//different QC-values for SNVs and Indels
	$type = "SNV";
	if(strlen($row[3]) > 1)	$type = "INDEL";
	$as = explode(",",$row[4]);
	foreach($as as $a)	//check for different alleles
	{
		if($type=="SNV" && strlen($a) > 1)	$type = "INDEL";
	}
	
	//extract info field
	$genotype = $row[8];
	$info = explode(";", $row[7]);
	$tmp = array();
	foreach($info as $entry)
	{
		if(strpos($entry,"=")!==FALSE)	list($key, $value) = explode("=", $entry, 2);
		else
		{
			$key = $entry;
			$value = NULL;
		}
		$tmp[$key] = $value;
	}
	$info = $tmp;
	
	//use different filter types
	$tmp_filters = array();
	$tmp_col_tum = $row[$tumor_col];
	$tmp_col_nor = null;
	if(!is_null($normal_col))	$tmp_col_nor = $row[$normal_col];
	if(in_array("somatic", $types))
	{
		filter_not_coding_splicing($filter, $info, $miso_terms_coding);
		filter_somatic($filter, $genotype,$tmp_col_tum,$tmp_col_nor,$type,$row[4],$contamination, $min_af, $var_caller);
	}
	if(in_array("somatic_ds", $types))
	{
		filter_not_coding_splicing($filter, $info, $miso_terms_coding);
		filter_somatic_ds($filter, $genotype,$tmp_col_tum,$tmp_col_nor,$type,$var_caller);
	}
	if(in_array("coding", $types))
	{
		filter_not_coding_splicing($filter, $info, $miso_terms_coding);
	}
	if(in_array("non_synonymous", $types))
	{
		filter_synonymous($filter, $info, $miso_terms_coding, $miso_terms_synonymous);
	}
	if(in_array("somatic_diag_capa", $types))
	{
		filter_not_coding_splicing($filter, $info, $miso_terms_coding);
		filter_synonymous($filter, $info, $miso_terms_coding, $miso_terms_synonymous);
		filter_somatic($filter, $genotype,$tmp_col_tum,$tmp_col_nor,$type,$row[4],$contamination, $min_af, $var_caller);
		filter_somatic_capa($filter, $info, $genotype, $tmp_col_tum, $tmp_col_nor, $type,$row[4],$var_caller);
		filter_off_target($filter, $row[0], $row[1], $tumor_id, $normal_id, $targets);
	}
	if(in_array("iVac", $types))
	{
		filter_not_coding_splicing($filter, $info, $miso_terms_coding);
		filter_synonymous($filter, $info, $miso_terms_coding, $miso_terms_synonymous);
		filter_somatic($filter, $genotype,$tmp_col_tum,$tmp_col_nor, $type,$row[4],$contamination, $min_af, $var_caller);
		filter_off_target($filter, $row[0], $row[1], $tumor_id, $normal_id, $targets);	//too slow - increase speed
	}
	if(in_array("all", $types))
	{
		filter_not_coding_splicing($filter, $info, $miso_terms_coding);
		filter_synonymous($filter, $info, $miso_terms_coding, $miso_terms_synonymous);
		filter_somatic($filter, $genotype,$tmp_col_tum,$tmp_col_nor,$type,$row[4],$contamination, $min_af, $var_caller);
		filter_somatic_capa($filter, $info, $genotype, $tmp_col_tum, $tmp_col_nor, $type,$row[4],$var_caller);
		filter_somatic_ds($filter, $genotype,$tmp_col_tum,$tmp_col_nor,$type,$row[4],$var_caller);
		filter_off_target($filter, $row[0], $row[1], $tumor_id, $normal_id, $targets);	//too slow - increase speed
	}
	
	//set PASS criterion if all filters were passed
	$tmp_filter = array();
	foreach($filter as $k => $f)
	{
		if(!isset($f["active"]) || !$f["active"])	continue;
		$tmp_filter[] = $prefix.$k;
	}
	$tmp_filter_old = explode(";",trim($row[6],"."));
	if($tmp_filter_old[0] == "PASS")	unset($tmp_filter_old[0]);
	$tmp_filter = array_filter(array_merge($tmp_filter,$tmp_filter_old));	//keep original filter
	if(empty($tmp_filter))	$tmp_filter[] = "PASS";
	$in_file->set($i,6,implode(";",	$tmp_filter));
}

//add comments and save files
//##FILTER=<ID=ID,Description="description">
$comments = $in_file->getComments();
foreach($filter as $id => $details)	//every line has the same filter criteria, so filter criteria of the las row is sufficient
{
	$comments[] = "#FILTER=<ID=".$prefix.$id.",Description=\"".$details["desc"]."\">";
}
$in_file->setComments($comments);

// freebayes: remove variants that are likely germline variants
if($var_caller=="freebayes")
{
	$filtered = new Matrix();
	$filtered->setComments($in_file->getComments());
	$filtered->setHeaders($in_file->getHeaders());
	$f = $in_file->getColumnIndex("FILTER", false, false);
	for($i=0;$i<$in_file->rows();++$i)
	{
		$row = $in_file->getRow($i);
		if(strpos($row[$f],$prefix."som_af_ratio")!==FALSE)	continue;
		if(strpos($row[$f],$prefix."som_no_depth")!==FALSE)	continue;
		$filtered->addRow($row);
	}
	$in_file = $filtered;	
}

if(!$keep)
{
	$filtered = new Matrix();
	$filtered->setComments($in_file->getComments());
	$filtered->setHeaders($in_file->getHeaders());
	for($i=0;$i<$in_file->rows();++$i)
	{
		$row = $in_file->getRow($i);
		if($row[6]!="PASS")	continue;
		$filtered->addRow($row);
	}
	$in_file = $filtered;
}


if(!$zipped)
{
	$in_file->toTSV($out);
}
else
{
	$tmp_out = $parser->tempFile("_filter.vcf.gz");
	$in_file->toTSV($tmp_out);

	//zip annotated VCF file
	if(strpos($out, ".gz")===FALSE)	$out = $out.".gz";
	$parser->exec("bgzip", "-c $tmp_out > $out", false); //no output logging, because Toolbase::extractVersion() does not return
	$parser->exec("tabix", "-f -p vcf $out", false); //no output logging, because Toolbase::extractVersion() does not return	
}


//filter functions
function filter_synonymous(&$filter, $info, $miso_terms_coding, $miso_terms_synonymous)
{	
	//get variant annotation
	$coding = false;
	$synonymous = true;
	if(isset($info["ANN"]))
	{
		$ann = explode(",", $info["ANN"]);
		foreach($ann as $a)
		{
			$parts = explode("|", $a);
			$variant_types[] = $parts[1];
		}

		//filter out bad rows
		foreach($variant_types as $variant_type)
		{
			if(in_array($variant_type,$miso_terms_coding))
			{
				$coding = true;
				if(!in_array($variant_type,$miso_terms_synonymous))
				{
					$synonymous = false;
				}
			}
		}
	}
	
	add_filter($filter, "syn_var", "Synonymous variant (filter_vcf).");
	if($coding && $synonymous) activate_filter($filter, "syn_var");
}


//filter functions
function filter_not_coding_splicing(&$filter, $info, $miso_terms_coding)
{
	$skip_variant = true;
	if(isset($info["ANN"]))
	{
		//get variant annotation
		$ann = explode(",", $info["ANN"]);
		foreach($ann as $a)
		{
			$parts = explode("|", $a);
			$variant_types[] = $parts[1];
		}

		//filter non coding variants
		foreach($variant_types as $variant_type)
		{
			foreach(explode("&",$variant_type) as $vt)
			{
				//skip other regions
				if(in_array($vt,$miso_terms_coding))	$skip_variant = false;
			}
		}
	}
	
	add_filter($filter, "not_cod_spli", "Not a coding or splicing variant (filter_vcf).");
	if($skip_variant) activate_filter ($filter, "not_cod_spli");
}

function filter_off_target(&$filter, $chr, $start, $tumor_id, $normal_id, $targets)
{
	//filter out bad rows
	$skip_variant = true;
	if(isset($targets[$chr]))
	{
		foreach($targets[$chr] as $regions)
		{
			list($s,$e) = $regions;
			if($e >= $start && $s <= $start )
			{
				$skip_variant = false;
				break;
			}
		}
	}
	add_filter($filter, "off_target", "Variant is off target (filter_vcf).");
	if($skip_variant) activate_filter($filter, "off_target");
}

function filter_somatic(&$filter, $genotype, $tumor, $normal, $type, $alt, $contamination, $min_af, $var_caller)
{
	if(is_null($tumor))	trigger_error("No normal or tumor column.",E_USER_ERROR);
	
	$tumor_only = false;
	if(is_null($normal))	$tumor_only = true;
	
	$td = NULL;
	$tf = NULL;
	$nd = NULL;
	$nf = NULL;
	if($var_caller=="strelka" && $type == "SNV")
	{
		list($td,$tf) = strelka_SNV($genotype,$tumor,$alt);
		list($nd,$nf) = strelka_SNV($genotype,$normal,$alt);
	}
	else if($var_caller=="strelka" && $type == "INDEL")
	{
		list($td,$tf) = strelka_INDEL($genotype,$tumor);
		list($nd,$nf) = strelka_INDEL($genotype,$normal);
	}
	else if($var_caller=="freebayes" && !$tumor_only)
	{
		list($td,$tf) = freebayes($genotype,$tumor);
		list($nd,$nf) = freebayes($genotype,$normal);
	}
	else if($tumor_only)
	{
		list($td,$tf) = freebayes($genotype,$tumor);
	}

	//filter out bad rows
	$skip_variant = false;
	$min_dp = 8;
	$min_tum_af = $min_af;	//min allele frequency in tumor required for reporting
	$noise_nor_af = 0.03;	//allowed max. mutant allele-frequency in normal
	$noise_nor_fa = 6;	//factor for comparison tumor and normal AF
	if($contamination > $noise_nor_af)
	{
		$noise_nor_af += $contamination;
		$noise_nor_fa = 3;
	}
	
	//filter
	add_filter($filter, "som_depth_tum", "Sequencing depth tumor is < $min_dp (filter_vcf).");
	add_filter($filter, "som_depth_nor", "Sequencing depth normal is < $min_dp (filter_vcf).");
	add_filter($filter, "som_all_freq_tum", "Allele frequency tumor is $min_tum_af (filter_vcf).");
	add_filter($filter, "som_lt_3_reads", "Less than 3 supporting tumor reads (filter_vcf).");
	add_filter($filter, "som_all_freq", "Difference in allele frequencies too low (filter_vcf).");
	add_filter($filter, "som_all_freq_nor", "Allele frequency normal too high (filter_vcf).");
	add_filter($filter, "som_tum_loh", "Loss of heterozygosity within tumor tissue (filter_vcf).");
	if($var_caller=="freebayes" && !$tumor_only)	add_filter($filter, "som_af_ratio", "Allele frequency ratio tumor/normal less-equal $noise_nor_fa; removed from list (filter_vcf).");
	if (!$tumor_only && $var_caller=="freebayes")	add_filter($filter, "som_no_depth", "No depth in tumor or normal; removed from list (filter_vcf).");
	
	//depth in reference too low
	if ($td<8)	activate_filter($filter, "som_depth_tum");
	if ($tf<$min_tum_af)	activate_filter($filter, "som_all_freq_tum");
	if ($td*$tf<2.9)	activate_filter($filter, "som_lt_3_reads");
	if(!$tumor_only)
	{
		if ($nd<8)	activate_filter($filter, "som_depth_nor");
		if ($tf<=$noise_nor_fa*$nf)	activate_filter($filter, "som_all_freq_nor");		// hom (WT) in reference, het in tumor (0.5-3.0% ref freq, but six times higher tumor freq)
		if ($nf>0.4 && $nf<0.6 && $tf>0.9)	activate_filter($filter, "som_tum_loh");
		if ($var_caller=="freebayes" && ($nd==0 || $td==0))	activate_filter($filter, "som_no_depth");
		if ($var_caller=="freebayes" && $nf>0 && $tf/$nf<=2.9)	activate_filter($filter, "som_af_ratio");	// tf/nf < 3 should be filtered to remove most of the germline variants
	}
}

function filter_somatic_ds(&$filter, $genotype,$tumor,$normal,$type,$alt,$var_caller)
{
	$tumor_only = true;
	$min_depth = 50;
	if(is_null($normal))	$tumor_only = true;
	
	//determine column indices
	$td = NULL;
	$tf = NULL;
	$nd = NULL;
	$nf = NULL;
	if($var_caller=="strelka" && $type == "SNV")
	{
		list($td,$tf) = strelka_SNV($genotype,$tumor,$alt);
		list($nd,$nf) = strelka_SNV($genotype,$normal,$alt);
	}
	else if($var_caller=="strelka" && $type == "INDEL")
	{
		list($td,$tf) = strelka_INDEL($genotype,$tumor);
		list($nd,$nf) = strelka_INDEL($genotype,$normal);
	}
	else if($var_caller=="freebayes" && !$tumor_only)
	{
		list($td,$tf) = freebayes($genotype,$tumor);
		list($nd,$nf) = freebayes($genotype,$normal);
	}
	else if($tumor_only)
	{
		list($td,$tf) = freebayes($genotype,$tumor);
	}

	//filter out bad rows => keep
	add_filter($filter, "somds_depth_too_low", "Depth of normal sample is < $min_depth (filter_vcf).");
	add_filter($filter, "somds_no_frequencies", "No allele frequency available (filter_vcf).");
	add_filter($filter, "somds_not_somatic", "Allel frequencies do not indicate a somatic variant (filter_vcf).");
	if($var_caller=="freebayes" && !$tumor_only)	add_filter($filter, "som_af_ratio", "Allele frequency ratio tumor/normal <= 3; removed from list (filter_vcf).");
	if (!$tumor_only && $var_caller=="freebayes")	add_filter($filter, "som_no_depth", "No depth in tumor or normal; removed from list (filter_vcf).");
	
	//
	if ($td<$min_depth || (!$tumor_only && $nd<$min_depth))	activate_filter($filter, "somds_depth_too_low");

	//no frequency => keep
	$keep = false;
	if (!$tumor_only && $var_caller=="freebayes" && $nf!=0 && $tf/$nf>3)	activate_filter($filter, "som_af_ratio");
	if (!$tumor_only && $var_caller=="freebayes" && ($nd==0 || $td==0))	activate_filter($filter, "som_no_depth");
	if ($nf=="n/a" || $nf=="n/a")	activate_filter($filter, "somds_no_frequencies");
	// hom (WT) in reference, het in tumor
	else if ($tf>0.03 && $nf<0.005)	$keep = true;
	// hom (WT) in reference, het in tumor (0.5-3.0% ref freq, but six times higher tumor freq)
	else if ($nf>=0.005 && $nf<=0.03 && $tf>=6*$nf)	$keep = true;
	else	activate_filter($filter, "somds_not_somatic");
}

function filter_somatic_capa(&$filter, $info, $genotype, $tumor, $normal, $type, $alt, $var_caller)
{
	$tumor_only = false;
	if(is_null($normal))	$tumor_only = true;

	//
	$min_td = 100;
	$min_nd = 100;
	$max_af = 0.01;
	
	//determine databases
	$tg = 0;
	if(isset($info["T1000GP_AF"]))	$tg = $info["T1000GP_AF"];	//T1000GP_AF=0.0379393;
	$ex = 0;
	if(isset($info["EXAC_AF"]))	$ex = $info["EXAC_AF"];	//EXAC_AF=0.347
	$kv = 0;
	if(isset($info["KAVIAR_AF"]))	$kv = $info["KAVIAR_AF"]; //KAVIAR_AF=0.2553
	
	//determine depth and frequency
	$td = NULL;
	$tf = NULL;
	$nd = NULL;
	$nf = NULL;
	if($var_caller=="strelka" && $type == "SNV")
	{
		list($td,$tf) = strelka_SNV($genotype,$tumor,$alt);
		list($nd,$nf) = strelka_SNV($genotype,$normal,$alt);
	}
	else if($var_caller=="strelka" && $type == "INDEL")
	{
		list($td,$tf) = strelka_INDEL($genotype,$tumor);
		list($nd,$nf) = strelka_INDEL($genotype,$normal);
	}
	else if($var_caller=="freebayes" && !$tumor_only)
	{
		list($td,$tf) = freebayes($genotype,$tumor);
		list($nd,$nf) = freebayes($genotype,$normal);
	}
	else if($tumor_only)
	{
		list($td,$tf) = freebayes($genotype,$tumor);
	}

	//filter
	add_filter($filter, "somca_depth_tum", "Sequencing depth tumor is too low < $min_td (filter_vcf).");
	if (!$tumor_only)	add_filter($filter, "somca_depth_norm", "Sequencing depth normal is too low < $min_nd (filter_vcf).");
	add_filter($filter, "somca_all_freq", "Allele frequencies (filter_vcf).");
	add_filter($filter, "somca_db_frequencies", "Allele frequencies in public databases are above $max_af (filter_vcf).");
	
	//skip variants with high MAF or low depth
	if ($td<$min_td)	activate_filter($filter, "somca_depth_tum");
	if (!$tumor_only && $nd<$min_nd)	activate_filter($filter, "somca_depth_norm");
	if(!$tumor_only &&  ($tf<0.05 || $nf>0.01))	activate_filter($filter, "somca_all_freq");
	if ($tg>$max_af || $ex>$max_af || $kv>$max_af) activate_filter($filter, "somca_db_frequencies");
}

function strelka_SNV($genotype, $column, $alt)
{
	$g = explode(":",$genotype);
	$index_depth = NULL;
	$index_TU = NULL;
	$index_AU = NULL;
	$index_CU = NULL;
	$index_GU = NULL;
	for($i=0;$i<count($g);++$i)
	{
		if($g[$i]=="DP")	$index_depth = $i;
		if($g[$i]=="TU")	$index_TU = $i;
		if($g[$i]=="AU")	$index_AU = $i;
		if($g[$i]=="CU")	$index_CU = $i;
		if($g[$i]=="GU")	$index_GU = $i;
	}
	if(is_null($index_depth) || is_null($index_TU) || is_null($index_AU) || is_null($index_CU) || is_null($index_GU))	trigger_error("Invalid strelka format; either field DP or A/C/G/T not available.",E_USER_ERROR);
	
	$d = explode(":",$column)[$index_depth];
	list($nuc_t,) = explode(",", explode(":",$column)[$index_TU]);
	list($nuc_a,) = explode(",", explode(":",$column)[$index_AU]);
	list($nuc_c,) = explode(",", explode(":",$column)[$index_CU]);
	list($nuc_g,) = explode(",", explode(":",$column)[$index_GU]);
	$o = 0;
	if($alt == "T") $o = $nuc_t;
	if($alt == "A") $o = $nuc_a;
	if($alt == "C") $o = $nuc_c;
	if($alt == "G") $o = $nuc_g;
	if(($nuc_a+$nuc_t+$nuc_c+$nuc_g)==0)	return array($d,"n/a");
	$f = number_format($o/($nuc_a+$nuc_t+$nuc_c+$nuc_g),4);

	return array($d,$f);
}

function strelka_INDEL($genotype, $column)
{
	$g = explode(":",$genotype);

	$index_depth = NULL;
	$index_TIR = NULL;
	$index_TAR = NULL;
	for($i=0;$i<count($g);++$i)
	{
		if($g[$i]=="DP")	$index_depth = $i;
		if($g[$i]=="TIR")	$index_TIR = $i;
		if($g[$i]=="TAR")	$index_TAR = $i;
	}
	
	if(is_null($index_depth) || is_null($index_TIR) || is_null($index_TAR))	trigger_error("Invalid strelka format; either field DP, TIR or TAR not available.",E_USER_ERROR);
	$d = explode(":",$column)[$index_depth];
	list($tir,) = explode(",", explode(":",$column)[$index_TIR]);
	list($tar,) = explode(",", explode(":",$column)[$index_TAR]);

	//tir and tar contain strong supportin reads, tor (not considered here) contains weak supportin reads like breakpoints
	//only strong supporting reads are used for filtering
	$f = "n/a";
	if(($tir+$tar) != 0)	$f = number_format($tir/($tir+$tar),4);
	
	return array($d,$f);
}

function freebayes($genotype, $column)
{
	$g = explode(":",$genotype);
	$index_DP = NULL;
	$index_AO = NULL;
	for($i=0;$i<count($g);++$i)
	{
		if($g[$i]=="DP")	$index_DP = $i;
		if($g[$i]=="AO")	$index_AO = $i;
	}

	if(is_null($index_DP) || is_null($index_AO))	trigger_error("Invalid freebayes format; either field DP or AO not available.",E_USER_ERROR);
	
	$d = explode(":",$column)[$index_DP];
	$f = "n/a";
	if($d>0)	$f = number_format(explode(":",$column)[$index_AO]/$d, 4);
	return array($d,$f);
}

function add_filter(&$filter,$id,$desc)
{
	if(isset($filter[$id]))	trigger_error("Filter $id already declared.",E_USER_ERROR);
	if(strpos($desc,",")!==FALSE || strpos($desc,"=")!==FALSE)	trigger_error("Error adding filter $id. ',=' currently not allowed in filter descriptions.",E_USER_ERROR);
	$filter[$id] = array("desc" => $desc, "active" => false);
}

function activate_filter(&$filter, $id)
{
	$filter[$id]["active"] = true;	
}
