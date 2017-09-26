<?php
/**
	@page filter_vcf
	@todo revisit somatic and somatic_ds (freebayes!) filter
	@todo rewrite - make filters selectable to the outside - not_off_target, not_coding_splicing, not_synonymous
	@todo additional parameters: min_tumor_dp, min_normal_dp, min_var_reads, min_var_af, conatmination, max_nor_af, max_var_pf
	@todo refactor somatic, somatic_ds, somatic_capa; add filter public db; add filter contamination (if > 0)
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("filter_vcf", "Filter VCF-files according to different filter criteria. This tool is designed to filter tumor/normal pair samples. This tools automatically chooses variant caller and tumor/normal sample from the vcf header.");
$parser->addInfile("in", "Input variant file in VCF format containing all necessary columns (s. below for each filter).", false);
$parser->addOutfile("out", "Output variant file in VCF format.", false);
$filter = array('not-coding-splicing', 'synonymous', 'off-target', 'somatic-lq');
$parser->addString("type", "Type(s) of variants that are supposed to be filtered out, can be comma delimited. Valid are: ".implode(",",$filter).".",false);
//optional
$parser->addFlag("keep", "Keep all variants. Otherwise only variants passing all filters will be kept");
$parser->addFlag("ignore_filter", "Ignore previous filter.");
$parser->addString("roi", "Target region BED file (is required by off-target filter).", true, "");
$parser->addFloat("contamination", "Estimated fraction of tumor cells in normal sample.",true,0.00);
$parser->addFloat("min_af", "Minimum variant allele frequency in tumor.",true,0.05);
// TODO add filter for min_af, min_dp, min_r, db_af, noise nor/cont / max af nor to remove somatic_XYZ filter sets
//$parser->addFloat("min_dp", "Minimum depth.",true,8);
//$parser->addFloat("min_r", "Minimum number reads.",true,3);
extract($parser->parse($argv));

//zipped
$zipped = false;
if(ends_with($in, ".gz"))	$zipped = true;	

//check selected filter types
$types = explode(",",$type);
foreach($types as $t)	if(!in_array($t,$filter) && $t!="all")	trigger_error("Unknown filter '".$t."' given.",E_USER_ERROR);
if(in_array("all",$types))	$types = $filter;
$in_file = Matrix::fromTSV($in);

// check filter column
$filter_column = "FILTER";
$fc = $in_file->getColumnIndex($filter_column, false, false);
if($fc!=6)	trigger_error("Wrong column index for filter column! Is ".(emtpy($fc)?"empty":$fc)." should be 6.",E_USER_ERROR);
$prefix = "";
$comments = $in_file->getComments();
$var_caller = NULL;

//extract tumor and normal column from pedigree tag
//e.g. ##PEDIGREE=<Tumor_DNA=GS0-9,Normal_DNA=GS0-9>
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

	if(strpos($comment,"#source=")===0)	//information about variant caller
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
if(is_null($normal_col))	trigger_error("No normal column given!",E_USER_NOTICE);
if(is_null($var_caller))	trigger_error("Variant caller not identified.",E_USER_ERROR);

//extract relevant MISO terms for filtering
$obo = Obo::fromOBO(repository_basedir()."/data/dbs/Ontologies/so-xp_3_0_0.obo");
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
$obo = Obo::fromOBO(repository_basedir()."/data/dbs/Ontologies/so-xp_3_0_0.obo");
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
if($roi!="") //skip enrichments without target
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
	$tmp_col_tum = $row[$tumor_col];
	$tmp_col_nor = null;
	if(!is_null($normal_col))	$tmp_col_nor = $row[$normal_col];
	//filter_basic($filter,$info,$tmp_col_tum,$tmp_col_nor$min_af,$min_dp,$min_r,$con,$var_caller);
	if(in_array("not-coding-splicing",$types))	filter_not_coding_splicing($filter, $info, $miso_terms_coding);
	if(in_array("synonymous",$types))	filter_synonymous($filter, $info, $miso_terms_coding, $miso_terms_synonymous);
	if(in_array("set-somatic",$types))	filter_somatic($filter, $genotype,$tmp_col_tum,$tmp_col_nor,$type,$row[4],$contamination, $min_af, $var_caller);
	if(in_array("off-target",$types))	filter_off_target($filter, $row[0], $row[1], $tumor_id, $normal_id, $targets);
	if(in_array("somatic-lq",$types))	filter_somatic_lq($filter, $info, $genotype, $tmp_col_tum, $tmp_col_nor, $type,$row[4],$var_caller,100);
		
	//set PASS criterion if all filters were passed
	$tmp_filter = array();
	foreach($filter as $k => $f)
	{
		if(!isset($f["active"]) || !$f["active"])	continue;
		$tmp_filter[] = $prefix.$k;
	}
	$tmp_filter_old = array_filter(explode(";",trim($row[6],".")));
	if(isset($tmp_filter_old[0]) && $tmp_filter_old[0]=="PASS")	unset($tmp_filter_old[0]);	
	
	// check if data was already annotated by filter_vcf or same filter
	if(count($dup = array_intersect($tmp_filter_old,$tmp_filter))>0)	trigger_error("Found '".implode(", ",$dup)."' filter multiple times. Looks like this files was already annotated by filter_vcf.",E_USER_ERROR);
	
	$tmp_filter = array_unique(array_merge($tmp_filter,$tmp_filter_old));	//keep original filter
	
	if(empty($tmp_filter))	$tmp_filter[] = "PASS";
	$in_file->set($i,6,implode(";",	$tmp_filter));
}

//add comments and save files
//##FILTER=<ID=ID,Description="description">
$comments = $in_file->getComments();
foreach($filter as $id => $details)	//every line has the same filter criteria, so filter criteria of the last row is sufficient
{
	$comments[] = "#FILTER=<ID=".$prefix.$id.",Description=\"".$details["desc"]."\">";
}
$in_file->setComments($comments);

// freebayes: remove variants that are likely germline variants
if($var_caller=="freebayes")
{
	trigger_error("Variants were called with freebayes. Candidate germline variants will be removed.",E_USER_NOTICE);
	$filtered = new Matrix();
	$filtered->setComments($in_file->getComments());
	$filtered->setHeaders($in_file->getHeaders());
	$f = $in_file->getColumnIndex("FILTER", false, false);
	for($i=0;$i<$in_file->rows();++$i)
	{
		$row = $in_file->getRow($i);
		if(strpos($row[$f],$prefix."som-af-ratio")!==FALSE)	continue;
		if(strpos($row[$f],$prefix."som-no-depth")!==FALSE)	continue;
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
	
	add_filter($filter, "syn-var", "Synonymous variant (filter_vcf).");
	if($coding && $synonymous) activate_filter($filter, "syn-var");
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
	
	add_filter($filter, "not-cod-spli", "Not a coding or splicing variant (filter_vcf).");
	if($skip_variant) activate_filter ($filter, "not-cod-spli");
}

function filter_off_target(&$filter, $chr, $start, $tumor_id, $normal_id, $targets)
{
	//filter out bad rows
	if(!empty($targets))
	{
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
		add_filter($filter, "off-target", "Variant is off target (filter_vcf).");
		if($skip_variant) activate_filter($filter, "off-target");
	}
	else	trigger_error("Cannot use off-target filter without targets.",E_USER_ERROR);
}

function filter_somatic_lq(&$filter, $info, $genotype, $tumor, $normal, $type, $alt, $var_caller, $min_dp = 100)
{
	$tumor_only = false;
	if(is_null($normal))	$tumor_only = true;

	//
	$min_td = 20;
	$min_nd = 20;
	$min_taf = 0.05;
	
	$max_af = 0.01;
	$max_naf = 0.01;
	
	//determine databases
	$tg = 0;
	if(isset($info["T1000GP_AF"]))	$tg = $info["T1000GP_AF"];	//T1000GP_AF=0.0379393;
	$ex = 0;
	if(isset($info["EXAC_AF"]))	$ex = $info["EXAC_AF"];	//EXAC_AF=0.347
	$gn = 0;
	if(isset($info["GNOMAD_AF"]))	$gn = $info["GNOMAD_AF"]; //GNOMAD_AF=0.2553
	
	//determine depth and frequency
	$td = NULL;
	$tf = NULL;
	$nd = NULL;
	$nf = NULL;
	if($var_caller=="strelka" && $type == "SNV")
	{
		list($td,$tf) = vcf_strelka_snv($genotype,$tumor,$alt);
		list($nd,$nf) = vcf_strelka_snv($genotype,$normal,$alt);
	}
	else if($var_caller=="strelka" && $type == "INDEL")
	{
		list($td,$tf) = vcf_strelka_indel($genotype,$tumor);
		list($nd,$nf) = vcf_strelka_indel($genotype,$normal);
	}
	else if($var_caller=="freebayes" && !$tumor_only)
	{
		list($td,$tf) = vcf_freebayes($genotype,$tumor);
		list($nd,$nf) = vcf_freebayes($genotype,$normal);
	}
	else if($tumor_only)
	{
		list($td,$tf) = vcf_freebayes($genotype,$tumor);
	}

	//filter
	add_filter($filter, "som-lt-3-reads", "Less than 3 supporting tumor reads (filter_vcf).");
	add_filter($filter, "som-depth-tum", "Sequencing depth tumor is too low < $min_td (filter_vcf).");
	if (!$tumor_only)	add_filter($filter, "som-depth-nor", "Sequencing depth normal is too low < $min_nd (filter_vcf).");
	add_filter($filter, "som-freq-tum", "Allele frequencies in tumor < $min_taf or allele frequencies in normal > $max_naf (filter_vcf).");
	if(!$tumor_only)	add_filter($filter, "som-freq-nor", "Allele frequencies in tumor < $min_taf or allele frequencies in normal > $max_naf (filter_vcf).");
	// freebayes add best practice filter to identify germline variants
	if($var_caller=="freebayes" && !$tumor_only)	add_filter($filter, "som-af-ratio", "Allele frequency ratio tumor/normal less equal 2.9; removed from list (filter_vcf).");
	if ($var_caller=="freebayes" && !$tumor_only)	add_filter($filter, "som-no-depth", "No depth in tumor or normal; removed from list (filter_vcf).");
	// tumor only use db filter
	if($tumor_only)	add_filter($filter, "som-db-frequencies", "Allele frequencies in public databases are > $max_af (filter_vcf).");
	
	//skip variants with high MAF or low depth
	if ($td*$tf<2.9)	activate_filter($filter, "som-lt-3-reads");
	if ($td<$min_td)	activate_filter($filter, "som-depth-tum");
	if (!$tumor_only && $nd<$min_nd)	activate_filter($filter, "som-depth-nor");
	if($tf<0.05)	activate_filter($filter, "som-freq-tum");
	if(!$tumor_only &&  $nf>1/6*$tf)	activate_filter($filter, "som-freq-nor");
	//
	if ($tumor_only && ($tg>$max_af || $ex>$max_af || $gn>$max_af)) activate_filter($filter, "somca-db-frequencies");
	if ($var_caller=="freebayes" && ($nd==0 || $td==0))	activate_filter($filter, "som-no-depth");
	if ($var_caller=="freebayes" && $nf>0 && $tf/$nf<=2.9)	activate_filter($filter, "som-af-ratio");	// tf/nf < 3 should be filtered to remove most of the germline variants
}


function add_filter(&$filter,$id,$desc)
{
	if(isset($filter[$id]))	trigger_error("Filter $id already declared.",E_USER_ERROR);
	if(strpos($desc,",")!==FALSE || strpos($desc,"=")!==FALSE)	trigger_error("Error adding filter $id. ',=' currently not allowed in filter descriptions.",E_USER_ERROR);
	$filter[$id] = array("desc" => $desc, "active" => false);
}

function activate_filter(&$filter, $id)
{
	if(!isset($filter[$id]))	trigger_error("Unknown filter $id.",E_USER_ERROR);
	$filter[$id]["active"] = true;	
}
