<?php
/**
	@page filter_tsv
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("filter_tsv", "Filter TSV-files according to different filters.");
$parser->addInfile("in",  "Input variant file in TSV format containing all necessary columns (s. below for each filter).", false);
$parser->addOutfile("out",  "Output variant file in TSV format.", false);
$filter = array('somatic', 'somatic_ds', 'coding', 'non_synonymous', 'somatic_diag_capa', 'iVac');
$parser->addString("type", "Filter set to use, can be comma delimited. Valid are: ".implode(",",$filter).".",false);
//optional
$parser->addInfile("roi", "Target region BED file (for off-target filter)", true);
$parser->addFlag("r", "Reduce variants to those passing the filter.");
extract($parser->parse($argv));

$in_file = Matrix::fromTSV($in);
$types = explode(",",$type);
foreach($types as $type)	if(!in_array($type,$filter))	trigger_error("Unknown filter '".$type."' used.",E_USER_ERROR);

//reset filter column
$filter_column = "filter_tsv";
$f = $in_file->getColumnIndex($filter_column, false, false);
$comments = $in_file->getComments();
if($f!==FALSE)
{
	$in_file->removeCol($f);
	$comments = array();
	
	//clean up column descriptions
	foreach($in_file->getComments() as $c)
	{
		if(strpos($c,"#DESCRIPTION=".$filter_column."=")===0)	continue;	//skip filter_tsv description
		$comments[] = $c;
	}
	$in_file->setComments($comments);
}
$comments[] = "#DESCRIPTION=$filter_column=Filter column that contains general filter options; filter(s) used: $type.";
$in_file->addCol(array_fill(0,$in_file->rows(),""),$filter_column);

if(in_array("somatic", $types))
{
	$in_file = filter_coding($in_file);
	$in_file = filter_somatic($in_file);
}
if(in_array("somatic_ds", $types))
{
	$in_file = filter_coding($in_file);
	$in_file = filter_somatic_ds($in_file);
}
if(in_array("coding", $types))
{
	$in_file = filter_coding($in_file);
}
if(in_array("non_synonymous", $types))
{
	$in_file = filter_synonymous($in_file);
}
if(in_array("somatic_diag_capa", $types))
{
	$in_file = filter_coding($in_file);
	$in_file = filter_synonymous($in_file);
	$in_file = filter_somatic($in_file);
	$in_file = filter_somatic_capa($in_file);
	if ($roi!="")
	{
		$in_file = filter_off_target($in_file,$in,$roi);
	}
}
if(in_array("iVac", $types))
{
	$in_file = filter_coding($in_file);
	$in_file = filter_synonymous($in_file);
	$in_file = filter_somatic($in_file);
	if ($roi!="")
	{
		$in_file = filter_off_target($in_file,$in,$roi);
	}
}

$f = $in_file->getColumnIndex($filter_column);
for($i=0;$i<$in_file->rows();++$i)
{
	$row = $in_file->getRow($i);
	$filter = $row[$f];
	if(empty($filter))	$filter = "PASS";
	$in_file->set($i,$f,$filter);
}

//filter
if($r)
{
	$f = $in_file->getColumnIndex($filter_column);
	$tmp_iff = new Matrix();
	$tmp_iff->setHeaders($in_file->getHeaders());
	$tmp_iff->setComments($in_file->getComments());
	for($i=0;$i<$in_file->rows();++$i)
	{
		$row = $in_file->getRow($i);
		if($row[$f] != "PASS")	continue;
		$tmp_iff->addRow($row);
	}
	$in_file = $tmp_iff;
}

//add comments and save files
$in_file->setComments($comments);
$in_file->toTSV($out);

//filter functions
function filter_synonymous(Matrix $data)
{
	
	//determine column indices
	$v = $data->getColumnIndex("variant_type");
	$f = $data->getColumnIndex("filter_tsv");
	
	//coding terms
	$obo = Obo::fromOBO(repository_basedir()."/data/dbs/Ontologies/so-xp_3_0_0.obo");
	$miso_terms_coding = array();
	$ids = array("SO:0001580","SO:0001568");
	foreach($ids as $id)
	{
		$miso_terms[] = $obo->getTermNameByID($id);
		$tmp = $obo->getTermChildrenNames($id,true);
		$miso_terms = array_merge($miso_terms,$tmp);
	}
	$miso_terms_coding = array_unique($miso_terms);
	
	//prepare variant_type filter (valid effects for coding & splicing)
	//children of synonymous_variant (SO:0001819)
	$obo = Obo::fromOBO(repository_basedir()."/data/dbs/Ontologies/so-xp_3_0_0.obo");
	$miso_terms = array();
	$ids = array("SO:0001819");
	foreach($ids as $id)
	{
		$miso_terms[] = $obo->getTermNameByID($id);
		$tmp = $obo->getTermChildrenNames($id,true);
		$miso_terms = array_merge($miso_terms,$tmp);
	}
	$miso_terms = array_unique($miso_terms);

	//filter out bad rows
	$output = new Matrix();
	$output->setHeaders($data->getHeaders());
	for ($i=0; $i<$data->rows(); ++$i)
	{
		$skip_variant = true;
		$row = $data->getRow($i);
		
		$variant_types = explode(",",str_replace("&", ",", $row[$v]));
		foreach($variant_types as $variant_type)
		{
			//skip other regions
			if(!in_array($variant_type,$miso_terms) && in_array($variant_type,$miso_terms_coding))	$skip_variant = false;
		}
		
		$filter = array();
		if($skip_variant)	$filter[]	= "synonymous_variant";	
		$row[$f] = set_filter($filter,$row[$f]);
		
		$output->addRow($row);
	}
	
	return $output;
}


//filter functions
function filter_coding(Matrix $data)
{
	
	//determine column indices
	$v = $data->getColumnIndex("variant_type");
	$f = $data->getColumnIndex("filter_tsv");
	
	//prepare variant_type filter (valid effects for coding & splicing)
	//children of coding_transcript_variant (SO:0001576)
	//=> SO:0001580	coding_sequence_variant
	//=> added SO:0001568 splicing_variant
	//=> skipped SO:0001969	coding_transcript_intron_variant || n.b.
	//=> skipped SO:0001622	UTR_variant
	$obo = Obo::fromOBO(repository_basedir()."/data/dbs/Ontologies/so-xp_3_0_0.obo");
	$miso_terms = array();
	$ids = array("SO:0001580","SO:0001568");
	foreach($ids as $id)
	{
		$miso_terms[] = $obo->getTermNameByID($id);
		$tmp = $obo->getTermChildrenNames($id,true);
		$miso_terms = array_merge($miso_terms,$tmp);
	}
	$miso_terms = array_unique($miso_terms);
	
	//filter out bad rows
	$output = new Matrix();
	$output->setHeaders($data->getHeaders());
	for ($i=0; $i<$data->rows(); ++$i)
	{
		$skip_variant = true;
		$filter = array();
		$row = $data->getRow($i);
		
		$variant_types = explode(",",str_replace("&", ",", $row[$v]));
		foreach($variant_types as $variant_type)
		{
			//skip other regions
			if(in_array($variant_type,$miso_terms))	$skip_variant = false;
		}
		
		if($skip_variant)	$filter[]	= "not_coding_or_splicing";	
		$row[$f] = set_filter($filter,$row[$f]);
		
		$output->addRow($row);
	}
	
	return $output;
}

function filter_off_target(Matrix $data, $in, $roi)
{
	//determine column indices
	$tumor_only = false;
	$tf = $data->getColumnIndex("tumor_af");
	$td = $data->getColumnIndex("tumor_dp");
	$nf = $data->getColumnIndex("normal_af",false,false);
	$nd = $data->getColumnIndex("normal_dp",false,false);
	$f = $data->getColumnIndex("filter_tsv");
	if($nf===FALSE)	$tumor_only = true;
	
	$ps_tid = "";
	$ps_nid = "";
	if($tumor_only)	$ps_tid = basename($in,".GSvar");
	if(!$tumor_only)
	{
		$pair = basename($in,".GSvar");
		list($ps_tid,$ps_nid) = explode("-",$pair);
	}
	
	//compare to bed-file
	$target_bed = Matrix::fromTSV($roi);
	$ts = array();
	for($i=0;$i<$target_bed->rows();++$i)
	{
		$row = $target_bed->getRow($i);
		if(!isset($ts[$row[0]]))	$ts[$row[0]] = array();
		$ts[$row[0]][] = array($row[1]+1,$row[2]);
	}
	
	//filter out bad rows
	$output = new Matrix();
	$output->setHeaders($data->getHeaders());
	for ($i=0; $i<$data->rows(); ++$i)
	{
		$row = $data->getRow($i);
		$filter = array();
		
		$skip_variant = true;
		if(isset($ts[$row[0]]))
		{
			foreach($ts[$row[0]] as $regions)
			{
				list($start,$end) = $regions;
				if($end >= $row[1] && $start <= $row[2] )	$skip_variant = false;
			}
		}
		if($skip_variant)	$filter[] = "off_target";
		$row[$f] = set_filter($filter,$row[$f]);

		$output->addRow($row);
	}

	return $output;
}

function filter_somatic(Matrix $data)
{
	//determine column indices
	$tumor_only = false;
	$tf = $data->getColumnIndex("tumor_af");
	$td = $data->getColumnIndex("tumor_dp");
	$nf = $data->getColumnIndex("normal_af",false,false);
	$nd = $data->getColumnIndex("normal_dp",false,false);
	$f = $data->getColumnIndex("filter_tsv");
	if($nf===FALSE)	$tumor_only = true;

	//filter out bad rows
	$output = new Matrix();
	$output->setHeaders($data->getHeaders());
	for ($i=0; $i<$data->rows(); ++$i)
	{
		$row = $data->getRow($i);
		$skip_variant = false;
		
		//depth in reference too low
		if ($row[$td]<8)	$filter[] = "depth_tumor_too_low";
		if (!$tumor_only && $row[$nd]<8)	$filter[] = "depth_normal_too_low";
		
		$filter = array();
		$min_af = 0.05;	//min allele frequency required for reporting
		if ($row[$tf]<$min_af)	$filter[] = "allele_frequencies_tumor";
		else if ($row[$td]*$row[$tf]<2.9)	$filter[] = "less_than_three_supporting_tumor_reads";
		else if (!$tumor_only && $row[$tf]<$min_af && $row[$nf]>0.005)	$filter[] = "allele_frequencies";		// hom (WT) in reference, het in tumor
		else if (!$tumor_only && $row[$nf]>=0.005 && $row[$nf]<=0.03 && $row[$tf]<=6*$row[$nf])	$filter[] = "allele_frequencies_normal";		// hom (WT) in reference, het in tumor (0.5-3.0% ref freq, but six times higher tumor freq)
		if(!$tumor_only && $row[$nf]>0.4 && $row[$nf]<0.6 && $row[$tf]>0.9)	$filter[] = "tumor_loh";
		
		$row[$f] = set_filter($filter,$row[$f]);

		$output->addRow($row);
	}

	return $output;
}

function filter_somatic_ds(Matrix $data)
{
	//determine column indices
	$filter = array();
	$tf = $data->getColumnIndex("ds_tum_freq");
	$td = $data->getColumnIndex("ds_tum_depth");
	$nf = $data->getColumnIndex("ds_ref_freq");
	$nd = $data->getColumnIndex("ds_ref_depth");
	$f = $data->getColumnIndex("filter_tsv");

	//filter out bad rows
	$output = new Matrix();
	$output->setHeaders($data->getHeaders());
	for ($i=0; $i<$data->rows(); ++$i)
	{
		$row = $data->getRow($i);
		
		//deep sequencing depth too low to analyze => keep
		if ($row[$nd]<50 || $row[$nd]<50)
		{
			$filter[] = "depth_too_low";
		}
		
		//no frequency => keep
		if ($row[$nf]=="n/a" || $row[$nf]=="n/a")
		{
			$filter[] = "no_frequencies";
		}		
		// hom (WT) in reference, het in tumor
		else if ($row[$tf]>0.03 && $row[$nf]<0.005)
		{
			$output->addRow($row);
		}
		// hom (WT) in reference, het in tumor (0.5-3.0% ref freq, but six times higher tumor freq)
		else if ($row[$nf]>=0.005 && $row[$nf]<=0.03 && $row[$tf]>=6*$row[$nf])
		{
			$output->addRow($row);
		}
		else
		{
			$filter[] = "not_somatic";
		}
	}
	
	//set filter row
	$row[$f] = set_filter($filter,$row[$f]);
	$output->addRow($row);	
	return $output;
}

function filter_somatic_capa(Matrix $data)
{
	//
	$min_td = 50;
	$min_nd = 20;
	
	//determine column indices
	$tumor_only = false;
	$tf = $data->getColumnIndex("tumor_af");
	$td = $data->getColumnIndex("tumor_dp");
	$nf = $data->getColumnIndex("normal_af",false,false);
	$nd = $data->getColumnIndex("normal_dp",false,false);
	$tg = $data->getColumnIndex("1000g", true);
	$ex = $data->getColumnIndex("ExAC", true);
	$gn = $data->getColumnIndex("gnomAD", true);
	$f = $data->getColumnIndex("filter_tsv");
	if($nf===FALSE)	$tumor_only = true;
	
	//filter out bad rows
	$output = new Matrix();
	$output->setHeaders($data->getHeaders());
	for ($i=0; $i<$data->rows(); ++$i)
	{
		$row = $data->getRow($i);
		$skip_variant = false;
		$filter = array();
		
		//skip variants with high MAF or low depth
		if ($row[$td]<$min_td)	$filter[] = "depth_tumor_too_low";
		if (!$tumor_only && $row[$nd]<$min_nd)	$filter[] = "depth_normal_too_low";
		if(!$tumor_only &&  ($row[$tf]<0.05 || $row[$nf]>0.01))	$filter[] = "allele_frequencies";
		if ($row[$tg]>0.05 || $row[$ex]>0.05 || $row[$gn]>0.05) $filter[] = "db_frequencies";
		
		$row[$f] = set_filter($filter,$row[$f]);
		$output->addRow($row);
	}
	
	return $output;
}

//skip duplicate filters
function set_filter($filter,$f_cell)
{
	$f_cell = explode(";",$f_cell);
	foreach($filter as $f)
	{
		if(!in_array($f,$f_cell))	$f_cell[] = $f;
	}
	
	return trim(implode(";",$f_cell),";");
}

