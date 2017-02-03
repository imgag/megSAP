<?php 
/** 
	@page vcf2gsvar
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("vcf2gsvar", "Converts an annotated VCF file from freebayes to a GSvar file.");
$parser->addInfile("in",  "Input file in VCF format.", false);
$parser->addOutfile("out", "Output file in GSvar format.", false);
extract($parser->parse($argv));

//helper functions
function extract_numeric($col, $data, $default, $decimal_places, $multi_strategy = "error")
{
	@$output = $data[$col];
	if ($output=="") return $default;
	
	//split remove "." elements
	$parts = explode(",", $output);
	$parts = array_diff($parts, array("."));
	if (count($parts)==0) return $default;
	if ($multi_strategy=="error")
	{
		if (count($parts)>1) trigger_error("Cannot extract numeric value from '$output' in column '$col'!", E_USER_ERROR);
		$output = $parts[0];
	}
	else
	{
		$output = $multi_strategy($parts);
	}
	
	//formatting (round, but indicate that the value is not 0 if it isn't) 
	$formatted = number_format($output, $decimal_places, ".", "");
	if ($output>0 && $formatted==0)
	{
		return substr($formatted, 0, -1)."1";
	}
	else
	{
		return $formatted;
	}
}

function extract_string($col, $data, $default)
{
	$output = ".";
	if (isset($data[$col])) $output = $data[$col];
	if ($output==".") $output = $default;
	
	return $output;
}

//determines if all the input genes are on the blacklist
function all_genes_blacklisted($genes)
{
	//init blacklist on first call
	static $blacklist = null;
	if ($blacklist === null)
	{
		$file = file(get_path("data_folder")."/gene_blacklist/blacklist.tsv");
		foreach($file as $line)
		{
			$line = trim($line);
			if ($line=="" || $line[0]=="#") continue;
			list($gene) = explode("\t", $line);
			$blacklist[$gene] = true;
		}
	}
  
	if (count($genes)==0)
	{
		return false;
	}
	
	foreach ($genes as $gene)
	{
		if (!isset($blacklist[$gene]))
		{
			return false;
		}
	}
	
	return true;
}

//init
$handle_out = fopen($out, "w");

//write column descriptions
$column_desc = array(
	array("filter", "Annotations for filtering and ranking variants."),
	array("quality", "Quality parameters - SNP quality (QUAL), depth (DP), allele frequency (AF), mean mapping quality of alternate allele (MQM), genotype, depth and allele frequency for child, mother and father (TRIO)."),
	array("gene", "Affected gene list (comma-separated)."),
	array("variant_type", "Variant type."),
	array("coding_and_splicing", "Coding and splicing details (Gene, NM number, type, impact, exon number, HGVS.c, HGVS.p)."),
	array("RepeatMasker", "RepeatMasker annotation."),
	array("dbSNP", "Identifier in dbSNP database."),
	array("1000g", "Allele frequency in all populations of 1000g project."),
	array("ExAC", "Allele frequency in all populations of ExAC project."),
	array("ExAC_hom", "Homoyzgous counts for populations ALL, NFE and AFR of ExAC project."),
	array("Kaviar", "Allele frequency in Kaviar database."),
	array("phyloP", "phyloP (100way vertebrate) annotation. Deleterious threshold > 1.6."),
	array("MetaLR", "MetaLR effect prediction: D=damaging, T=tolerated."),
	array("Sift", "Sift effect prediction: D=damaging, T=tolerated."),
	array("PP2_HVAR", "Polyphen2 (HVAR) effect prediction: D=probably damaging, P=possibly damaging, B=benign."),
	array("PP2_HDIV", "Polyphen2 (HDIV) effect prediction: D=probably damaging, P=possibly damaging, B=benign."),
	array("OMIM", "OMIM database annotation."),
	array("ClinVar", "ClinVar database annotation."),
	array("HGMD", "HGMD database annotation."),
	array("COSMIC", "COSMIC somatic variant database anntotation."),
);
foreach($column_desc as $entry)
{
	fwrite($handle_out, "##DESCRIPTION=".$entry[0]."=".$entry[1]."\n");
}
function write_header_line($handle, $column_desc)
{
	fwrite($handle, "#chr\tstart\tend\tref\tobs\tgenotype");
	foreach($column_desc as $entry)
	{
		fwrite($handle, "\t".$entry[0]);
	}
	fwrite($handle, "\n");
}

//write filter descriptions
$filter_desc = array(
	array("low_DP", "Depth less than 20."),
	array("low_MQM", "Mean mapping quality of alternate allele less than Q50."),
	array("low_QUAL", "Variant quality less than Q30."),
	array("phylop_high", "PhyloP value higher than 1.6 i.e. highly-conserved base and possibly pathognic."),
	array("pred_pathogenic", "Variant predicted to be pathogenic by one or more prediction tools."),
	array("gene_blacklist", "The gene(s) are contained on the blacklist of unreliable genes."),
	array("anno_pathogenic_clinvar", "Variant annotated to be pathogenic by ClinVar."),
	array("anno_pathogenic_hgmd", "Variant annotated to be pathogenic by HGMD."),
	array("anno_high_impact", "Variant annotated to have high impact by SnpEff."),
);
foreach($filter_desc as $entry)
{
	fwrite($handle_out, "##FILTER=".$entry[0]."=".$entry[1]."\n");
}

//parse input
$in_header = true;
$handle = fopen($in, "r");
while(!feof($handle))
{
	$line = nl_trim(fgets($handle));
	if ($line=="") continue;
	
	//write filter descriptions
	if ($line[0]=="#") 
	{
		if (starts_with($line, "##FILTER=<ID="))
		{
			$parts = explode(",Description=\"", substr(trim($line), 13, -2));
			fwrite($handle_out, "##FILTER=".$parts[0]."=".$parts[1]."\n");
		}
		continue;
	}
	//after last header line, write our header
	else if ($in_header) 
	{
		write_header_line($handle_out, $column_desc);
		$in_header = false;
	}
	
	//write content lines
	$cols = explode("\t", $line);
	if (count($cols)<9) trigger_error("VCF file line contains less than 10 columns:\n$line", E_USER_ERROR);
	list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $sample) = $cols;
	if ($filter=="" || $filter=="." || $filter=="PASS")
	{
		$filter = array();
	}
	else
	{
		$filter = explode(";", $filter);
	}
	
	//parse data from VCF
	if(chr_check($chr, 22, false) === FALSE) continue; //skip bad chromosomes
	$start = $pos;
	$end = $pos;
	$ref = strtoupper($ref);
	$alt = strtoupper($alt);
	if(strlen($ref)>1 || strlen($alt)>1) //correct indels
	{
		list($start, $end, $ref, $alt) = correct_indel($start, $ref, $alt);
	}
	$info = explode(";", $info);
	$tmp = array();
	foreach($info as $entry)
	{
		if (!contains($entry, "="))
		{
			$tmp[$entry] = "";
		}
		else
		{
			list($key, $value) = explode("=", $entry, 2);
			$tmp[$key] = $value;
		}
	}
	$info = $tmp;
	$sample = array_combine(explode(":", $format), explode(":", $sample));
	
	//convert genotype to TSV format
	if (!isset($sample["GT"])) trigger_error("VCF sample column does not contain GT value!", E_USER_ERROR);
	$genotype = strtr($sample["GT"], array("."=>"0", "/"=>"|"));
	//skip wildtype
	if ($genotype=="0|0") continue;			
	//convert genotype to TSV format
	else if ($genotype=="0|1" || $genotype=="1|0") $genotype = "het";
	else if ($genotype=="1|1") $genotype = "hom";
	else trigger_error("Unknown VCF genotype value '$genotype' in line '$line'.", E_USER_ERROR);

	//quality
	$quality = array();
	$qual = intval($qual);
	$quality[] = "QUAL=".$qual;
	if ($qual<30) $filter[] = "low_QUAL";
	if (isset($sample["DP"]))
	{
		$quality[] = "DP=".$sample["DP"];
		if ($sample["DP"]<20) $filter[] = "low_DP";
	}
	if (isset($sample["AO"]) && isset($sample["DP"]))
	{
		//@todo check regularly if this bug is fixed in vcflib. Example: NA12878_03 chr2:215632255 (AF~0.5) and chr2:215632256 (AF~1.0)
		/* 
		if (contains($sample["AO"], ",")) 
		{
			$sample["AO"] = array_sum(explode(",", $sample["AO"]));
		}
		$quality[] = "AF=".number_format(min(1.0, $sample["AO"]/$sample["DP"]), 2);
		*/
		$quality[] = "AF=".number_format($sample["AO"]/$sample["DP"], 2);
	}
	if (isset($info["MQM"])) 
	{
		$quality[] = "MQM=".intval($info["MQM"]);
		if ($info["MQM"]<50) $filter[] = "low_MQM";
	}
	if (isset($sample["TRIO"]))
	{
		$quality[] = "TRIO=".$sample["TRIO"];
	}
	
	//variant details
	//ANN field: Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO
	$genes = array();
	$variant_details = array();
	$coding_and_splicing_details = array();
	$ann = extract_string("ANN", $info, "");
	if ($ann!="")
	{
		$anns = explode(",", $ann);
		foreach($anns as $entry)
		{
			$parts = explode("|", $entry);
			$details = strtr($parts[1], array("_variant"=>""));
			$details = strtr($details, array("splice_acceptor&splice_region&intron"=>"splice_acceptor", "splice_donor&splice_region&intron"=>"splice_donor", "splice_acceptor&intron"=>"splice_acceptor", "splice_donor&intron"=>"splice_donor", "_prime_"=>"'"));
			if ($details=="intragenic") continue; //skip these details
			$variant_details[] = $details;
			$genes[] = $parts[3];
			$exon = $parts[8];
			if ($exon!="") $exon = "exon".$exon;
			$coding_and_splicing_details[] = $parts[3].":".$parts[6].":$details:".$parts[2].":$exon:".$parts[9].":".$parts[10];
		}
	}
	$genes = array_unique($genes);
	if (all_genes_blacklisted($genes))
	{
		$filter[] = "gene_blacklist";
	}
	$variant_details = implode(",", array_unique($variant_details));
	$coding_and_splicing_details =  implode(",", $coding_and_splicing_details);
	if(contains($coding_and_splicing_details, ":HIGH:")) $filter[] = "anno_high_impact";

	//RepeatMasker
	$repeatmasker = extract_string("REPEATMASKER", $info, "");
	
	//Frequency databases
	$dbsnp = "";
	if ($id!=".") $dbsnp = $id;
	$kgenomes = extract_numeric("T1000GP_AF", $info, "0.0000", 4, "max");
	$kaviar = extract_numeric("KAVIAR_AF", $info, "0.0000", 4);
	$exac = extract_numeric("EXAC_AF", $info, "0.0000", 4);
	$exac_hom_all = extract_numeric("EXAC_AC_Hom", $info, "0", 0);
	$exac_hom_nfe = extract_numeric("EXAC_Hom_NFE", $info, "0", 0);
	$exac_hom_afr = extract_numeric("EXAC_Hom_AFR", $info, "0", 0);
	$exac_hom = "$exac_hom_all/$exac_hom_nfe/$exac_hom_afr";
	//effect predicions
	$phylop = extract_numeric("dbNSFP_phyloP100way_vertebrate", $info, "", 4, "max");
	if ($phylop>=1.6)
	{
		$filter[] = "phylop_high";
	}
	$metalr = extract_string("dbNSFP_MetaLR_pred", $info, "");
	$sift = extract_string("dbNSFP_SIFT_pred", $info, "");
	$pp2v = extract_string("dbNSFP_Polyphen2_HDIV_pred", $info, "");
	$pp2d = extract_string("dbNSFP_Polyphen2_HVAR_pred", $info, "");
	$pp_count = contains($metalr, "D") + contains($sift, "D") + contains($pp2v, "D") + contains($pp2d, "D");
	if ($pp_count>0)
	{
		$filter[] = "pred_pathogenic";
	}
	
	//OMIM
	$omim = strtr(extract_string("OMIM", $info, ""), "_", " ");
	if ($omim!="") $omim .= ";";
	
	//ClinVar
	$clin_acc = explode("|", extract_string("CLINVAR_ACC", $info, ""));
	$clin_sig = explode("|", extract_string("CLINVAR_SIG", $info, ""));
	if (count($clin_acc)!=count($clin_sig)) trigger_error("Clinvar field counts do not match:\n".implode("|",$clin_acc)."\n".implode("|",$clin_sig)."" , E_USER_ERROR);
	$clinvar = "";
	for($i=0; $i<count($clin_acc); ++$i)
	{
		if (trim($clin_acc[$i]=="")) continue;
		$clinvar .= $clin_acc[$i]." [".strtr($clin_sig[$i], "_", " ")."]; ";
	}
	if (contains($clinvar, "pathogenic")) $filter[] = "anno_pathogenic_clinvar";  //matches "pathogenic" and "likely pathogenic"
	
	//HGMD
	$hgmd_id = explode("|", extract_string("HGMD_ID", $info, ""));
	$hgmd_class = explode("|", extract_string("HGMD_CLASS", $info, ""));
	$hgmd_mut = explode("|", extract_string("HGMD_MUT", $info, ""));
	$hgmd_gene = explode("|", extract_string("HGMD_GENE", $info, ""));
	$hgmd_phen = explode("|", extract_string("HGMD_PHEN", $info, ""));
	if (count($hgmd_id)!=count($hgmd_class) || count($hgmd_id)!=count($hgmd_mut) || count($hgmd_id)!=count($hgmd_phen) || count($hgmd_id)!=count($hgmd_gene)) trigger_error("HGMD filed counts do not match:\n".implode("|",$hgmd_id)."\n".implode("|",$hgmd_class)."\n".implode("|",$hgmd_mut)."\n".implode("|",$hgmd_gene)."\n".implode("|",$hgmd_phen)."" , E_USER_ERROR);
	$hgmd = "";
	for($i=0; $i<count($hgmd_id); ++$i)
	{
		if (trim($hgmd_id[$i]=="")) continue;
		$hgmd .= $hgmd_id[$i]." [CLASS=".$hgmd_class[$i]." MUT=".$hgmd_mut[$i]." PHEN=".strtr($hgmd_phen[$i], "_", " ")." GENE=".$hgmd_gene[$i]."]; ";
	}
	if (contains($hgmd, "CLASS=DM")) $filter[] = "anno_pathogenic_hgmd"; //matches both "DM" and "DM?"
		
	//COSMIC
	$cosmic = strtr(extract_string("COSMIC_ID", $info, ""), array(","=>", "));

	//write data
	fwrite($handle_out, "$chr\t$start\t$end\t$ref\t$alt\t$genotype\t".implode(";", $filter)."\t".implode(";", $quality)."\t".implode(",", $genes)."\t$variant_details\t$coding_and_splicing_details\t$repeatmasker\t$dbsnp\t$kgenomes\t$exac\t$exac_hom\t$kaviar\t$phylop\t$metalr\t$sift\t$pp2v\t$pp2d\t$omim\t$clinvar\t$hgmd\t$cosmic\n");
}

//if no variants are present, we need to write the header line after the loop
if ($in_header) 
{
	write_header_line($handle_out, $column_desc);
}

fclose($handle);
fclose($handle_out);

?>
