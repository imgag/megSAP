<?php 
/** 
	@page vcf2gsvar
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("vcf2gsvar", "Converts an annotated VCF file from freebayes to a GSvar file.");
$parser->addInfile("in",  "Input file in VCF format.", false);
$parser->addOutfile("out", "Output file in GSvar format.", false);
$parser->addString("build", "The genome build to use.", false);
//optional
$parser->addFlag("multi", "Enable multi-sample mode.");
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
		$file = file(repository_basedir()."/data/gene_lists/blacklist.tsv");
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


//write header line
function write_header_line($handle, $column_desc, $filter_desc)
{
	foreach($column_desc as $entry)
	{
		fwrite($handle, "##DESCRIPTION=".$entry[0]."=".$entry[1]."\n");
	}
	
	foreach($filter_desc as $entry)
	{
		fwrite($handle, "##FILTER=".$entry[0]."=".$entry[1]."\n");
	}

	fwrite($handle, "#chr\tstart\tend\tref\tobs");
	foreach($column_desc as $entry)
	{
		fwrite($handle, "\t".$entry[0]);
	}
	fwrite($handle, "\n");
}

//write column descriptions
$column_desc = array(
	array("filter", "Annotations for filtering and ranking variants."),
	array("quality", "Quality parameters - SNP quality (QUAL), depth (DP), allele frequency (AF), mean mapping quality of alternate allele (MQM)."),
	array("gene", "Affected gene list (comma-separated)."),
	array("variant_type", "Variant type."),
	array("coding_and_splicing", "Coding and splicing details (Gene, ENST number, type, impact, exon number, HGVS.c, HGVS.p)."),
	array("RepeatMasker", "RepeatMasker annotation."),
	array("dbSNP", "Identifier in dbSNP database."),
	array("1000g", "Allele frequency in all populations of 1000g project."),
	array("ExAC", "Allele frequency in all populations of ExAC project."),
	array("ExAC_hom", "Homoyzgous counts for populations ALL, NFE and AFR of ExAC project."),
	array("ExAC_sub", "Sub-population allele frequenciens for populations AFR,AMR,EAS,NFE,SAS of ExAC project."),
	array("gnomAD", "Allele frequency in gnomAD database."),
	array("phyloP", "phyloP (100way vertebrate) annotation. Deleterious threshold > 1.6."),
	array("Sift", "Sift effect prediction: D=damaging, T=tolerated."),
	array("MetaLR", "MetaLR effect prediction: D=damaging, T=tolerated."),
	array("PolyPhen2", "PolyPhen2 (HVAR) effect prediction: D=probably damaging, P=possibly damaging, B=benign."),
	array("FATHMM", "FATHMM effect prediction: D=damaging, T=tolerated."),
	array("CADD", "CADD pathogenicity prediction scores (scaled phred-like). Deleterious threshold > 15-20."),
	array("OMIM", "OMIM database annotation."),
	array("ClinVar", "ClinVar database annotation."),
	array("HGMD", "HGMD database annotation."),
	array("COSMIC", "COSMIC somatic variant database anntotation."),
);
if (!$multi)
{
	array_unshift($column_desc, array("genotype", "Genotype of variant in sample."));	
}

//write filter descriptions
$filter_desc = array(
	array("low_DP", "Depth less than 20."),
	array("low_MQM", "Mean mapping quality of alternate allele less than Q50."),
	array("low_QUAL", "Variant quality less than Q30."),
	array("pred_pathogenic", "Variant predicted to be pathogenic by one or more tools (conservation or effect prediction)."),
	array("pred_pathogenic_3", "Variant predicted to be pathogenic by three or more tools (conservation or effect prediction)."),
	array("gene_blacklist", "The gene(s) are contained on the blacklist of unreliable genes."),
	array("anno_pathogenic_clinvar", "Variant annotated to be pathogenic by ClinVar."),
	array("anno_pathogenic_hgmd", "Variant annotated to be pathogenic by HGMD."),
	array("anno_high_impact", "Variant annotated to have high impact by SnpEff."),
	array("anno_omim", "Variant annotated with information from OMIM."),
);

//load gencode basic transcripts
$gencode_basic = array();
if ($build=="GRCh37")
{
	$file = file(repository_basedir()."/data/dbs/Ensembl/gencode_basic.txt");
	foreach($file as $line)
	{
		$line = trim($line);
		if ($line=="" || $line[0]=="#") continue;
		$gencode_basic[$line] = true;
	}
}
else
{
	trigger_error("No gencode basic transcripts available for genome '$build', skipping the gencode basic filter.", E_USER_WARNING);
}

//parse input
$multi_cols = array();
$in_header = true;
$handle = fopen($in, "r");
$handle_out = fopen($out, "w");
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
		
		if (starts_with($line, "##SAMPLE="))
		{
			$line = trim($line);
			fwrite($handle_out, $line."\n");
			if ($multi)
			{
				list($name) = explode(",", substr($line, 13, -1));
				$multi_cols[] = $name;
				array_splice($column_desc, count($multi_cols)-1, 0, array(array($name, "genotype of sample $name")));
			}
		}
		
		if (starts_with($line, "##ANALYSISTYPE="))
		{
			fwrite($handle_out, trim($line)."\n");
		}
		
		continue;
	}
	//after last header line, write our header
	else if ($in_header) 
	{
		write_header_line($handle_out, $column_desc, $filter_desc);
		$in_header = false;
	}
	
	//write content lines
	$cols = explode("\t", $line);
	if (count($cols)<10) trigger_error("VCF file line contains less than 10 columns:\n$line", E_USER_ERROR);
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
	
	//convert genotype information to TSV format
	if ($multi)
	{
		if (!isset($sample["MULTI"])) 
		{
			trigger_error("VCF sample column does not contain MULTI value!", E_USER_ERROR);
		}
		
		//extract GT/DP/AO info
		$tmp = array();
		$tmp2 = array();
		$tmp3 = array();
		$parts = explode(",", $sample["MULTI"]);
		foreach($parts as $part)
		{
			list($name, $gt, $dp, $ao) = explode("=", strtr($part, "|", "=")."=");
			$tmp[$name] = $gt;
			$tmp2[$name] = $dp;
			$tmp3[$name] = $ao;
		}
		
		//recombine GT/DP/AO in the correct order
		$genotype = array();
		$depth = array();
		$ao = array();
		foreach($multi_cols as $col)
		{
			$genotype[] = $tmp[$col];
			$depth[] = $tmp2[$col];
			$ao[] = $tmp3[$col];
		}
		$genotype = implode("\t", $genotype);
		$sample["DP"] = implode(",", $depth);
		$sample["AO"] = implode(",", $ao);
	}
	else
	{
		if (!isset($sample["GT"])) 
		{
			trigger_error("VCF sample column does not contain GT value!", E_USER_ERROR);
		}
		$genotype = vcfgeno2human($sample["GT"]);
		
		//skip wildtype
		if ($genotype=="wt") continue;
	}

	//quality
	$quality = array();
	$qual = intval($qual);
	$quality[] = "QUAL=".$qual;
	if ($qual<30) $filter[] = "low_QUAL";
	if (isset($sample["DP"]))
	{
		$quality[] = "DP=".$sample["DP"];
		if (min(explode(",", $sample["DP"]))<20) //comma-separated values in case of multi-sample data
		{
			$filter[] = "low_DP";
		}
	}
	if (isset($sample["AO"]) && isset($sample["DP"]))
	{
		//comma-separated values in case of multi-sample data
		$afs = array();
		$aos = explode(",", $sample["AO"]);
		$dps = explode(",", $sample["DP"]);
		for($i=0; $i<count($dps); ++$i)
		{
			if (is_numeric($aos[$i]) && is_numeric($dps[$i]))
			{
				$afs[] = number_format($aos[$i]/$dps[$i], 2);
			}
		}
		if (count($afs)>0)
		{
			$quality[] = "AF=".implode(",", $afs);
		}
	}
	if (isset($info["MQM"])) 
	{
		$quality[] = "MQM=".intval($info["MQM"]);
		if ($info["MQM"]<50) $filter[] = "low_MQM";
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
			
			//skip sequence_feature
			if ($details=="sequence_feature") continue; 

			//skip empty gene names entries (TF-binding site, etc)
			$gene = trim($parts[3]);
			if ($gene=="") continue; 
			
			//split genes
			if ($details=="intergenic_region")
			{
				$gene_list = explode("-", $gene);
				foreach($gene_list as $gene)
				{
					$genes[] = $gene;
				}
			}
			else
			{
				$genes[] = $gene;
			}
			
			//skip transcripts that are not flagged as "gencode basic"
			if (count($gencode_basic)>0 && $details!="intergenic_region")
			{
				$transcript = trim($parts[6]);
				if (!isset($gencode_basic[$transcript])) continue;
			}

			$variant_details[] = $details;
			$exon = $parts[8];
			if ($exon!="") $exon = "exon".$exon;
			$coding_and_splicing_details[] = $gene.":".$parts[6].":$details:".$parts[2].":$exon:".$parts[9].":".$parts[10];
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
	$gnomad = max(extract_numeric("GNOMAD_AF", $info, "0.0000", 4), extract_numeric("GNOMAD_GENOME_AF", $info, "0.0000", 4));
	$exac = extract_numeric("EXAC_AF", $info, "0.0000", 4);
	$exac_hom_all = extract_numeric("EXAC_AC_Hom", $info, "0", 0);
	$exac_hom_nfe = extract_numeric("EXAC_Hom_NFE", $info, "0", 0);
	$exac_hom_afr = extract_numeric("EXAC_Hom_AFR", $info, "0", 0);
	$exac_hom = "$exac_hom_all,$exac_hom_nfe,$exac_hom_afr";
	$exac_afr = extract_numeric("EXAC_AF_AFR", $info, "0", 4);
	$exac_amr = extract_numeric("EXAC_AF_AMR", $info, "0", 4);
	$exac_eas = extract_numeric("EXAC_AF_EAS", $info, "0", 4);
	$exac_nfe = extract_numeric("EXAC_AF_NFE", $info, "0", 4);
	$exac_sas = extract_numeric("EXAC_AF_SAS", $info, "0", 4);
	$exac_sub = "$exac_afr,$exac_amr,$exac_eas,$exac_nfe,$exac_sas";
	
	//effect predicions
	$phylop = extract_numeric("dbNSFP_phyloP100way_vertebrate", $info, "", 4, "max");
	$sift = extract_string("dbNSFP_SIFT_pred", $info, "");
	$metalr = extract_string("dbNSFP_MetaLR_pred", $info, "");
	$pp2 = extract_string("dbNSFP_Polyphen2_HVAR_pred", $info, "");
	$fathmm = extract_string("dbNSFP_FATHMM_pred", $info, "");
	$cadd = extract_numeric("dbNSFP_CADD_phred", $info, "", 2, "max");
	if ($cadd!="") $cadd = number_format($cadd, 2);
	$pp_count = contains($metalr, "D") + contains($sift, "D") + contains($pp2, "D") + contains($fathmm, "D") + ($cadd!="" && $cadd>20) + ($phylop>=1.6);
	if ($pp_count>0)
	{
		$filter[] = "pred_pathogenic";
	}
	if ($pp_count>2)
	{
		$filter[] = "pred_pathogenic_3";
	}
	
	//OMIM
	$omim = strtr(extract_string("OMIM", $info, ""), "_", " ");
	if ($omim!="")
	{
		$omim .= ";";	
		$filter[] = "anno_omim";
	}
	
	//ClinVar
	$clin_acc = explode("|", extract_string("CLINVAR_ACC", $info, ""));
	$clin_sig = explode("|", extract_string("CLINVAR_SIG", $info, ""));
	$clin_dis = explode("|", extract_string("CLINVAR_DISEASE", $info, ""));
	if (count($clin_acc)!=count($clin_sig)) trigger_error("ClinVar field counts do not match: ACC:".count($clin_acc)." SIG:".count($clin_sig)." DISEASE:".count($clin_dis) , E_USER_ERROR);
	while (count($clin_dis)<count($clin_acc))
	{
		$clin_dis[] = "";
	}
	$clinvar = "";
	for($i=0; $i<count($clin_acc); ++$i)
	{
		if (trim($clin_acc[$i]=="")) continue;
		$disease = trim($clin_dis[$i]);
		if ($disease!="") $disease = " DISEASE=".$disease;
		$clinvar .= $clin_acc[$i]." [".strtr($clin_sig[$i], "_", " ").$disease."]; ";
	}
	if (contains($clinvar, "pathogenic") && !contains($clinvar, "conflicting")) $filter[] = "anno_pathogenic_clinvar";  //matches "pathogenic" and "likely pathogenic"
	
	//HGMD
	$hgmd_id = explode("|", extract_string("HGMD_ID", $info, ""));
	$hgmd_class = explode("|", extract_string("HGMD_CLASS", $info, ""));
	$hgmd_mut = explode("|", extract_string("HGMD_MUT", $info, ""));
	$hgmd_gene = explode("|", extract_string("HGMD_GENE", $info, ""));
	$hgmd_phen = explode("|", extract_string("HGMD_PHEN", $info, ""));
	if (count($hgmd_id)!=count($hgmd_class) || count($hgmd_id)!=count($hgmd_mut) || count($hgmd_id)!=count($hgmd_phen) || count($hgmd_id)!=count($hgmd_gene)) trigger_error("HGMD field counts do not match:\n".implode("|",$hgmd_id)."\n".implode("|",$hgmd_class)."\n".implode("|",$hgmd_mut)."\n".implode("|",$hgmd_gene)."\n".implode("|",$hgmd_phen)."" , E_USER_ERROR);
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
	fwrite($handle_out, "$chr\t$start\t$end\t$ref\t$alt\t$genotype\t".implode(";", $filter)."\t".implode(";", $quality)."\t".implode(",", $genes)."\t$variant_details\t$coding_and_splicing_details\t$repeatmasker\t$dbsnp\t$kgenomes\t$exac\t$exac_hom\t$exac_sub\t$gnomad\t$phylop\t$sift\t$metalr\t$pp2\t$fathmm\t$cadd\t$omim\t$clinvar\t$hgmd\t$cosmic\n");
}

//if no variants are present, we need to write the header line after the loop
if ($in_header) 
{
	write_header_line($handle_out, $column_desc, $filter_desc);
}

fclose($handle);
fclose($handle_out);

?>
