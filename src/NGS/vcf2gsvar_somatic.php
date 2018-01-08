<?php 
/** 
	@page vcf2gsvar
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("vcf2gsvar_somatic", "Converts an annotated VCF file to GSvar file.");
$parser->addInfile("in", "Input file in VCF format.", false);
$parser->addOutfile("out", "Output file in VCF format.", false);
$parser->addString("t_col", "Column name of tumor sample.",false);
$parser->addString("n_col", "Column name of tumor sample.",true,"na");
//optional
$parser->addFlag("strand","Field indicating read distribution of sequenced strand (INFO column field 'fs') available.");
extract($parser->parse($argv));

//extract values from vcf-columns
function extract_from_info_field($field_name, $info_column, $default = "", $error = TRUE)
{
	//do a type check based on the info field information
	
	$delimiter = ";";
	$fs = explode($delimiter, $info_column);

	foreach($fs as $f)
	{
		//return field values
		if(!contains($f, "=")) continue;	//strelka has fields without values
		list($f_name, $f_values) = explode("=", $f);
		if($f_name==$field_name)
		{
			$tmp_values = array();
			foreach(explode(",",$f_values) as $f_value)
			{
				if($f_value==".")	$f_value = $default;
				$tmp_values[] =$f_value;
			}
			return implode(",",$tmp_values);
		}
	}

	if(!$error)	return $default;
	trigger_error("Could not find field '$field_name' in info column.", E_USER_ERROR);
}

//extract from genotype field (sample column)
function extract_from_genotype_field($genotype_column, $genotype_values, $field_name, $delimiter = ":")
{
	$fs = explode($delimiter, $genotype_column);
	$vs = explode($delimiter, $genotype_values);
	for($j=0; $j<count($fs);++$j)
	{
		if($fs[$j] == $field_name)
		{
			return $vs[$j];
		}
	}
	trigger_error("Could not find field '$field_name' in sample column.", E_USER_ERROR);
}

//get sample_ids and column indices
//suggested columns:
//multiple af,qualities (SNP,MAP),filter,gene,variant_type,conding_and_splicing,dbSNP,1000g,ExAC,ESP,phylop,Metlr,sift,pp2_HVAR,pp2_HDIV,omim,clinvar,hgmd,cosmic
//annotated later: inhouse somatic/germline

$tumor_only = false;
if(empty($n_col) || $n_col == "na")	$tumor_only = true;

//init
$anno_cols = array();
$anno_cols[] = array("tumor_af", "Mutant allele frequency in tumor (Sample ".$t_col.").");
$anno_cols[] = array("tumor_dp", "Tumor Depth (Sample ".$t_col.").");
if(!$tumor_only)	$anno_cols[] = array("normal_af", "Mutant allele frequency in normal (Sample ".$n_col.").");
if(!$tumor_only)	$anno_cols[] = array("normal_dp", "Normal depth (Sample ".$n_col.").");
$anno_cols[] = array("filter", "Filter criteria from vcf file.");
$anno_cols[] = array("snp_q", "Quality parameters - SNP quality (QUAL).");
if($strand)	$anno_cols[] = array("strand_tumor","Strand information. Format: [mutation_plus]|[mutation_minus]|[mutation_unkown],[wildtype_plus]|[wildtype_minus]|[wildtype_unkown],[amplicons_plus]|[amplicons_minus].", $t_col."_fs");
if($strand && !$tumor_only)	$anno_cols[] = array("strand_normal","Strand information. Format: [mutation_plus]|[mutation_minus]|[mutation_unkown],[wildtype_plus]|[wildtype_minus]|[wildtype_unkown],[amplicons_plus]|[amplicons_minus].", $n_col."_fs");
//$anno_cols[] = array("map_q", "Quality parameters - MAP quality (QUAL).");
$anno_cols[] = array("gene", "Affected gene list (comma-separated).", "genes");
$anno_cols[] = array("variant_type", "Variant type.", "variant_details");
$anno_cols[] = array("coding_and_splicing", "Coding and splicing details (Gene, ENST number, type, impact, exon number, HGVS.c, HGVS.p).");
$anno_cols[] = array("interpro", "Interpro domains.");
$anno_cols[] = array("RepeatMasker", "RepeatMasker annotation.", "repeatmasker");
$anno_cols[] = array("dbSNP", "Identifier in dbSNP database.", "dbsnp");
$anno_cols[] = array("1000g", "Allele frequency in all populations of 1000g project.","kgenomes");
$anno_cols[] = array("ExAC", "Allele frequency in all populations of ExAC project.", "exac");
$anno_cols[] = array("gnomAD", "Allele frequency in gnomAD database.", "gnomad");
$anno_cols[] = array("COSMIC", "COSMIC somatic variant database anntotation.", "cosmic");
$anno_cols[] = array("OMIM", "OMIM database annotation.", "omim");
$anno_cols[] = array("ClinVar", "ClinVar database annotation.", "clinvar");
$anno_cols[] = array("HGMD", "HGMD database annotation.", "hgmd");
$anno_cols[] = array("phyloP", "phyloP (100way vertebrate) annotation. Deleterious threshold > 1.6.", "phylop");
$anno_cols[] = array("MetaLR", "MetaLR effect prediction: D=damaging, T=tolerated.", "metalr");
$anno_cols[] = array("Sift", "Sift effect prediction: D=damaging, T=tolerated.", "sift");
$anno_cols[] = array("PolyPhen2", "Polyphen2 (HVAR) effect prediction: D=probably damaging, P=possibly damaging, B=benign.", "pp2");
$anno_cols[] = array("FATHMM", "FATHMM effect prediction: D=damaging, T=tolerated.", "fathmm");
$anno_cols[] = array("CADD", "CADD pathogenicity prediction scores (scaled phred-like). Deleterious threshold > 15-20.", "cadd");
$anno_cols[] = array("ihdb_hom", "Homozygous variant counts in NGSD for the same processing system.");
$anno_cols[] = array("ihdb_het", "Heterozyous variant counts in NGSD for the same processing system.");
$anno_cols[] = array("ihdb_wt", "Wildtype variant counts in NGSD for the same processing system.");
$anno_cols[] = array("ihdb_allsys_hom", "Homozygous variant counts in NGSD independent of the processing system.");
$anno_cols[] = array("ihdb_allsys_het", "Heterozygous variant counts in NGSD independent of the processing system.");
$anno_cols[] = array("classification", "Classification from the NGSD.");
$anno_cols[] = array("classification_comment", "Classification comment from the NGSD.");
$anno_cols[] = array("validated", "Validation information from the NGSD. Validation results of other samples are listed in brackets!");
$anno_cols[] = array("comment", "Comments from the NGSD. Comments of other samples are listed in brackets!");
$anno_cols[] = array("gene_info", "General gene information from the NGSD.");
$anno_cols[] = array("som_ihdb_c", "Somatic variant count within NGSD.");
$anno_cols[] = array("som_ihdb_p", "Projects with somatic variant in NGSD.");

//write column descriptions
$handle_out = fopen($out, "w");
foreach($anno_cols as $entry)
{
	fwrite($handle_out, "##DESCRIPTION=".$entry[0]."=".$entry[1]."\n");
}

//read vcf header, extract variant caller, add filter information
$handle = gzopen($in, "r");
$line = "";
$count  =0;
$header = array();
$has_header = FALSE;
$var_caller = NULL;
$filters = array();
$samples = array();
while($has_header===FALSE && !feof($handle))
{
	$line = nl_trim(fgets($handle));
	if(!empty($line))	++$count;
	
	if(strpos($line,"##source=")===0)	//information about variant caller
	{
		$tmp_vcaller = "";
		
		if(stripos($line,"freebayes")!==FALSE)	$tmp_vcaller = "freebayes";
		else if(stripos($line,"strelka")!==FALSE)	$tmp_vcaller = "strelka";
		else if(stripos($line,"Torrent Variant Caller")!==FALSE)	$tmp_vcaller = "iontorrent";
		else	trigger_error("Could not identify variant caller '$line'!",E_USER_ERROR);

		if(!is_null($var_caller) && $var_caller!=$tmp_vcaller)	trigger_error("Two different source informations identified ($var_caller != $tmp_vcaller!)",E_USER_ERROR);
		if(!is_null($var_caller) && $var_caller==$tmp_vcaller)	trigger_error("Source information used twice!",E_USER_NOTICE);
		
		$var_caller = $tmp_vcaller;
		
		if(empty($var_caller))	trigger_error("Could not identify variant caller. Source information available in vcf file?",E_USER_ERROR);
	}
	
	if(strpos($line,"##FILTER")===0)	//filter columns
	{
		$id = NULL;
		$desc = NULL;
		list(,$m) = explode("=",$line,2);
		$m = trim($m,'><');
		$fs = explode(",",$m);
		foreach($fs as $f)
		{
			list($n,$v) = explode("=",$f);
			if($n=="ID")	$id = $v;
			if($n=="Description")	$desc = trim($v,"\"");
		}
		
		if(is_null($id) || is_null($desc))	trigger_error("Could not identify Filter in line $line",E_USER_ERROR);		
		$filters[] = "##FILTER=$id=$desc";
	}
	
	if(starts_with($line,"##SAMPLE="))	//sample columns
	{
		$samples[] = $line;
	}
	
	if(strpos($line,"#CHROM")!==FALSE)	//final header line
	{
		$header = explode("\t",$line);
		$has_header=TRUE;
	}
}
if($count==0)	trigger_error("No variants found.",E_USER_WARNING);
if($count>0)
{
	if($has_header===FALSE)	trigger_error("VCF does not contain a header line.", E_USER_ERROR);
	$tumor_idx = array_search($t_col,$header);
	if($tumor_idx===FALSE)	trigger_error("Could not identify tumor column '$t_col'.",E_USER_ERROR);
	$normal_idx = null;
	if(!$tumor_only)
	{
		$normal_idx = array_search($n_col,$header);
		if($normal_idx===FALSE)	trigger_error("Could not identify normal column '$n_col'.",E_USER_ERROR);
	}
}

//write header
fwrite($handle_out, implode("\n",$filters)."\n");
fwrite($handle_out, implode("\n",$samples)."\n");
fwrite($handle_out, "##ANALYSISTYPE=SOMATIC_".($tumor_only ? "SINGLESAMPLE" : "PAIR")."\n");
fwrite($handle_out, "#chr\tstart\tend\tref\tobs");
foreach($anno_cols as $entry)
{
	fwrite($handle_out, "\t".$entry[0]);
}
fwrite($handle_out, "\n");

//
while(!feof($handle))
{
	$line = nl_trim(fgets($handle));
	++$count;
	if ($line=="" || $line[0]=="#") continue;	
	
	$cols = explode("\t", $line);
	if (count($cols)<9) trigger_error("VCF file line contains less than 10 columns: $line\n", E_USER_ERROR);
	list($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format) = $cols;
	
	// use gt field to identify reference and alternative allele (if available); mainly for use with iontorrent variant-calling
	$idx_al2 = null;
	if(in_array("GT",explode(":",$cols[8])))
	{
		$gt = extract_from_genotype_field($cols[8], $cols[$tumor_idx], "GT");
		
		if($gt=="0/0" || $gt=="0|0" || $gt=="./.")	continue;	// continue if only reference or no allele present
		
		$sep = "/";
		if(strpos($gt,"|")!==FALSE)	$sep = "|";
		
		if(count(explode($sep,$gt))>2)	trigger_error("Sample not diploid.", E_USER_ERROR);
		$idx_al1 = explode($sep,$gt)[0];
		$idx_al2 = explode($sep,$gt)[1];
		$al = array_merge(explode(",",$ref), explode(",",$alt));
		
		if($al[$idx_al1]!=0 && $al[$idx_al1]!=$al[$idx_al2])	trigger_error("Multiallelic variant.",E_USER_ERROR);
		$alt = $al[$idx_al2];
	}
	
	if(chr_check($chr, 22, false) === FALSE) continue; //skip bad chromosomes
	if($alt==".")	continue;	//skip missing alternative Alleles

	//allele independent annotations
	//filter column
	if($filter=="PASS")	$filter = "";	//passed variants were not filtered
	
	//RepeatMasker
	$repeatmasker = extract_from_info_field($info, "REPEATMASKER", "",FALSE);

	//effect predictions
	$phylop = extract_from_info_field("dbNSFP_phyloP100way_vertebrate", $info, "", FALSE);
	$metalr = extract_from_info_field("dbNSFP_MetaLR_pred", $info, "", FALSE);
	$sift = extract_from_info_field("dbNSFP_SIFT_pred", $info, "", FALSE);
	$pp2 = extract_from_info_field("dbNSFP_Polyphen2_HVAR_pred", $info, "", FALSE);
	$fathmm = extract_from_info_field("dbNSFP_FATHMM_pred", $info, "", FALSE);
	$cadd = max(explode(",",extract_from_info_field("dbNSFP_CADD_phred", $info, "", FALSE)));
	if ($cadd!="") $cadd = number_format($cadd, 2);

	$interpro = extract_from_info_field("dbNSFP_Interpro_domain", $info, "",FALSE);
		
	//OMIM
	$omim = strtr(extract_from_info_field("OMIM", $info, "", FALSE), "_", " ");
	if ($omim!="") $omim .= ";";
	
	//ClinVar
	$clin_acc = explode("|", extract_from_info_field("CLINVAR_ACC", $info, "", FALSE));
	$clin_sig = explode("|", extract_from_info_field("CLINVAR_SIG", $info, "", FALSE));
	$clin_dis = explode("|", extract_from_info_field("CLINVAR_DISEASE", $info, "", FALSE));
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
	
	//HGMD
	$hgmd_id = explode("|", extract_from_info_field("HGMD_ID", $info, "",FALSE));
	$hgmd_class = explode("|", extract_from_info_field("HGMD_CLASS", $info, "",FALSE));
	$hgmd_mut = explode("|", extract_from_info_field("HGMD_MUT", $info, "",FALSE));
	$hgmd_gene = explode("|", extract_from_info_field("HGMD_GENE", $info, "",FALSE));
	$hgmd_phen = explode("|", extract_from_info_field("HGMD_PHEN", $info, "",FALSE));
	if (count($hgmd_id)!=count($hgmd_class) || count($hgmd_id)!=count($hgmd_mut) || count($hgmd_id)!=count($hgmd_phen) || count($hgmd_id)!=count($hgmd_gene)) trigger_error("HGMD filed counts do not match:\n".implode("|",$hgmd_id)."\n".implode("|",$hgmd_class)."\n".implode("|",$hgmd_mut)."\n".implode("|",$hgmd_gene)."\n".implode("|",$hgmd_phen)."" , E_USER_ERROR);
	$hgmd = "";
	for($i=0; $i<count($hgmd_id); ++$i)
	{
		if (trim($hgmd_id[$i]=="")) continue;
		$hgmd .= $hgmd_id[$i]." [CLASS=".$hgmd_class[$i]." MUT=".$hgmd_mut[$i]." PHEN=".strtr($hgmd_phen[$i], "_", " ")." GENE=".$hgmd_gene[$i]."]; ";
	}
	
	//COSMIC
	$cosmic = strtr(extract_from_info_field("COSMIC_ID", $info, "",FALSE), array(","=>", "));

	//NGSD
	$ihdb_hom = extract_from_info_field("ihdb_hom", $info, "", FALSE);
	$ihdb_het = extract_from_info_field("ihdb_het", $info, "", FALSE);
	$ihdb_wt = extract_from_info_field("ihdb_wt", $info, "", FALSE);
	$ihdb_allsys_hom = extract_from_info_field("ihdb_allsys_hom", $info, "", FALSE);
	$ihdb_allsys_het = extract_from_info_field("ihdb_allsys_het", $info, "", FALSE);
	$classification = extract_from_info_field("classification", $info, "", FALSE);
	$classification_comment = extract_from_info_field("classification_comment", $info, "", FALSE);
	$validated = extract_from_info_field("validated", $info, "", FALSE);
	$comment = extract_from_info_field("comment", $info, "", FALSE);
	$gene_info = extract_from_info_field("gene_info", $info, "", FALSE);
	$som_ihdb_c = extract_from_info_field("som_ihdb_c", $info, "", FALSE);
	$som_ihdb_p = extract_from_info_field("som_ihdb_p", $info, "", FALSE);

	//allele dependent annotations
	//variant details
	//ANN field: Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO
	//allele dependent annotation
	//variant type, allele frequencies and depth
	$out_lines = array();
	$as = explode(",",$alt);
	for($i=0;$i<count($as);++$i)
	{
		$a = $as[$i];
		
		//variant type
		$variant = "SNV";
		if(strlen($ref) > 1 || strlen($a) > 1)	$variant = "INDEL";	//different QC-values for SNVs and Indels
	
		//chr,start,end,ref,obs = alt
		$start = $pos;
		$end = $start;
		if($variant == "INDEL" && $ref!="-" && $a!="-")	list($start, $end, $ref, $a) = correct_indel($start, $ref, $a);	//correct indels

		//quality
		$snp_q = $cols[5];	//column 5 is normally empty in strelka vcfs
		if($variant == "SNV" && $var_caller=="strelka")	$snp_q = extract_from_info_field("QSS", $info);
		else if($variant == "INDEL" && $var_caller=="strelka")	$snp_q = extract_from_info_field("QSI", $info);		
		
		//strand
		if($strand)
		{
			$strand_tumor = extract_from_info_field((!empty($t_col)?$t_col."_fs":"fs"), $info);
			if(!$tumor_only)	$strand_normal = extract_from_info_field((!empty($n_col)?$n_col."_fs":"fs"), $info);
		}
		
		$genes = array();
		$variant_details = array();
		$coding_and_splicing = array();
		$ann = extract_from_info_field("ANN",$info,"", FALSE);
		if ($ann!="")
		{
			$anns = explode(",", $ann);
			foreach($anns as $entry)
			{
				$parts = explode("|", $entry);
				
				if($parts[0]!=$as[$i])	continue;	//different alternative allele, use original vcf allele (esp. for indels)
				
				$details = $parts[1];
				$variant_details[] = $details;
				$genes[] = $parts[3];
				$exon = $parts[8];
				if ($exon!="") $exon = "exon".$exon;
				$coding_and_splicing[] = $parts[3].":".$parts[6].":$details:".$parts[2].":$exon:".$parts[9].":".$parts[10];
			}
		}
		$genes = implode(",", array_unique($genes));
		$variant_details = implode(",", array_unique($variant_details));
		$coding_and_splicing =  implode(",", $coding_and_splicing);

		//Frequency databases
		$dbsnp = "";
		if ($id!=".") $dbsnp = $id;	//allele independent
		$tmp = explode(",",extract_from_info_field("T1000GP_AF", $info, "0.0000", FALSE));
		$kgenomes = $tmp[0];	//no T1000GP_AF available
		if(count($tmp)>1) $kgenomes = $tmp[$i];
		$tmp = explode(",",extract_from_info_field("EXAC_AF", $info, "0.0000", FALSE));
		$exac = $tmp[0];
		if(count($tmp)>1) $exac = $tmp[$i];
		$tmp = explode(",",extract_from_info_field("GNOMAD_AF", $info, "0.0000", FALSE));
		$gnomad = $tmp[0];
		if(count($tmp)>1) $gnomad = $tmp[$i];
		
		// calculate frequencies
		$tumor_dp = extract_from_genotype_field($cols[8], $cols[$tumor_idx], "DP");
		$normal_dp = "na";
		if(!$tumor_only)	$normal_dp = extract_from_genotype_field($cols[8], $cols[$normal_idx], "DP");
		$tumor_af = "na";
		$normal_af = "na";
		if($var_caller=="strelka")	//strelka
		{
			if($variant == "SNV")
			{
				list($tumor_dp,$tumor_af) = vcf_strelka_snv($cols[8],$cols[$tumor_idx],$a);
				if(!$tumor_only)	list($normal_dp,$normal_af) = vcf_strelka_snv($cols[8],$cols[$normal_idx],$a);
			}
			else if($variant == "INDEL")
			{
				list($tumor_dp,$tumor_af) = vcf_strelka_indel($cols[8],$cols[$tumor_idx]);
				if(!$tumor_only)	list($normal_dp,$normal_af) = vcf_strelka_indel($cols[8],$cols[$normal_idx]);
			}
		}
		else if($var_caller=="freebayes")	//freebayes
		{
			list($tumor_dp, $tumor_af) = vcf_freebayes($cols[8], $cols[$tumor_idx]);
			if(!$tumor_only)
			{
				list($normal_dp, $normal_af) = vcf_freebayes($cols[8], $cols[$normal_idx]);
			}
		}
		else if($var_caller=="iontorrent")	//iontorrent
		{
			list($tumor_dp, $tumor_af) = vcf_iontorrent($cols[8], $cols[$tumor_idx],$idx_al2-1);
			if(!$tumor_only)	list($normal_dp, $normal_af) = vcf_freebayes($cols[8], $cols[$idx_al2]-1);
		}

		$out = array($chr,$start,$end,$ref,$a);
		foreach($anno_cols as $col)
		{
			$var = $col[0];
			if(isset($col[2]) && !isset($$var))	$var = $col[2];
			$out[] = $$var;
		}
		$out_lines[] = implode("\t",$out)."\n";
	}
	
	// write line to out
	fwrite($handle_out, implode("",$out_lines));
}

gzclose($handle);
fclose($handle_out);

?>
