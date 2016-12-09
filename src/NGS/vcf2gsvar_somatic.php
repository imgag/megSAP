<?php 
/** 
	@page vcf2gsvar
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("vcf2gsvar_somatic", "\$Rev: 874 $", "Converts an annotated VCF file to GSvar file.");
$parser->addInfile("in", "Input file in VCF format.", false);
$parser->addOutfile("out", "Output file in VCF format.", false);
$parser->addString("t_col", "Column name of tumor sample.",false);
$parser->addString("n_col", "Column name of tumor sample.",true,"na");
//optional
extract($parser->parse($argv));

//extract values from vcf-columns
function extract_from_info_field($field_name, $info_column, $default = "", $error = TRUE)
{
	//do a type check based on the info field information
	
	$delimiter = ";";
	$fs = explode($delimiter, $info_column);

	foreach($fs as $f)
	{
//		print $f."\n";
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
//$anno_cols[] = array("map_q", "Quality parameters - MAP quality (QUAL).");
$anno_cols[] = array("gene", "Affected gene list (comma-separated).");
$anno_cols[] = array("variant_type", "Variant type.");
$anno_cols[] = array("coding_and_splicing", "Coding and splicing details (Gene, NM number, type, impact, exon number, HGVS.c, HGVS.p).");
$anno_cols[] = array("interpro", "Interpro domains.");
$anno_cols[] = array("RepeatMasker", "RepeatMasker annotation.");
$anno_cols[] = array("dbSNP", "Identifier in dbSNP database.");
$anno_cols[] = array("1000g", "Allele frequency in all populations of 1000g project.");
$anno_cols[] = array("ExAC", "Allele frequency in all populations of ExAC project.");
$anno_cols[] = array("Kaviar", "Allele frequency in Kaviar database.");
$anno_cols[] = array("COSMIC", "COSMIC somatic variant database anntotation.");
$anno_cols[] = array("OMIM", "OMIM database annotation.");
$anno_cols[] = array("ClinVar", "ClinVar database annotation.");
$anno_cols[] = array("HGMD", "HGMD database annotation.");
$anno_cols[] = array("phyloP", "phyloP (100way vertebrate) annotation. Deleterious threshold > 1.6.");
$anno_cols[] = array("MetaLR", "MetaLR effect prediction: D=damaging, T=tolerated.");
$anno_cols[] = array("Sift", "Sift effect prediction: D=damaging, T=tolerated.");
$anno_cols[] = array("PP2_HVAR", "Polyphen2 (HVAR) effect prediction: D=probably damaging, P=possibly damaging, B=benign.");
$anno_cols[] = array("PP2_HDIV", "Polyphen2 (HDIV) effect prediction: D=probably damaging, P=possibly damaging, B=benign.");

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
while($has_header===FALSE && !feof($handle))
{
	$line = nl_trim(fgets($handle));
	if(!empty($line))	++$count;
	
	if(strpos($line,"##source")===0)	//information about variant caller
	{
		$tmp_vcaller = "";
		if(stripos($line,"freebayes")!==FALSE)	$tmp_vcaller = "freebayes";
		if(stripos($line,"strelka")!==FALSE)	$tmp_vcaller = "strelka";
		if(!is_null($var_caller) && $var_caller!=$tmp_vcaller)	trigger_error("Two different source informations identified ($var_caller != $tmp_vcaller!)",E_USER_ERROR);
		if(!is_null($var_caller) && $var_caller==$tmp_vcaller)	trigger_error("Source information used twice!",E_USER_NOTICE);
		$var_caller = $tmp_vcaller;
	}
	
	if(strpos($line,"##FILTER")===0)	//filter column
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
	$normal_idx = -1;
	if(!empty($normal_idx) && $normal_idx!="na")	$normal_idx = array_search($n_col,$header);
	//if($sample_col==-1)	trigger_error("VCF does not contain a column $col.", E_USER_ERROR);
}

//write header
fwrite($handle_out, implode("\n",$filters)."\n");
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
	
	//allele independent annotations
	//chr,start,end,ref,obs = alt
	if(chr_check($chr, 22, false) === FALSE) continue; //skip bad chromosomes
	if($alt==".")	continue;	//skip missing alternative Alleles
	$start = $cols[1];
	$end = $start;

	//quality
	$snp_q = $cols[5];	//column 5 is normally empty in strelka vcfs
	$variant = "SNV";
	$as = explode(",",$alt);
	if(strlen($ref) > 1)	$variant = "INDEL";	//different QC-values for SNVs and Indels
	foreach($as as $a)	//check for different alleles
	{
		if($variant=="SNV" && strlen($a) > 1)	$variant = "INDEL";	//different QC-values for SNVs and Indels
	}
	if($variant == "SNV" && $var_caller=="strelka")	$snp_q = extract_from_info_field("QSS", $info);
	else if($variant == "INDEL" && $var_caller=="strelka")	$snp_q = extract_from_info_field("QSI", $info);		

	//filter column
	if($filter=="PASS")	$filter = "";	//passed variants were not filtered
	
	//RepeatMasker
	$repeatmasker = extract_from_info_field($info, "REPEATMASKER", "",FALSE);

	//effect predicions
	$phylop = extract_from_info_field("dbNSFP_phyloP100way_vertebrate", $info, "", FALSE);
	$metalr = extract_from_info_field("dbNSFP_MetaLR_pred", $info, "", FALSE);
	$sift = extract_from_info_field("dbNSFP_SIFT_pred", $info, "", FALSE);
	$pp2v = extract_from_info_field("dbNSFP_Polyphen2_HDIV_pred", $info, "", FALSE);
	$pp2d = extract_from_info_field("dbNSFP_Polyphen2_HVAR_pred", $info, "", FALSE);

	$interpro = extract_from_info_field("dbNSFP_Interpro_domain", $info, "",FALSE);
		
	//OMIM
	$omim = strtr(extract_from_info_field("OMIM", $info, "", FALSE), "_", " ");
	if ($omim!="") $omim .= ";";
	
	//ClinVar
	$clin_acc = explode("|", extract_from_info_field("CLINVAR_ACC", $info, "", FALSE));
	$clin_sig = explode("|", extract_from_info_field("CLINVAR_SIG", $info, "", FALSE));
	if (count($clin_acc)!=count($clin_sig)) trigger_error("Clinvar field counts do not match:\n".implode("|",$clin_acc)."\n".implode("|",$clin_sig)."" , E_USER_ERROR);
	$clinvar = "";
	for($i=0; $i<count($clin_acc); ++$i)
	{
		if (trim($clin_acc[$i]=="")) continue;
		$clinvar .= $clin_acc[$i]." [".strtr($clin_sig[$i], "_", " ")."]; ";
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
		$out_line = "";
		
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
				
				if($parts[0]!=$a)	continue;	//different alternative allele
				
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
		$tmp = explode(",",extract_from_info_field("KAVIAR_AF", $info, "0.0000", FALSE));
		$kaviar = $tmp[0];
		if(count($tmp)>1) $kaviar = $tmp[$i];
		
		if($variant == "INDEL")	list($start, $end, $ref, $a) = correct_indel($start, $ref, $a);	//correct indels
		$tumor_depth = extract_from_genotype_field($cols[8], $cols[$tumor_idx], "DP");
		$normal_depth = "na";
		if(!$tumor_only)	$normal_depth = extract_from_genotype_field($cols[8], $cols[$normal_idx], "DP");
		$tumor_af = "na";
		$normal_af = "na";
		if($var_caller=="strelka")	//strelka
		{
			if($variant == "SNV")
			{
				//tumor
				list($nuc_t,) = explode(",", extract_from_genotype_field($cols[8], $cols[$tumor_idx], "TU"));
				list($nuc_a,) = explode(",", extract_from_genotype_field($cols[8], $cols[$tumor_idx], "AU"));
				list($nuc_c,) = explode(",", extract_from_genotype_field($cols[8], $cols[$tumor_idx], "CU"));
				list($nuc_g,) = explode(",", extract_from_genotype_field($cols[8], $cols[$tumor_idx], "GU"));
				$o = 0;
				if($a == "T") $o = $nuc_t;
				if($a == "A") $o = $nuc_a;
				if($a == "C") $o = $nuc_c;
				if($a == "G") $o = $nuc_g;
				$tumor_af = number_format($o/($nuc_a+$nuc_t+$nuc_c+$nuc_g),4);

				//normal
				if(!$tumor_only)
				{
					list($nuc_t,) = explode(",", extract_from_genotype_field($cols[8], $cols[$normal_idx], "TU"));
					list($nuc_a,) = explode(",", extract_from_genotype_field($cols[8], $cols[$normal_idx], "AU"));
					list($nuc_c,) = explode(",", extract_from_genotype_field($cols[8], $cols[$normal_idx], "CU"));
					list($nuc_g,) = explode(",", extract_from_genotype_field($cols[8], $cols[$normal_idx], "GU"));
					$o = 0;
					if($a == "T") $o = $nuc_t;
					if($a == "A") $o = $nuc_a;
					if($a == "C") $o = $nuc_c;
					if($a == "G") $o = $nuc_g;
					if(($nuc_a+$nuc_t+$nuc_c+$nuc_g)!=0)	$normal_af = number_format($o/($nuc_a+$nuc_t+$nuc_c+$nuc_g),4);
				}
			}
			else if($variant == "INDEL")
			{
				list($tir,) = explode(",", extract_from_genotype_field($cols[8], $cols[$tumor_idx], "TIR"));
				list($tar,) = explode(",", extract_from_genotype_field($cols[8], $cols[$tumor_idx], "TAR"));
				if(($tir+$tar)!=0)	$tumor_af = number_format($tir/($tir+$tar),4);
				
				//normal
				if(!$tumor_only)
				{
					list($tir,) = explode(",", extract_from_genotype_field($cols[8], $cols[$normal_idx], "TIR"));
					list($tar,) = explode(",", extract_from_genotype_field($cols[8], $cols[$normal_idx], "TAR"));
					if(($tir+$tar)!=0)	$normal_af = number_format($tir/($tir+$tar),4);
				}
			}
		}
		else if($var_caller=="freebayes")	//freebayes
		{
			$tumor_af = number_format(extract_from_genotype_field($cols[8], $cols[$tumor_idx], "AO")/extract_from_genotype_field($cols[8], $cols[$tumor_idx], "DP"), 4);
			if(!$tumor_only)
			{
				$normal_af = "n/a";
				$d = extract_from_genotype_field($cols[8], $cols[$normal_idx], "DP");
				if($d > 0)	$normal_af = number_format(extract_from_genotype_field($cols[8], $cols[$normal_idx], "AO")/$d, 4);
			}
		}
		else
		{
			trigger_error("Could not identify variant caller! Source information of vcf file empty?",E_USER_ERROR);
		}
		
		$out = "$chr\t$start\t$end\t$ref\t$a\t$tumor_af\t$tumor_depth\t";
		if(!$tumor_only)	$out .= "$normal_af\t$normal_depth\t";
		$out .= "$filter\t$snp_q\t$genes\t$variant_details\t$coding_and_splicing\t$interpro\t";
		$out .= "$repeatmasker\t$dbsnp\t$kgenomes\t$exac\t$kaviar\t$cosmic\t$omim\t$clinvar\t$hgmd\t$phylop\t$metalr\t$sift\t$pp2v\t$pp2d\n";
		$out_lines[] = $out;
	}
	
	// write line to out
	fwrite($handle_out, implode("",$out_lines));
}

gzclose($handle);
fclose($handle_out);

?>
