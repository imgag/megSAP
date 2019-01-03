<?php 
/** 
	@page vcf2gsvar
	
	@todo replace blacklist genes by blacklist regions:
			- /mnt/share/data/dbs/ABB/ABB_075.bed
			- /mnt/share/data/blacklist_regions/blacklist.bed 
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("vcf2gsvar", "Converts an annotated VCF file from freebayes to a GSvar file.");
$parser->addInfile("in",  "Input file in VCF format.", false);
$parser->addOutfile("out", "Output file in GSvar format.", false);
//optional
$parser->addEnum("genotype_mode", "Genotype handling mode.", true, array("single", "multi", "skip"), "single");
$parser->addFlag("updown", "Don't discard up- or downstream anntations (5000 bases around genes).");
$parser->addFlag("blacklist", "Annotate variants in blacklisted genes with 'gene_blacklist' in filter column.");
$parser->addFlag("wgs", "Enables WGS mode: MODIFIER variants with a AF>2% are skipped to reduce the numer of variants to a manageable size.");
extract($parser->parse($argv));

//skip common MODIFIER variants in WGS mode
function skip_in_wgs_mode($chr, $coding_and_splicing_details, $kg, $gnomad, $clinvar, $hgmd)
{
	//don't skip mito variants
	if ($chr=='chrMT') return false;
	
	//don't skip exonic variants
	if (contains($coding_and_splicing_details, ":LOW:") || contains($coding_and_splicing_details, ":MODERATE:") || contains($coding_and_splicing_details, ":HIGH:")) return false;
	
	//don't skip variants annotated to be (likely) pathognic
	if (contains($hgmd, "CLASS=DM") || (contains($clinvar, "pathogenic") && !contains($clinvar, "conflicting"))) return false;	
	
	//skip common variants >2%AF
	if ($kg!="" && $kg>0.02) return true;
	if ($gnomad!="" && $gnomad>0.02) return true;
	
	return false; //non-exonic but rare
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

//get index of columnn in QSC header.
function index_of($cols, $name, $error_if_missing = true)
{
	$index = array_search($name, $cols);
	if ($index===FALSE && $error_if_missing)
	{
		trigger_error("Could not find column '$name' in VEP QSC annotation. Valid column names are: ".implode(", ", array_values($cols)), E_USER_ERROR);
	}
	return $index;
}

//translate value (and throw error if not valid)
function translate($error_name, $value, $dict)
{
	$value = trim($value);
	
	if (!isset($dict[$value]))
	{
		trigger_error("Cannot translate value '{$value}' for '{$error_name}'. Valid terms are: '".implode("','", array_keys($dict)), E_USER_ERROR);
	}
	
	return $dict[$value];
}

//collapse several values to a single value
function collapse($error_name, $values, $mode, $decimal_places = null)
{
	for($i=0; $i<count($values); ++$i)
	{
		$v = trim($values[$i]);
		if (!is_null($decimal_places) && $v!="")
		{
			if (!is_numeric($v))
			{
				trigger_error("Invalid numeric value '{$v}' in mode '{$mode}' while collapsing '{$error_name}'!", E_USER_ERROR);
			}
			$v = number_format($v, $decimal_places, ".", "");
			if ($values[$i]>0 && $v==0)
			{
				$v = substr($v, 0, -1)."1";
			}
		}
		$values[$i] = $v;
	}

	if ($mode=="one")
	{
		$values = array_unique($values);
		if (count($values)>1)
		{
			trigger_error("Several values '".implode("','", $values)."' in mode '{$mode}' while collapsing '{$error_name}'!", E_USER_ERROR);
		}
		else if (count($values)==0)
		{
			return "";
		}
		return $values[0];
	}
	else if ($mode=="max")
	{
		return max($values);
	}
	else if ($mode=="unique")
	{
		return array_unique($values);
	}
	else
	{
		trigger_error("Invalid mode '{$mode}' while collapsing '{$error_name}'!", E_USER_ERROR);
	}
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
	array("quality", "Quality parameters - variant quality (QUAL), depth (DP), allele frequency (AF), mean mapping quality of alternate allele (MQM)."),
	array("gene", "Affected gene list (comma-separated)."),
	array("variant_type", "Variant type."),
	array("coding_and_splicing", "Coding and splicing details (Gene, ENST number, type, impact, exon/intron number, HGVS.c, HGVS.p, Pfam domain)."),
	array("regulatory", "Regulatory consequence details."),
	array("OMIM", "OMIM database annotation."),
	array("ClinVar", "ClinVar database annotation."),
	array("HGMD", "HGMD database annotation."),
	array("RepeatMasker", "RepeatMasker annotation."),
	array("dbSNP", "Identifier in dbSNP database."),
	array("1000g", "Allele frequency in 1000 genomes project."),
	array("gnomAD", "Allele frequency in gnomAD project."),
	array("gnomAD_hom_hemi", "Homoyzgous counts and hemizygous counts of gnomAD project (genome data)."),
	array("gnomAD_sub", "Sub-population allele frequenciens (AFR,AMR,EAS,NFE,SAS) in gnomAD project."),
	array("ESP_sub", "Sub-population allele frequency (EA,AA) in NHLBI Exome Sequencing project."),
	array("phyloP", "phyloP (100way vertebrate) annotation. Deleterious threshold > 1.6."),
	array("Sift", "Sift effect prediction for each transcript: D=damaging, T=tolerated."),
	array("PolyPhen", "PolyPhen (humVar) effect prediction for each transcript: D=probably damaging, P=possibly damaging, B=benign."),
	array("fathmm-MKL", "fathmm-MKL score (for coding/non-coding regions). Deleterious threshold > 0.5."),
	array("CADD", "CADD pathogenicity prediction scores (scaled phred-like). Deleterious threshold > 10-20."),
	array("REVEL", "REVEL pathogenicity prediction score. Deleterious threshold > 0.5."),
	array("MaxEntScan", "MaxEntScan splicing prediction (reference bases score/alternate bases score)."),
	array("GeneSplicer", "GeneSplicer splicing prediction (state/type/coordinates/confidence/score)."),
	array("dbscSNV", "dbscSNV splicing prediction (ADA/RF score)."),
	array("COSMIC", "COSMIC somatic variant database anntotation."),
);
if ($genotype_mode=="single")
{
	array_unshift($column_desc, array("genotype", "Genotype of variant in sample."));	
}

//write filter descriptions
$filter_desc = array();
if ($blacklist) $filter_desc[] = array("gene_blacklist", "The gene(s) are contained on the blacklist of unreliable genes.");

//parse input
$c_written = 0;
$c_skipped_wgs = 0;
$multi_cols = array();
$in_header = true;
$handle = fopen($in, "r");
$handle_out = fopen($out, "w");
while(!feof($handle))
{
	$line = nl_trim(fgets($handle));
	if ($line=="" || trim($line)=="") continue;
	
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
			if ($genotype_mode=="multi")
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
		
		if (starts_with($line, "##PIPELINE="))
		{
			fwrite($handle_out, trim($line)."\n");
		}
		
		//get annotation indices in CSQ field
		if (starts_with($line, "##INFO=<ID=CSQ,"))
		{
			$cols = explode("|", substr($line, 0, -2));
			$cols[0] = "Allele";
			$i_consequence = index_of($cols, "Consequence");
			$i_impact = index_of($cols, "IMPACT");
			$i_symbol = index_of($cols, "SYMBOL");
			$i_feature = index_of($cols, "Feature");
			$i_featuretype = index_of($cols, "Feature_type");
			$i_biotype = index_of($cols, "BIOTYPE");
			$i_exon = index_of($cols, "EXON");
			$i_intron = index_of($cols, "INTRON");
			$i_hgvsc = index_of($cols, "HGVSc");
			$i_hgvsp = index_of($cols, "HGVSp");
			$i_domains = index_of($cols, "DOMAINS");
			$i_sift = index_of($cols, "SIFT");
			$i_polyphen = index_of($cols, "PolyPhen");
			$i_phylop = index_of($cols, "PHYLOP", false);
			$i_cadd = index_of($cols, "CADD_PHRED", false);
			$i_revel = index_of($cols, "REVEL", false);
			$i_fathmm_c = index_of($cols, "FATHMM_MKL_C");
			$i_fathmm_nc = index_of($cols, "FATHMM_MKL_NC");
			$i_polyphen = index_of($cols, "PolyPhen");
			$i_existingvariation = index_of($cols, "Existing_variation");
			$i_af_kg = index_of($cols, "AF");
			$i_af_gnomad = index_of($cols, "gnomAD_AF");
			$i_af_gnomad_genome = index_of($cols, "gnomADg_AF");
			$i_hom_gnomad_genome = index_of($cols, "gnomADg_Hom");
			$i_hemi_gnomad_genome = index_of($cols, "gnomADg_Hemi");
			$i_af_gnomad_afr = index_of($cols, "gnomAD_AFR_AF");
			$i_af_gnomad_amr = index_of($cols, "gnomAD_AMR_AF");
			$i_af_gnomad_eas = index_of($cols, "gnomAD_EAS_AF");
			$i_af_gnomad_nfe = index_of($cols, "gnomAD_NFE_AF");
			$i_af_gnomad_sas = index_of($cols, "gnomAD_SAS_AF");
			$i_af_esp_ea = index_of($cols, "EA_AF");
			$i_af_esp_aa = index_of($cols, "AA_AF");
			$i_repeat = index_of($cols, "REPEATMASKER");
			$i_clinvar = index_of($cols, "CLINVAR");
			$i_clinvar_details = index_of($cols, "CLINVAR_DETAILS");
			$i_omim = index_of($cols, "OMIM", false);
			$i_hgmd = index_of($cols, "HGMD", false);
			$i_hgmd_class = index_of($cols, "HGMD_CLASS", false);
			$i_hgmd_mut = index_of($cols, "HGMD_MUT", false);
			$i_hgmd_gene = index_of($cols, "HGMD_GENE", false);
			$i_hgmd_phen = index_of($cols, "HGMD_PHEN", false);
			$i_maxes_ref = index_of($cols, "MaxEntScan_ref");
			$i_maxes_alt = index_of($cols, "MaxEntScan_alt");
			$i_genesplicer = index_of($cols, "GeneSplicer");
			$i_dbscsnv_ada = index_of($cols, "ada_score");
			$i_dbscsnv_rf = index_of($cols, "rf_score");
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
	if (count($cols)<10) trigger_error("VCF file line contains less than 10 columns: '$line'", E_USER_ERROR);
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
	if ($genotype_mode=="multi")
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
		$genotype = "\t".implode("\t", $genotype);
		$sample["DP"] = implode(",", $depth);
		$sample["AO"] = implode(",", $ao);
	}
	else if ($genotype_mode=="single")
	{
		if (!isset($sample["GT"])) 
		{
			trigger_error("VCF sample column does not contain GT value!", E_USER_ERROR);
		}
		$genotype = vcfgeno2human($sample["GT"]);
		
		//skip wildtype
		if ($genotype=="wt") continue;
		
		$genotype = "\t".$genotype;
	}
	else if ($genotype_mode=="skip")
	{
		$genotype = "";
	}
	else
	{
		trigger_error("Invalid mode '{$genotype_mode}'!", E_USER_ERROR);
	}

	//quality
	$quality = array();
	$qual = intval($qual);
	$quality[] = "QUAL=".$qual;
	if (isset($sample["DP"]))
	{
		$quality[] = "DP=".$sample["DP"];
	}
	if (isset($sample["AO"]) && isset($sample["DP"]))
	{
		//comma-separated values in case of multi-sample data
		$afs = array();
		$aos = explode(",", $sample["AO"]);
		$dps = explode(",", $sample["DP"]);
		for($i=0; $i<count($dps); ++$i)
		{
			if (is_numeric($aos[$i]) && is_numeric($dps[$i]) && $dps[$i]>0)
			{
				$afs[] = number_format($aos[$i]/$dps[$i], 2);
			}
			else if ($genotype_mode=="multi")
			{
				$afs[] = "";
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
	}
	
	//for dragen VCFs
	if (isset($info["MQ"])) 
	{
		$quality[] = "MQM=".intval($info["MQ"]);
	}
	if (isset($sample["AF"]))
	{
		$quality[] = "AF=".$sample["AF"];
	}

	//variant details
	$sift = array();
	$polyphen = array();
	$phylop = array();
	$fathmm = array();
	$cadd = array();
	$revel = array();
	$dbsnp = array();
	$cosmic = array();
	$genes = array();
	$variant_details = array();
	$coding_and_splicing_details = array();
	$af_kg = array();
	$af_gnomad = array();
	$af_gnomad_genome = array();
	$af_gnomad_afr = array();
	$af_gnomad_amr = array();
	$af_gnomad_eas = array();
	$af_gnomad_nfe = array();
	$af_gnomad_sas = array();
	$af_esp_ea = array();
	$af_esp_aa = array();
	$hom_gnomad = array();
	$hemi_gnomad = array();
	$repeat = array();
	$clinvar = array();
	$omim = array();
	$hgmd = array();
	$maxentscan = array();
	$genesplicer = array();
	$dbscsnv = array();
	$regulatory = array();
	
	//variant details (up/down-stream)
	$variant_details_updown = array();
	$genes_updown = array();
	$sift_updown = array();
	$polyphen_updown = array();
	$coding_and_splicing_details_updown = array();
	if (isset($info["CSQ"]))
	{
		$anns = explode(",", $info["CSQ"]);
		foreach($anns as $entry)
		{			
			$parts = explode("|", $entry);
			
			//######################### general information (not transcript-specific) #########################
			
			//AFs
			$af_kg[] = max(explode("&", trim($parts[$i_af_kg]))); //some variants are annotated with two values, e.g. chr7:130297119 GA>G
			$af_gnomad[] = trim($parts[$i_af_gnomad]);
			$gnomad_value = trim($parts[$i_af_gnomad_genome]);
			$af_gnomad_genome[] = $gnomad_value=="." ? "" : $gnomad_value;
			$af_gnomad_afr[] = trim($parts[$i_af_gnomad_afr]);
			$af_gnomad_amr[] = trim($parts[$i_af_gnomad_amr]);
			$af_gnomad_eas[] = trim($parts[$i_af_gnomad_eas]);
			$af_gnomad_nfe[] = trim($parts[$i_af_gnomad_nfe]);
			$af_gnomad_sas[] = trim($parts[$i_af_gnomad_sas]);
			$af_esp_ea[] = trim($parts[$i_af_esp_ea]);
			$af_esp_aa[] = trim($parts[$i_af_esp_aa]);
			$hom_gnomad[] = trim($parts[$i_hom_gnomad_genome]);
			$hemi_gnomad[] = trim($parts[$i_hemi_gnomad_genome]);
			
			//dbSNP, COSMIC
			$ids = explode("&", $parts[$i_existingvariation]);
			foreach($ids as $id)
			{
				if (starts_with($id, "rs"))
				{
					$dbsnp[] = $id;
				}
				if (starts_with($id, "COSM") || starts_with($id, "COSN"))
				{
					$cosmic[] = $id;
				}
			}
			
			//RepeatMasker
			$repeat[] = strtr(trim($parts[$i_repeat]), array("[s]"=>" "));
			
			//ClinVar
			$clin_accs = explode("&", $parts[$i_clinvar]);
			$clin_details = explode("&", $parts[$i_clinvar_details]);
			for ($i=0; $i<count($clin_accs); ++$i)
			{
				if ($clin_accs[$i]!="")
				{
					$clinvar[] = $clin_accs[$i]." [".strtr($clin_details[$i], array("[c]"=>",", "[p]"=>" DISEASE=", "_"=>" "))."];";
				}
			}
			
			//OMIM
			if ($i_omim!==FALSE)
			{
				$text = trim(strtr($parts[$i_omim], array("[c]"=>",", "&"=>",", "_"=>" ")));
				if ($text!="") $text .= ";";
				$omim[] = $text;
			}
			
			//HGMD
			if ($i_hgmd!==FALSE)
			{
				$hgmd_id = explode("|", $parts[$i_hgmd]);
				$hgmd_class = explode("|", $parts[$i_hgmd_class]);
				$hgmd_mut = explode("|", $parts[$i_hgmd_mut]);
				$hgmd_gene = explode("|", $parts[$i_hgmd_gene]);
				$hgmd_phen = explode("|", $parts[$i_hgmd_phen]);
				if (count($hgmd_id)!=count($hgmd_class) || count($hgmd_id)!=count($hgmd_mut) || count($hgmd_id)!=count($hgmd_phen) || count($hgmd_id)!=count($hgmd_gene)) trigger_error("HGMD field counts do not match:\n".implode("|",$hgmd_id)."\n".implode("|",$hgmd_class)."\n".implode("|",$hgmd_mut)."\n".implode("|",$hgmd_gene)."\n".implode("|",$hgmd_phen)."" , E_USER_ERROR);
				$text = "";
				for($i=0; $i<count($hgmd_id); ++$i)
				{
					if (trim($hgmd_id[$i]=="")) continue;
					$text .= $hgmd_id[$i]." [CLASS=".$hgmd_class[$i]." MUT=".$hgmd_mut[$i]." PHEN=".strtr($hgmd_phen[$i], "_", " ")." GENE=".$hgmd_gene[$i]."]; ";
				}
				$hgmd[] = trim($text);
			}
			
			//MaxEntScan
			if ($parts[$i_maxes_ref]!="")
			{
				$result = number_format($parts[$i_maxes_ref], 2).">".number_format($parts[$i_maxes_alt], 2);
				$maxentscan[] = $result;
			}
			
			//GeneSplicer
			if ($parts[$i_genesplicer]!="")
			{
				$genesplicer[] = $parts[$i_genesplicer];
			}
			
			//dbscSNV
			if ($parts[$i_dbscsnv_ada]!="" || $parts[$i_dbscsnv_rf]!="")
			{
				$ada = $parts[$i_dbscsnv_ada]=="" ? "" : number_format($parts[$i_dbscsnv_ada],3);
				$rf = $parts[$i_dbscsnv_rf]=="" ? "" : number_format($parts[$i_dbscsnv_rf],3);
				$dbscsnv[] = "{$ada}/{$rf}";
			}
			
			//pathogenicity predictions (not transcript-specific)
			$phylop_parts = explode("&", trim($parts[$i_phylop]));
			if (count($phylop_parts)>1) //deletions are annotated for each base => use maximum q
			{
				$phylop[] = max($phylop_parts);
			}
			else
			{			
				$phylop[] = trim($parts[$i_phylop]);
			}
			$cadd[] = trim($parts[$i_cadd]);
			$revel_score = trim($parts[$i_revel]);
			if ($revel_score!="") $revel[] = $revel_score;
			$fathmm_c = trim($parts[$i_fathmm_c]);
			$fathmm_nc = trim($parts[$i_fathmm_nc]);
			$fathmm[] = ($fathmm_c=="" && $fathmm_nc=="") ? "" : number_format($fathmm_c, 2).",".number_format($fathmm_nc, 2);
			
			//######################### transcript-specific information #########################
			
			//only transcripts
			$feature_type = $parts[$i_featuretype];
			if ($feature_type=="Transcript" || $feature_type=="")
			{
				$transcript_id = $parts[$i_feature];
				
				//extract variant type
				$variant_type = strtr($parts[$i_consequence], array("_variant"=>""));
				$variant_type = strtr($variant_type, array("splice_acceptor&splice_region&intron"=>"splice_acceptor", "splice_donor&splice_region&intron"=>"splice_donor", "splice_acceptor&intron"=>"splice_acceptor", "splice_donor&intron"=>"splice_donor", "_prime_"=>"'"));
				$is_updown = $variant_type=="upstream_gene" || $variant_type=="downstream_gene";
				if (!$is_updown)
				{
					$variant_details[] = $variant_type;
				}
				else
				{
					$variant_details_updown[] = $variant_type;
				}
				
				//split genes
				$gene = $parts[$i_symbol];
				if (!$is_updown)
				{
					$genes[] = $gene;
				}
				else
				{
					$genes_updown[] = $gene;
				}
				
				//pathogenicity predictions (transcript-specific)
				$sift_entry = translate("Sift", $parts[$i_sift], array(""=>" ", "deleterious"=>"D", "tolerated"=>"T", "tolerated_low_confidence"=>"T", "deleterious_low_confidence"=>"D"));
				if (!$is_updown)
				{
					$sift[] = $sift_entry;
				}
				else
				{
					$sift_updown[] = $sift_entry;
				}
				$polyphen_entry = translate("PolyPhen", $parts[$i_polyphen], array(""=>" ", "unknown"=>" ",  "probably_damaging"=>"D", "possibly_damaging"=>"P", "benign"=>"B"));
				if (!$is_updown)
				{
					$polyphen[] = $polyphen_entry;
				}
				else
				{
					$polyphen_updown[] = $polyphen_entry;
				}
				
				//exon
				$exon = trim($parts[$i_exon]);
				if ($exon!="") $exon = "exon".$exon;
				$intron = trim($parts[$i_intron]);
				if ($intron!="") $intron = "exon".$intron;
				
				//hgvs
				$hgvs_c = trim($parts[$i_hgvsc]);
				if ($hgvs_c!="") $hgvs_c = explode(":", $hgvs_c)[1];			
				$hgvs_p = trim($parts[$i_hgvsp]);
				if ($hgvs_p!="") $hgvs_p = explode(":", $hgvs_p)[1];
				$hgvs_p = str_replace("%3D", "=", $hgvs_p);
				
				//domain
				$domain = "";
				$domains = explode("&", $parts[$i_domains]);
				foreach($domains as $entry)
				{
					if(starts_with($entry, "Pfam_domain:"))
					{
						$domain = explode(":", $entry, 2)[1];
					}
				}
				
				$transcript_entry = "{$gene}:{$transcript_id}:".$parts[$i_consequence].":".$parts[$i_impact].":{$exon}{$intron}:{$hgvs_c}:{$hgvs_p}:{$domain}";
				if (!$is_updown)
				{
					$coding_and_splicing_details[] = $transcript_entry;
				}
				else
				{
					$coding_and_splicing_details_updown[] = $transcript_entry;
				}
			}
			else if ($feature_type=="RegulatoryFeature")
			{
				$regulatory[] = $parts[$i_consequence].":".$parts[$i_biotype];
			}
			else if ($feature_type=="MotifFeature")
			{
				$regulatory[] = $parts[$i_consequence];
			}
			else
			{				
				trigger_error("Unknown VEP feature type '{$feature_type}' for variant $chr:$pos!", E_USER_ERROR);
			}
		}
	}
	
	//add up/down-stream variants if requested (or no other transcripts exist)
	if ($updown || count($coding_and_splicing_details)==0)
	{
		$variant_details = array_merge($variant_details, $variant_details_updown);
		$genes = array_merge($genes, $genes_updown);
		$sift = array_merge($sift, $sift_updown);
		$polyphen = array_merge($polyphen, $polyphen_updown);
		$coding_and_splicing_details = array_merge($coding_and_splicing_details, $coding_and_splicing_details_updown);
	}
	
	$genes = array_unique($genes);
	if ($blacklist && all_genes_blacklisted($genes))
	{
		$filter[] = "gene_blacklist";
	}
	$variant_details = implode(",", array_unique($variant_details));
	$coding_and_splicing_details =  implode(",", $coding_and_splicing_details);
	
	//regulatory
	$regulatory = implode(",", collapse("Regulatory", $regulatory, "unique"));
	
	//RepeatMasker
	$repeatmasker = collapse("RepeatMasker", $repeat, "one");
	
	//AFs
	$dbsnp = implode(",", collapse("dbSNP", $dbsnp, "unique"));	
	$kg = collapse("1000g", $af_kg, "one", 4);
	$gnomad = collapse("gnomAD", $af_gnomad, "one", 4);
	$gnomad_genome = collapse("gnomAD genome", $af_gnomad_genome, "one", 4);
	$gnomad = max($gnomad, $gnomad_genome);
	$gnomad_hom_hemi = collapse("gnomAD Hom", $hom_gnomad, "one").",".collapse("gnomAD Hemi", $hemi_gnomad, "one");
	if ($gnomad_hom_hemi==",") $gnomad_hom_hemi = "";
	$gnomad_sub = collapse("gnomAD AFR", $af_gnomad_afr, "one", 4).",".collapse("gnomAD AMR", $af_gnomad_amr, "one", 4).",".collapse("gnomAD EAS", $af_gnomad_eas, "one", 4).",".collapse("gnomAD NFE", $af_gnomad_nfe, "one", 4).",".collapse("gnomAD SAS", $af_gnomad_sas, "one", 4);
	if (str_replace(",", "", $gnomad_sub)=="") $gnomad_sub = "";
	$esp_sub = collapse("ESP ea", $af_esp_ea, "one", 4).",".collapse("ESP aa", $af_esp_aa, "one", 4);
	if (str_replace(",", "", $esp_sub)=="") $esp_sub = "";
	
	//effect predicions
	$phylop = collapse("phyloP", $phylop, "one", 4);
	$sift = implode(",", $sift);
	if (trim(strtr($sift, ",", " "))=="") $sift = "";
	$polyphen = implode(",", $polyphen);
	if (trim(strtr($polyphen, ",", " "))=="") $polyphen = "";
	$fathmm = collapse("fathmm-MKL", $fathmm, "one");
	$cadd = collapse("CADD", $cadd, "one", 2);
	$revel = empty($revel) ? "" : collapse("REVEL", $revel, "max", 2);
	
	//OMIM
	$omim = collapse("OMIM", $omim, "one");

	//ClinVar
	$clinvar = implode(" ", collapse("ClinVar", $clinvar, "unique"));
	
	//HGMD
	$hgmd = collapse("HGMD", $hgmd, "one");

	//MaxEntScan
	$maxentscan = implode(",", collapse("MaxEntScan", $maxentscan, "unique"));
	
	//GeneSplicer
	$genesplicer = implode(",", collapse("GeneSplicer", $genesplicer, "unique"));
	
	//dbscSNV
	$dbscsnv = empty($dbscsnv) ? "" : collapse("dbscSNV", $dbscsnv, "one");
	
	//COSMIC
	$cosmic = implode(",", collapse("COSMIC", $cosmic, "unique"));
	
	//skip common MODIFIER variants in WGS mode
	if ($wgs && skip_in_wgs_mode($chr, $coding_and_splicing_details, $kg, $gnomad, $clinvar, $hgmd))
	{
		++$c_skipped_wgs;
		continue;
	}
	
	//write data
	++$c_written;
	fwrite($handle_out, "$chr\t$start\t$end\t$ref\t{$alt}{$genotype}\t".implode(";", $filter)."\t".implode(";", $quality)."\t".implode(",", $genes)."\t$variant_details\t$coding_and_splicing_details\t$regulatory\t$omim\t$clinvar\t$hgmd\t$repeatmasker\t$dbsnp\t$kg\t$gnomad\t$gnomad_hom_hemi\t$gnomad_sub\t$esp_sub\t$phylop\t$sift\t$polyphen\t$fathmm\t$cadd\t$revel\t$maxentscan\t$genesplicer\t$dbscsnv\t$cosmic\n");
}

//if no variants are present, we need to write the header line after the loop
if ($in_header) 
{
	write_header_line($handle_out, $column_desc, $filter_desc);
}

fclose($handle);
fclose($handle_out);

//print debug output
print "Variants written: {$c_written}\n";
if ($wgs)
{
	print "Variants skipped because WGS mode is enabled: {$c_skipped_wgs}\n";
}


?>
