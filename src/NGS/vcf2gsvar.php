<?php 
/** 
	@page vcf2gsvar
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("vcf2gsvar", "Converts an annotated VCF file from freebayes to a GSvar file.");
$parser->addInfile("in",  "Input file in VCF format.", false);
$parser->addOutfile("out", "Output file in GSvar format.", false);
//optional
$parser->addEnum("genotype_mode", "Genotype handling mode.", true, array("single", "multi", "skip"), "single");
$parser->addFlag("updown", "Don't discard up- or downstream annotations (5000 bases around genes).");
$parser->addFlag("wgs", "Enables WGS mode: MODIFIER variants with a AF>2% are skipped to reduce the number of variants to a manageable size.");
extract($parser->parse($argv));

//skip common MODIFIER variants in WGS mode
function skip_in_wgs_mode($chr, $coding_and_splicing_details, $kg, $gnomad, $clinvar, $hgmd, $ngsd_clas)
{	
	//don't skip mito variants
	if ($chr=='chrMT') return false;
	
	//don't skip exonic variants
	if (contains($coding_and_splicing_details, ":LOW:") || contains($coding_and_splicing_details, ":MODERATE:") || contains($coding_and_splicing_details, ":HIGH:")) return false;
	
	//don't skip variants annotated to be (likely) pathognic
	if (contains($hgmd, "CLASS=DM") || (contains($clinvar, "pathogenic") && !contains($clinvar, "conflicting"))) return false;	
	
	//don't skip variants of class 4/5 in NGSD
	if ($ngsd_clas=='4' || $ngsd_clas=='5') return false;
	
	//skip common variants >2%AF
	if ($kg!="" && $kg>0.02) return true;
	if ($gnomad!="" && $gnomad>0.02) return true;
	
	return false; //non-exonic but rare
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
function collapse(&$tag, $error_name, $values, $mode, $decimal_places = null)
{
	for($i=0; $i<count($values); ++$i)
	{
		$v = trim($values[$i]);
		if (!is_null($decimal_places) && $v!="")
		{
			if (!is_numeric($v))
			{
				trigger_error("Invalid numeric value '{$v}' in mode '{$mode}' while collapsing '{$error_name}' in variant '{$tag}'!", E_USER_ERROR);
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
			trigger_error("Several values '".implode("','", $values)."' in mode '{$mode}' while collapsing '{$error_name}' in variant '{$tag}'!", E_USER_ERROR);
		}
		else if (count($values)==0)
		{
			return "";
		}
		return $values[0];
	}
	else if ($mode=="max")
	{
		if (count($values)==0) return "";
		return max($values);
	}
	else if ($mode=="unique")
	{
		return array_unique($values);
	}
	else
	{
		trigger_error("Invalid mode '{$mode}' while collapsing '{$error_name}' in variant '{$tag}'!", E_USER_ERROR);
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

//load replaced/removed Pfam IDs
function load_pfam_replacements()
{
	$pfam_filepath = repository_basedir()."/data/misc/pfam_replacements.tsv";
	if (!is_readable($pfam_filepath))
	{
		trigger_error("Pfam replacement file '$pfam_filepath' is not readable!", E_USER_ERROR);
	}
	$pfam_list = file($pfam_filepath);
	$pfam_replacements = array();
	foreach($pfam_list as $line)
	{
		// ignore comments
		if (starts_with($line, '#'))
		{
			continue;
		}

		$split_line = explode("\t",$line);
		if (count($split_line) < 2)
		{
			trigger_error("Error parsing Pfam file '$pfam_filepath'!", E_USER_ERROR);
		}
		$pfam_replacements[trim($split_line[0])] = trim($split_line[1]);
	}
	return $pfam_replacements;
}
$pfam_replacements = load_pfam_replacements();

//load Pfam description
function load_pfam_description()
{
	$pfam_filepath = repository_basedir()."/data/misc/pfam_description.tsv";
	if (!is_readable($pfam_filepath))
	{
		trigger_error("Pfam description file '$pfam_filepath' is not readable!", E_USER_ERROR);
	}
	$pfam_list = file($pfam_filepath);
	$pfam_description = array();
	foreach($pfam_list as $line)
	{
		// ignore comments
		if (starts_with($line, '#'))
		{
			continue;
		}

		$split_line = explode("\t",$line);
		if (count($split_line) < 2)
		{
			trigger_error("Error parsing Pfam file '$pfam_filepath'!", E_USER_ERROR);
		}
		$description_string = trim($split_line[1]);
		$description_string = strtr($description_string, array(":" => " ", "," => "", "[" => "(", "]" => ")"));
		$pfam_description[trim($split_line[0])] = $description_string;
	}
	return $pfam_description;
}
$pfam_description = load_pfam_description();

//load HGNC data
function load_hgnc_db()
{
	$output = array();
	
	//parse approved genes
	$filename = get_path("data_folder")."/dbs/HGNC/hgnc_complete_set.tsv";
	foreach (file($filename) as $line)
	{
		$line = trim($line);
		if ($line=="" || starts_with($line, "hgnc_id\t")) continue;
		list($id, $symbol, $name, $locus_group, $locus_type, $status) = explode("\t", $line);
		
		$id = trim(strtr($id, array("HGNC:"=>"")));
		$symbol = trim($symbol);
		
		$status = trim($status);
		if ($status!="Approved") continue;
		
		$output[$id] = $symbol;
	}
	
	//try to replace withdrawn symbols by current symbols
	$withdrawn = array();
	$filename = get_path("data_folder")."/dbs/HGNC/hgnc_withdrawn.tsv";
	foreach (file($filename) as $line)
	{
		$line = nl_trim($line);
		if ($line=="" || starts_with($line, "HGNC_ID\t")) continue;
		list($id, $status, $symbol, $merged_into) = explode("\t", $line);
		
		//skip all but approved merged
		$status = trim($status);
		if ($status!="Merged/Split") continue;		
		if(contains($merged_into, ",")) continue; 
		if(!contains($merged_into, "Approved")) continue;
		
		$id = trim(strtr($id, array("HGNC:"=>"")));
		$id_new = explode("|", $merged_into)[0];
		$id_new = trim(strtr($id_new, array("HGNC:"=>"")));
		$withdrawn[$id] = $id_new;
	}
	
	foreach($withdrawn as $id_old => $new_id)
	{
		if(isset($output[$new_id]))
		{
			//print "replaced: $id_old > $new_id (".$output[$new_id].")\n";
			$output[$id_old] = $output[$new_id];
		}
	}
	
	return $output;
}

$hgnc = load_hgnc_db();

//write column descriptions
$column_desc = array(
	array("filter", "Annotations for filtering and ranking variants."),
	array("quality", "Quality parameters - variant quality (QUAL), depth (DP), allele frequency (AF), mean mapping quality of alternate allele (MQM), probability of strand bias for alternate bases as phred score (SAP), probability of allele ballance as phred score (ABP)"),
	array("gene", "Affected gene list (comma-separated)."),
	array("variant_type", "Variant type."),
	array("coding_and_splicing", "Coding and splicing details (Gene, ENST number, type, impact, exon/intron number, HGVS.c, HGVS.p, Pfam domain)."),
	array("coding_and_splicing_refseq", "Coding and splicing details based on RefSeq (Gene, RefSeq id, type, impact, exon/intron number, HGVS.c, HGVS.p, Pfam domain)."),
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
	array("phyloP", "phyloP (100way vertebrate) annotation. Deleterious threshold > 1.6."),
	array("Sift", "Sift effect prediction and score for each transcript: D=damaging, T=tolerated."),
	array("PolyPhen", "PolyPhen (humVar) effect prediction and score for each transcript: D=probably damaging, P=possibly damaging, B=benign."),
	array("CADD", "CADD pathogenicity prediction scores (scaled phred-like). Deleterious threshold > 10-20."),
	array("REVEL", "REVEL pathogenicity prediction score. Deleterious threshold > 0.5."),
	array("MaxEntScan", "MaxEntScan splicing prediction (reference bases score/alternate bases score)."),
	array("COSMIC", "COSMIC somatic variant database anntotation."),
	array("SpliceAI", "SpliceAI prediction of splice-site variations. Probability of the variant being splice-altering (range from 0-1). The score is the maximum value of acceptor/donor gain/loss of all effected genes."),
	array("PubMed", "PubMed ids to publications on the given variant.")
);

// optional NGSD somatic header description if vcf contains NGSD somatic information
$column_desc_ngsd_som = array(
	array("NGSD_som_c", "Somatic variant count in the NGSD."),
	array("NGSD_som_p", "Project names of project containing this somatic variant in the NGSD."),
	array("NGSD_som_vicc_interpretation", "Somatic variant interpretation according VICC standard in the NGSD."),
	array("NGSD_som_vicc_comment", "Somatic VICC interpretation comment in the NGSD.")
);

// optional NGSD header description if vcf contains NGSD information
$column_desc_ngsd = array(
	array("NGSD_hom", "Homozygous variant count in NGSD."),
	array("NGSD_het", "Heterozygous variant count in NGSD."),
	array("NGSD_group", "Homozygous / heterozygous variant count in NGSD with the same disease group."),
	array("classification", "Classification from the NGSD."),
	array("classification_comment", "Classification comment from the NGSD."),
	array("validation", "Validation information from the NGSD. Validation results of other samples are listed in brackets!"),
	array("comment", "Variant comments from the NGSD."),
	array("gene_info", "Gene information from NGSD (inheritance mode, gnomAD o/e scores).")
);

//optional COSMIC CMC header description if vcf contains COSMIC CMC information
$colum_desc_cosmic_cmc = array(
	array("CMC_genes", "Gene symbol from COSMIC Cancer Mutation Census (CMC)."),
	array("CMC_MUTATION_ID", "COSV identifier of variant from COSMIC Cancer Mutation Census (CMC)."),
	array("CMC_disease", "diseases with > 1% of samples mutated where disease = Primary site(tissue) / Primary histology / Sub-histology = Samples mutated / Samples tested = Frequency from COSMIC Cancer Mutation Census (CMC)."),
	array("CMC_DNDS_ratio", "diseases with significant dn/ds ratio (q-value < 5%) from COSMIC Cancer Mutation Census (CMC)."),
	array("CMC_mutation_significance", "Significance tier of the variant from COSMIC Cancer Mutation Census (CMC).")
);

//optional CANCERHOTSPOTS header description if vcf contains Cancerhotspots.org information
$column_desc_cancerhotspots = array(
	array("CANCERHOTSPOTS_AA_CHANGE", "Amino acid change as in original cancerhotspots.org file"),
	array("CANCERHOTSPOTS_TOTAL_MUT", "Total mutation count in cancerhotspots.org at certain amino acid position."),
	array("CANCERHOTSPOTS_ALT_COUNT", "Count of specific amino acid alteration at same position in cancerhotspots.org.")
);


if ($genotype_mode=="single")
{
	array_unshift($column_desc, array("genotype", "Genotype of variant in sample."));	
}

//write filter descriptions
$filter_desc = array();
$filter_desc[] = array("low_conf_region", "Low confidence region for small variant calling based on gnomAD AC0/RF filters and IMGAG trio/twin data.");

//parse input
$c_written = 0;
$c_skipped_wgs = 0;
$multi_cols = array();
$hgnc_messages = array();
$in_header = true;
$handle = fopen2($in, "r");
$handle_out = fopen2($out, "w");
fwrite($handle_out, "##GENOME_BUILD=GRCh38\n");
$skip_ngsd = true; // true as long as no NGSD header is found
$skip_ngsd_som = true; // true as long as no NGSD somatic header is found
$skip_cosmic_cmc = true; //true as long as no COSMIC Cancer Mutation Census (CMC) header is found.
$skip_cancerhotspots = true; //true as long as no CANCERHOTSPOTS header is found.

//write date (of input file)
fwrite($handle_out, "##CREATION_DATE=".date("Y-m-d", filemtime($in))."\n");

while(!feof($handle))
{
	$line = nl_trim(fgets($handle));
	if ($line=="" || trim($line)=="") continue;
	
	//write header
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
			list($name) = explode(",", substr($line, 13, -1));
			if ($genotype_mode=="single")
			{
				if ($column_desc[0][0]!="genotype")
				{
					trigger_error("Several sample header lines in 'single' mode!", E_USER_ERROR);
				}
				$column_desc[0][0] = $name;
			}
			else if ($genotype_mode=="multi")
			{
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
			$i_hgnc_id = index_of($cols, "HGNC_ID", false);
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
			$i_revel = index_of($cols, "REVEL", false);
			$i_polyphen = index_of($cols, "PolyPhen");
			$i_existingvariation = index_of($cols, "Existing_variation");
			$i_af_kg = index_of($cols, "AF");
			$i_af_gnomad = index_of($cols, "gnomAD_AF");
			$i_af_gnomad_afr = index_of($cols, "gnomAD_AFR_AF");
			$i_af_gnomad_amr = index_of($cols, "gnomAD_AMR_AF");
			$i_af_gnomad_eas = index_of($cols, "gnomAD_EAS_AF");
			$i_af_gnomad_nfe = index_of($cols, "gnomAD_NFE_AF");
			$i_af_gnomad_sas = index_of($cols, "gnomAD_SAS_AF");
			$i_maxes_ref = index_of($cols, "MaxEntScan_ref");
			$i_maxes_alt = index_of($cols, "MaxEntScan_alt");
			$i_pubmed = index_of($cols, "PUBMED"); 
		}

		//get annotation indices in CSQ_refseq field
		if (starts_with($line, "##INFO=<ID=CSQ_refseq,"))
		{
			$cols = explode("|", substr($line, 0, -2));
			$cols[0] = "Allele";
			$i_consequence_refseq = index_of($cols, "Consequence");
			$i_impact_refseq = index_of($cols, "IMPACT");
			$i_symbol_refseq = index_of($cols, "SYMBOL");
			$i_hgnc_id_refseq = index_of($cols, "HGNC_ID", false);
			$i_feature_refseq = index_of($cols, "Feature");
			$i_featuretype_refseq = index_of($cols, "Feature_type");
			$i_exon_refseq = index_of($cols, "EXON");
			$i_intron_refseq = index_of($cols, "INTRON");
			$i_hgvsc_refseq = index_of($cols, "HGVSc");
			$i_hgvsp_refseq = index_of($cols, "HGVSp");
			$i_domains_refseq = index_of($cols, "DOMAINS");
		}


		// detect NGSD header lines
		if (starts_with($line, "##INFO=<ID=NGSD_"))
		{
			$skip_ngsd = false;
			if (starts_with($line, "##INFO=<ID=NGSD_SOM_"))
			{
				$skip_ngsd_som = false;
			}
		}
		
		//COSMIC CMC (Cancer Mutation Census) header line
		if (starts_with($line, "##INFO=<ID=COSMIC_CMC,") )
		{
			$skip_cosmic_cmc = false;
		
			//trim to " (from file "... and split by "|"
			$cols = explode("|", substr($line, 0, strpos($line, " (from file")) );
			$cols[0] = "GENE_NAME";
			
			$i_cosmic_cmc_gene_name = index_of($cols, "GENE_NAME");
			$i_cosmic_cmc_mut_id = index_of($cols, "GENOMIC_MUTATION_ID");
			$i_cosmic_cmc_disease = index_of($cols, "DISEASE");
			$i_cosmic_cmc_dnds_disease = index_of($cols, "DNDS_DISEASE_QVAL");
			$i_cosmic_cmc_mut_sign_tier = index_of($cols, "MUTATION_SIGNIFICANCE_TIER");
		}
		
		//Cancerhotspots.org header line
		if( starts_with($line, "##INFO=<ID=CANCERHOTSPOTS,") )
		{
			$skip_cancerhotspots = false;
			$cols = explode("|", substr($line, 0,-2) );
			$cols[0] = "GENE_SYMBOL";
			
			$i_cancerhotspots_gene_symbol = index_of($cols, "GENE_SYMBOL");
			$i_cancerhotspots_transcript_id = index_of($cols, "ENSEMBL_TRANSCRIPT");
			$i_cancerhotspots_aa_pos = index_of($cols, "AA_POS");
			$i_cancerhotspots_aa_ref = index_of($cols, "AA_REF");
			$i_cancerhotspots_aa_alt = index_of($cols, "AA_ALT");
			$i_cancerhotspots_total_count = index_of($cols, "TOTAL_COUNT");
			$i_cancerhotspots_alt_count = index_of($cols, "ALT_COUNT");
		}
		
		continue;
	}
	//after last header line, write our header
	else if ($in_header) 
	{
		if (!$skip_ngsd_som)
		{
			// append optional NGSD somatic header
			$column_desc = array_merge($column_desc, $column_desc_ngsd_som);
		}
		if (!$skip_ngsd)
		{
			// append optional NGSD header
			$column_desc = array_merge($column_desc, $column_desc_ngsd);
		}
		if (!$skip_cosmic_cmc)
		{
			//append optional COSMIC CMC header
			$column_desc = array_merge($column_desc, $colum_desc_cosmic_cmc);
		}
		
		if( !$skip_cancerhotspots )
		{
			//append optional cancerhotspots header
			$column_desc = array_merge($column_desc, $column_desc_cancerhotspots);
		}
		
		write_header_line($handle_out, $column_desc, $filter_desc);
		$in_header = false;
	}
	
	//write content lines
	$cols = explode("\t", $line);
	if (count($cols)<10) trigger_error("VCF file line contains less than 10 columns: '$line'", E_USER_ERROR);
	list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $sample) = $cols;
	$tag = "{$chr}:{$pos} {$ref}>{$alt}";
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
		$genotypes = array();
		$depths = array();
		$aos = array();
		foreach($multi_cols as $col)
		{
			$gt = $tmp[$col];
			$dp = $tmp2[$col];
			$ao = $tmp3[$col];
			if ($dp<3) $gt = "n/a";
			$genotypes[] = $gt;
			$depths[] = $dp;
			$aos[] = $ao;
		}
		$genotype = "\t".implode("\t", $genotypes);
		$sample["DP"] = implode(",", $depths);
		$sample["AO"] = implode(",", $aos);
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
	if (isset($info["SAP"])) 
	{
		$quality[] = "SAP=".intval($info["SAP"]);
	}
	if (isset($info["ABP"])) 
	{
		$quality[] = "ABP=".intval($info["ABP"]);
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
	$hom_gnomad = array();
	$hemi_gnomad = array();
	$clinvar = array();
	$hgmd = array();
	$maxentscan = array();
	$regulatory = array();
	$pubmed = array();
	
	//variant details (up/down-stream)
	$variant_details_updown = array();
	$genes_updown = array();
	$sift_updown = array();
	$polyphen_updown = array();
	$coding_and_splicing_details_updown = array();

	// RefSeq transcript annotation
	$coding_and_splicing_details_refseq = array();
	$coding_and_splicing_details_updown_refseq = array();

	
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
			$af_gnomad_afr[] = trim($parts[$i_af_gnomad_afr]);
			$af_gnomad_amr[] = trim($parts[$i_af_gnomad_amr]);
			$af_gnomad_eas[] = trim($parts[$i_af_gnomad_eas]);
			$af_gnomad_nfe[] = trim($parts[$i_af_gnomad_nfe]);
			$af_gnomad_sas[] = trim($parts[$i_af_gnomad_sas]);
			
			//dbSNP, COSMIC
			$ids = explode("&", $parts[$i_existingvariation]);
			foreach($ids as $id)
			{
				if (starts_with($id, "rs"))
				{
					$dbsnp[] = $id;
				}
				if (starts_with($id, "COSM") || starts_with($id, "COSN") || starts_with($id, "COSV"))
				{
					$cosmic[] = $id;
				}
			}
			
			//MaxEntScan
			if ($parts[$i_maxes_ref]!="")
			{
				$result = number_format($parts[$i_maxes_ref], 2).">".number_format($parts[$i_maxes_alt], 2);
				$maxentscan[] = $result;
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
			$revel_score = trim($parts[$i_revel]);
			if ($revel_score!="") $revel[] = $revel_score;

			//PubMed ids
			$pubmed = array_merge($pubmed, explode("&", $parts[$i_pubmed]));

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
				
				//determine gene name (update if neccessary)
				$gene = trim($parts[$i_symbol]);
				if ($i_hgnc_id!==false && $gene!="")
				{
					$hgnc_id = $parts[$i_hgnc_id];
					$hgnc_id = trim(strtr($hgnc_id, array("HGNC:"=>"")));
					if (isset($hgnc[$hgnc_id]))
					{
		
						$hgnc_gene = $hgnc[$hgnc_id];
						if ($gene!=$hgnc_gene)
						{
							//@$hgnc_messages["gene name '$gene' with ID '$hgnc_id' replaced by '$hgnc_gene'"] += 1;
							$gene = $hgnc_gene;
						}
					}
					else if ($hgnc_id!="")
					{
						@$hgnc_messages["ID '$hgnc_id' not valid or withdrawn for gene '$gene'"] += 1;
					}
				}
				if (!$is_updown)
				{
					$genes[] = $gene;
				}
				else
				{
					$genes_updown[] = $gene;
				}
				
				//pathogenicity predictions (transcript-specific)
				list($sift_type, $sift_score) = explode("(", strtr($parts[$i_sift], ")", "(")."(");
				if ($sift_type!="")
				{
					$sift_type = translate("Sift", $sift_type, array("deleterious"=>"D", "tolerated"=>"T", "tolerated_low_confidence"=>"T", "deleterious_low_confidence"=>"D"));
					$sift_score = number_format($sift_score, 2);
					$sift_entry = $sift_type."($sift_score)";
				}
				else
				{
					$sift_entry = " ";
				}
				if (!$is_updown)
				{
					$sift[] = $sift_entry;
				}
				else
				{
					$sift_updown[] = $sift_entry;
				}
				
				list($pp_type, $pp_score) = explode("(", strtr($parts[$i_polyphen], ")", "(")."(");
				if ($pp_type!="")
				{
					$pp_type = translate("PolyPhen", $pp_type, array("unknown"=>" ",  "probably_damaging"=>"D", "possibly_damaging"=>"P", "benign"=>"B"));
					$pp_score = number_format($pp_score, 2);
					if (!is_numeric($pp_score))
					{
						print "ERROR: ".$parts[$i_polyphen]." ".$pp_score."\n";
					}
					$polyphen_entry = $pp_type."($pp_score)";
				}
				else
				{
					$polyphen_entry = " ";
				}
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
				if ($intron!="") $intron = "intron".$intron;
				
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
					if(starts_with($entry, "Pfam:"))
					{
						$domain = explode(":", $entry, 2)[1];
					}
				}

				// extend domain ID by description
				$domain_description = "";
				if ($domain != "")
				{
					// update Pfam ID 
					if (array_key_exists($domain, $pfam_replacements))
					{
						if ($pfam_replacements[$domain] == "")
						{
							$domain_description = "removed";
						}
						else
						{
							$domain_description = "(new id of $domain) ";
							$domain = $pfam_replacements[$domain];
						}
					}
					// append description
					if (array_key_exists($domain, $pfam_description))
					{
						$domain_description = $domain_description.$pfam_description[$domain];
					}

					// throw error if Pfam id is neither found in replacement data nor in description data
					if ($domain_description == "")
					{
						trigger_error("No description found for '$domain'!", E_USER_WARNING);
					}

					// combine decription and id
					$domain = "$domain [$domain_description]";
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

	// annotation of RefSeq transcripts
	if (isset($info["CSQ_refseq"]))
	{
		$anns = explode(",", $info["CSQ_refseq"]);
		foreach($anns as $entry)
		{			
			$parts = explode("|", $entry);
			
			//only transcripts
			$feature_type = $parts[$i_featuretype_refseq];
			if ($feature_type=="Transcript" || $feature_type=="")
			{
				$transcript_id = $parts[$i_feature_refseq];
				
				//extract variant type
				$variant_type = strtr($parts[$i_consequence_refseq], array("_variant"=>""));
				$variant_type = strtr($variant_type, array("splice_acceptor&splice_region&intron"=>"splice_acceptor", "splice_donor&splice_region&intron"=>"splice_donor", "splice_acceptor&intron"=>"splice_acceptor", "splice_donor&intron"=>"splice_donor", "_prime_"=>"'"));
				$is_updown = $variant_type=="upstream_gene" || $variant_type=="downstream_gene";
				
				//determine gene name (update if neccessary)
				$gene = trim($parts[$i_symbol_refseq]);
				
				//exon
				$exon = trim($parts[$i_exon_refseq]);
				if ($exon!="") $exon = "exon".$exon;
				$intron = trim($parts[$i_intron_refseq]);
				if ($intron!="") $intron = "intron".$intron;
				
				//hgvs
				$hgvs_c = trim($parts[$i_hgvsc_refseq]);
				if ($hgvs_c!="") $hgvs_c = explode(":", $hgvs_c)[1];			
				$hgvs_p = trim($parts[$i_hgvsp_refseq]);
				if ($hgvs_p!="") $hgvs_p = explode(":", $hgvs_p)[1];
				$hgvs_p = str_replace("%3D", "=", $hgvs_p);
				
				//domain
				$domain = "";
				$domains = explode("&", $parts[$i_domains_refseq]);
				foreach($domains as $entry)
				{
					if(starts_with($entry, "Pfam_domain:"))
					{
						$domain = explode(":", $entry, 2)[1];
					}
				}

				// extend domain ID by description
				$domain_description = "";
				if ($domain != "")
				{
					// update Pfam ID 
					if (array_key_exists($domain, $pfam_replacements))
					{
						if ($pfam_replacements[$domain] == "")
						{
							$domain_description = "removed";
						}
						else
						{
							$domain_description = "(new id of $domain) ";
							$domain = $pfam_replacements[$domain];
						}
					}
					// append description
					if (array_key_exists($domain, $pfam_description))
					{
						$domain_description = $domain_description.$pfam_description[$domain];
					}

					// throw error if Pfam id is neither found in replacement data nor in description data
					if ($domain_description == "")
					{
						trigger_error("No description found for '$domain'!", E_USER_WARNING);
					}

					// combine decription and id
					$domain = "$domain [$domain_description]";
				}
				
				$transcript_entry = "{$gene}:{$transcript_id}:".$parts[$i_consequence_refseq].":".$parts[$i_impact_refseq].":{$exon}{$intron}:{$hgvs_c}:{$hgvs_p}:{$domain}";
				if (!$is_updown)
				{
					$coding_and_splicing_details_refseq[] = $transcript_entry;
				}
				else
				{
					$coding_and_splicing_details_updown_refseq[] = $transcript_entry;
				}
			
			}
		}
	}

	// gnomAD_genome
	if (isset($info["gnomADg_AF"]))
	{
		$gnomad_value = trim($info["gnomADg_AF"]);
		if (strpos($gnomad_value, "&") !== false)
		{
			// special handling of the rare case that 2 gnomAD AF values exist for this variant
			$gnomad_values = explode("&", $gnomad_value);
			foreach($gnomad_values as $value)
			{
				if ($value != ".")
				{
					$af_gnomad_genome[] = $value;
				}
			}		
		}
		else
		{
			if ($gnomad_value != ".")
			{
				$af_gnomad_genome[] = $gnomad_value;
			}
			
		}
		
		
		
	}
	if (isset($info["gnomADg_Hom"])) $hom_gnomad[] = trim($info["gnomADg_Hom"]);
	if (isset($info["gnomADg_Hemi"])) $hemi_gnomad[] = trim($info["gnomADg_Hemi"]);

	//ClinVar
	if (isset($info["CLINVAR_ID"]) && isset($info["CLINVAR_DETAILS"]))
	{
		$clin_accs = explode("&", trim($info["CLINVAR_ID"]));
		$clin_details = explode("&", trim($info["CLINVAR_DETAILS"]));
		for ($i=0; $i<count($clin_accs); ++$i)
		{
			if ($clin_accs[$i]!="")
			{
				$clinvar[] = $clin_accs[$i]." [".strtr(vcf_decode_url_string($clin_details[$i]), array("_"=>" "))."];";
			}
		}
	}
	

	//HGMD
	if (isset($info["HGMD_ID"]))
	{
		$hgmd_id = trim($info["HGMD_ID"]);
		$hgmd_class = trim($info["HGMD_CLASS"]);
		$hgmd_mut = trim($info["HGMD_MUT"]);
		$hgmd_gene = trim($info["HGMD_GENE"]);
		$hgmd_phen = trim($info["HGMD_PHEN"]);
		
		//TODO: Possible to have more than one match per variant?
		$text = $hgmd_id." [CLASS=".$hgmd_class." MUT=".$hgmd_mut." PHEN=".strtr($hgmd_phen, "_", " ")." GENE=".$hgmd_gene."]; ";
		$hgmd[] = trim($text);
	}

	//AFs
	$dbsnp = implode(",", collapse($tag, "dbSNP", $dbsnp, "unique"));	
	$kg = collapse($tag, "1000g", $af_kg, "one", 4);
	$gnomad = collapse($tag, "gnomAD", $af_gnomad, "one", 4);
	$gnomad_genome = collapse($tag, "gnomAD genome", $af_gnomad_genome, "max", 4);
	$gnomad = max($gnomad, $gnomad_genome);
	$gnomad_hom_hemi = collapse($tag, "gnomAD Hom", $hom_gnomad, "one").",".collapse($tag, "gnomAD Hemi", $hemi_gnomad, "one");
	if ($gnomad_hom_hemi==",") $gnomad_hom_hemi = "";
	$gnomad_sub = collapse($tag, "gnomAD AFR", $af_gnomad_afr, "one", 4).",".collapse($tag, "gnomAD AMR", $af_gnomad_amr, "one", 4).",".collapse($tag, "gnomAD EAS", $af_gnomad_eas, "one", 4).",".collapse($tag, "gnomAD NFE", $af_gnomad_nfe, "one", 4).",".collapse($tag, "gnomAD SAS", $af_gnomad_sas, "one", 4);
	if (str_replace(",", "", $gnomad_sub)=="") $gnomad_sub = "";

	//PubMed
	$pubmed = implode(",", collapse($tag, "PubMed", $pubmed, "unique"));	


	if (!$skip_ngsd)
	{
		// extract NGSD somatic counts
		if (!$skip_ngsd_som)
		{
			if (isset($info["NGSD_SOM_C"]))
			{
				$ngsd_som_counts = intval(trim($info["NGSD_SOM_C"]));
			}
			else
			{
				$ngsd_som_counts = "0";
			}

			if (isset($info["NGSD_SOM_P"]))
			{
				$ngsd_som_projects = vcf_decode_url_string(trim($info["NGSD_SOM_P"]));
			}
			else
			{
				$ngsd_som_projects = "";
			}
			
			if (isset($info["NGSD_SOM_VICC"]))
			{
				$ngsd_som_vicc = trim($info["NGSD_SOM_VICC"]);
			}
			else
			{
				$ngsd_som_vicc = "";
			}
			
			if (isset($info["NGSD_SOM_VICC_COMMENT"]))
			{
				$ngsd_som_vicc_comment = vcf_decode_url_string(trim($info["NGSD_SOM_VICC_COMMENT"]) );
			}
			else
			{
				$ngsd_som_vicc_comment = "";
			}
		}

		//NGSD
		if (isset($info["NGSD_HAF"]) || $gnomad >= 0.05 || $kg >= 0.05)
		{
			$ngsd_hom = "n/a (AF>5%)";
			$ngsd_het = "n/a (AF>5%)";
			$ngsd_group = "n/a (AF>5%)";
		}
		elseif(isset($info["NGSD_COUNTS"]))
		{
			$ngsd_counts = explode(",", trim($info["NGSD_COUNTS"]));
			$ngsd_hom = $ngsd_counts[0];
			$ngsd_het = $ngsd_counts[1];
			if (isset($info["NGSD_GROUP"]))
			{
				$ngsd_group_raw = explode(",", trim($info["NGSD_GROUP"]));
				$ngsd_group = intval($ngsd_group_raw[0])." / ".intval($ngsd_group_raw[1]);
			}
			else
			{
				$ngsd_group = "0 / 0";
			}
		}
		else
		{
			$ngsd_hom = "0";
			$ngsd_het = "0";
			$ngsd_group = "0 / 0";
		}

		if (isset($info["NGSD_CLAS"]))
		{
			$ngsd_clas = trim($info["NGSD_CLAS"]);
		}
		else
		{
			$ngsd_clas = "";
		}

		if (isset($info["NGSD_CLAS_COM"]))
		{
			$ngsd_clas_com = vcf_decode_url_string(trim($info["NGSD_CLAS_COM"]));
		}
		else
		{
			$ngsd_clas_com = "";
		}

		if (isset($info["NGSD_COM"]))
		{
			$ngsd_com = vcf_decode_url_string(trim($info["NGSD_COM"]));
		}
		else
		{
			$ngsd_com = "";
		}

		if (isset($info["NGSD_VAL"]))
		{
			$ngsd_val = trim($info["NGSD_VAL"]);
		}
		else
		{
			$ngsd_val = "";
		}

		if (isset($info["NGSD_GENE_INFO"]))
		{
			$ngsd_gene_info = trim(str_replace("&", ", ", vcf_decode_url_string($info["NGSD_GENE_INFO"])));
		}
		else
		{
			$ngsd_gene_info = "";
		}
	}
	
	//SpliceAI
	$spliceai = "";
	if (isset($info["SpliceAI"]))
	{
		$splice_number = null;
		$spliceai_info = trim($info["SpliceAI"]);
		$spliceai_values = array();

		$entries = explode(",", $spliceai_info);
		foreach($entries as $entry)
		{
			$delta_scores = explode("|", $entry);
			if(sizeof($delta_scores) == 10)
			{
				$tmp_score = max(floatval($delta_scores[2]), floatval($delta_scores[3]), floatval($delta_scores[4]), floatval($delta_scores[5]));
				if(is_null($splice_number)) $splice_number = $tmp_score;
				$splice_number = max($splice_number, $tmp_score);
			}
			else
			{
				trigger_error("Wrong SpliceAI annotation in line: ${line} in SpliceAI annotation: ${spliceai_info}! Delimiter for several genes must be ','.", E_USER_WARNING);
			}
		}

		if(!is_null($splice_number))
		{
			$spliceai = $splice_number;
		}

	}

	// CADD
	$cadd_scores = array();
	if (isset($info["CADD_SNV"]))
	{
		$cadd_scores = array_map(function($score){return number_format($score, 2, ".", "");}, explode("&", $info["CADD_SNV"]));
	}
	if (isset($info["CADD_INDEL"]))
	{
		$cadd_scores = array_map(function($score){return number_format($score, 2, ".", "");}, explode("&", $info["CADD_INDEL"]));
	}
	if (count(array_unique($cadd_scores)) == 0)
	{
		//No CADD score available
		$cadd = "";
	}
	else if (count(array_unique($cadd_scores)) > 1)
	{
		//trigger_error("Multiple values for CADD score for variant $chr:$pos! Choosing max value.", E_USER_WARNING);
		$cadd = max($cadd_scores);
	}
	else
	{
		$cadd = $cadd_scores[0];
	}
	
	
	// COSMIC CMC
	if ( !$skip_cosmic_cmc && isset($info["COSMIC_CMC"]) )
	{
		$anns = explode("&", $info["COSMIC_CMC"]);
		
		$cmc_gene = array();
		$cmc_mut_id = array();
		$cmc_disease = array();
		$cmc_dnds_disease = array();
		$cmc_mut_sign_tier = array();
		

		foreach($anns as $entry)
		{
			$parts = explode("|", $entry);
			
			$cmc_gene[] = vcf_decode_url_string( $parts[$i_cosmic_cmc_gene_name] );
			$cmc_mut_id[] = vcf_decode_url_string( $parts[$i_cosmic_cmc_mut_id] ) ;
			$cmc_disease[] = vcf_decode_url_string( $parts[$i_cosmic_cmc_disease] );
			$cmc_dnds_disease[] = vcf_decode_url_string( $parts[$i_cosmic_cmc_dnds_disease] );
			$cmc_mut_sign_tier[] = vcf_decode_url_string( $parts[$i_cosmic_cmc_mut_sign_tier] );
			
		}
	}
	
	// CANCERHOTSPOTS
	if( !$skip_cancerhotspots && isset($info["CANCERHOTSPOTS"]) )
	{
		$anns = explode(",", $info["CANCERHOTSPOTS"] );
		
		$cancerhotspots_protein_change = array();
		$cancerhotspots_total_count = array();
		$cancerhotspots_alt_count = array();
		
		foreach($anns as $entry)
		{
			$parts = explode("|", $entry);
			$cancerhotspots_protein_change[] = ($parts[$i_cancerhotspots_transcript_id]) .":p." . aa1_to_aa3($parts[$i_cancerhotspots_aa_ref]) .$parts[$i_cancerhotspots_aa_pos] .  aa1_to_aa3($parts[$i_cancerhotspots_aa_alt]);
			$cancerhotspots_total_count[] =  $parts[$i_cancerhotspots_total_count];
			$cancerhotspots_alt_count[] = $parts[$i_cancerhotspots_alt_count];
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
	
	$variant_details = implode(",", array_unique($variant_details));
	$coding_and_splicing_details =  implode(",", $coding_and_splicing_details);
	
	//regulatory
	$regulatory = implode(",", collapse($tag, "Regulatory", $regulatory, "unique"));

	//RepeatMasker
	$repeatmasker = "";
	if (isset($info["REPEATMASKER"]))
	{
		$repeatmasker = trim(str_replace("&", ", ", vcf_decode_url_string($info["REPEATMASKER"])));
	}
	
	//effect predicions
	$phylop = collapse($tag, "phyloP", $phylop, "one", 4);
	$sift = implode(",", $sift);
	if (trim(strtr($sift, ",", " "))=="") $sift = "";
	$polyphen = implode(",", $polyphen);
	if (trim(strtr($polyphen, ",", " "))=="") $polyphen = "";
	$fathmm = collapse($tag, "fathmm-MKL", $fathmm, "one");
	
	$revel = empty($revel) ? "" : collapse($tag, "REVEL", $revel, "max", 2);
	
	//OMIM
	$omim = "";
	if (isset($info["OMIM"]))
	{
		$omim = trim(vcf_decode_url_string($info["OMIM"]));
	}

	//ClinVar
	$clinvar = implode(" ", collapse($tag, "ClinVar", $clinvar, "unique"));
	
	//HGMD
	$hgmd = collapse($tag, "HGMD", $hgmd, "one");

	//MaxEntScan
	$maxentscan = implode(",", collapse($tag, "MaxEntScan", $maxentscan, "unique"));
	
	//COSMIC
	$cosmic = implode(",", collapse($tag, "COSMIC", $cosmic, "unique"));
	
	//skip common MODIFIER variants in WGS mode
	if ($wgs && skip_in_wgs_mode($chr, $coding_and_splicing_details, $kg, $gnomad, $clinvar, $hgmd, $ngsd_clas))
	{
		++$c_skipped_wgs;
		continue;
	}
	
	//write data
	++$c_written;
	$genes = array_unique($genes);
	fwrite($handle_out, "$chr\t$start\t$end\t$ref\t{$alt}{$genotype}\t".implode(";", $filter)."\t".implode(";", $quality)."\t".implode(",", $genes)."\t$variant_details\t$coding_and_splicing_details\t".implode(",", $coding_and_splicing_details_refseq)."\t$regulatory\t$omim\t$clinvar\t$hgmd\t$repeatmasker\t$dbsnp\t$kg\t$gnomad\t$gnomad_hom_hemi\t$gnomad_sub\t$phylop\t$sift\t$polyphen\t$cadd\t$revel\t$maxentscan\t$cosmic\t$spliceai\t$pubmed");
	if (!$skip_ngsd_som)
	{
		fwrite($handle_out, "\t$ngsd_som_counts\t$ngsd_som_projects\t$ngsd_som_vicc\t$ngsd_som_vicc_comment");
	}
	if (!$skip_ngsd)
	{
		fwrite($handle_out, "\t$ngsd_hom\t$ngsd_het\t$ngsd_group\t$ngsd_clas\t$ngsd_clas_com\t$ngsd_val\t$ngsd_com\t$ngsd_gene_info");
	}
	
	if ( !$skip_cosmic_cmc && isset($info["COSMIC_CMC"]) )
	{
		fwrite($handle_out, "\t". implode(",",$cmc_gene) ."\t". implode(",",$cmc_mut_id) ."\t". implode(",",$cmc_disease) ."\t" . implode(",",$cmc_dnds_disease) ."\t".implode(",",$cmc_mut_sign_tier));
	}
	elseif( !$skip_cosmic_cmc)
	{
		fwrite( $handle_out, str_repeat("\t", count($colum_desc_cosmic_cmc)) );
	}
	
	if( !$skip_cancerhotspots &&isset($info["CANCERHOTSPOTS"]) )
	{
		fwrite( $handle_out, "\t" . implode(",",  $cancerhotspots_protein_change) . "\t" . implode(",", $cancerhotspots_total_count) ."\t" . implode(",", $cancerhotspots_alt_count) ); 
	}
	elseif( !$skip_cancerhotspots )
	{
		fwrite( $handle_out, str_repeat("\t", count($column_desc_cancerhotspots) ) );
	}

	fwrite($handle_out, "\n");
}

//print HGNC messages
foreach($hgnc_messages as $message => $c)
{
	$parser->log("HGNC: {$message} ({$c}x)");
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
