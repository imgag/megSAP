<?php 
/** 
	@page vcf2gsvar
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("vcf2gsvar", "Converts an annotated VCF file from freebayes to a GSvar file.");
$parser->addInfile("in",  "Input file in VCF or VCF.GZ format.", false);
$parser->addOutfile("out", "Output file in GSvar format.", false);
//optional
$parser->addEnum("genotype_mode", "Genotype handling mode.", true, array("single", "multi", "skip"), "single");
$parser->addFlag("updown", "Don't discard up- or downstream annotations (5000 bases around genes).");
$parser->addFlag("wgs", "Enables WGS mode: MODIFIER variants with a AF>2% are skipped to reduce the number of variants to a manageable size.");
$parser->addFlag("longread", "Add additional columns for long-read samples (e.g. pahsing information)");
$parser->addString("custom", "Settings key name for custom column definitions.", true, "");
$parser->addFlag("test", "Run in test mode. Skips replacing sample headers with NGSD information.");
extract($parser->parse($argv));

$custom_columns = [];
if ($custom!="")
{
	$tmp = get_path($custom, false);
	if (is_array($tmp)) $custom_columns = $tmp;
}

//skip common MODIFIER variants in WGS mode
function skip_in_wgs_mode($chr, $coding_and_splicing_details, $gnomad, $clinvar, $hgmd, $ngsd_clas)
{	
	//don't skip mito variants
	if ($chr=='chrMT') return false;
	
	//don't skip exonic/splicing variants
	if (contains($coding_and_splicing_details, ":LOW:") || contains($coding_and_splicing_details, ":MODERATE:") || contains($coding_and_splicing_details, ":HIGH:")) return false;
	
	//don't skip variants annotated to be (likely) pathogenic
	if (contains($hgmd, "CLASS=DM") || (contains($clinvar, "pathogenic") && !contains($clinvar, "conflicting"))) return false;	
	
	//don't skip variants of class 4/5/M in NGSD
	if ($ngsd_clas=='4' || $ngsd_clas=='5'|| $ngsd_clas=='M') return false;
	
	//skip common variants >2%AF
	if ($gnomad!="" && $gnomad>0.02) return true;
	
	return false; //non-exonic but rare
}

//get index of columnn in QSC header.
function index_of($cols, $name, $label, $optonal = true)
{
	$index = array_search($name, $cols);
	if ($index===FALSE)
	{
		if ($optonal) return -1;
		trigger_error("Could not find column '$name' in annotation field '{$label}'. Valid column names are: ".implode(", ", array_values($cols)), E_USER_ERROR);
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
	//remove '' and '.'
	$values = array_map('trim', $values);
	$tmp = [];
	foreach($values as $v)
	{
		if ($v=="" || $v==".") continue;
		$tmp [] = $v;
	}
	$values = $tmp;
	
	//set decimal places
	if (!is_null($decimal_places))
	{
		$tmp = [];
		foreach($values as $v)
		{
			if (!is_numeric($v)) trigger_error("Invalid numeric value '{$v}' in mode '{$mode}' while collapsing '{$error_name}' in variant '{$tag}'!", E_USER_ERROR);
			
			$v_new = number_format($v, $decimal_places, ".", "");
			if ($v>0 && $v_new==0)
			{
				$v = substr($v, 0, -1)."1";
			}
			$tmp [] = $v_new;
		}
		$values = $tmp;
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
	//add optional header depending on annotation
	global $skip_ngsd_som;
	global $column_desc_ngsd_som;
	if (!$skip_ngsd_som) $column_desc = array_merge($column_desc, $column_desc_ngsd_som);
	global $skip_ngsd;
	global $column_desc_ngsd;
	if (!$skip_ngsd) $column_desc = array_merge($column_desc, $column_desc_ngsd);
	global $skip_cosmic_cmc;
	global $colum_desc_cosmic_cmc;
	if (!$skip_cosmic_cmc) $column_desc = array_merge($column_desc, $colum_desc_cosmic_cmc);
	global $skip_cancerhotspots;
	global $column_desc_cancerhotspots;
	if(!$skip_cancerhotspots) $column_desc = array_merge($column_desc, $column_desc_cancerhotspots);
	global $skip_short_read_overlap_annotation;
	global $column_desc_short_read_overlap_annotation;
	if(!$skip_short_read_overlap_annotation) $column_desc = array_merge($column_desc, $column_desc_short_read_overlap_annotation);
	global $column_desc_custom;
	$column_desc = array_merge($column_desc, $column_desc_custom);
	
	//write out header
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
	$filename = get_path("data_folder")."/dbs/HGNC/hgnc_complete_set_2025-09-02.tsv";
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
	$filename = get_path("data_folder")."/dbs/HGNC/hgnc_withdrawn_2025-09-02.tsv";
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
	array("quality", "Quality parameters - variant quality (QUAL), depth (DP), quality divided by depth (QD), allele frequency (AF), mean mapping quality of alternate allele (MQM), probability of strand bias for alternate bases as phred score (SAP), probability of allele ballance as phred score (ABP)"),
	array("gene", "Affected gene list (comma-separated)."),
	array("variant_type", "Variant type."),
	array("coding_and_splicing", "Coding and splicing details (Gene, ENST number, type, impact, exon/intron number, HGVS.c, HGVS.p, Pfam domain)."),
	array("regulatory", "Regulatory consequence details."),
	array("OMIM", "OMIM database annotation."),
	array("ClinVar", "ClinVar database annotation."),
	array("HGMD", "HGMD database annotation."),
	array("RepeatMasker", "RepeatMasker annotation."),
	array("dbSNP", "Identifier in dbSNP database."),
	array("gnomAD", "Allele frequency in gnomAD project."),
	array("gnomAD_sub", "Sub-population allele frequenciens (AFR,AMR,EAS,NFE,SAS) in gnomAD project."),
	array("gnomAD_hom_hemi", "Homoyzgous/hemizygous case count of gnomAD project (genome data)."),
	array("gnomAD_het", "Heterozygous allele count of the gnomAD project (genome data)."),
	array("gnomAD_wt", "Wildtype allele count of the gnomAD project (genome data)."),
	array("phyloP", "phyloP (100way vertebrate) annotation. Deleterious threshold > 1.6."),
	array("CADD", "CADD pathogenicity prediction scores (scaled phred-like). Deleterious threshold > 10-20."),
	array("REVEL", "REVEL pathogenicity prediction score. Deleterious threshold > 0.5."),
	array("AlphaMissense", "AlphaMissense pathogenicity score. Deleterious threshold > 0.564."),
	array("MaxEntScan", "MaxEntScan reference score and alternate score for (1) native splice site, (2) acceptor gain and (3) donor gain. Comma-separated list if there are different predictions for several transcripts."),
	array("SpliceAI", "SpliceAI prediction. These include delta scores (DS) and delta positions (DP) for acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). Format: GENE|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL."),
	array("PubMed", "PubMed ids to publications on the given variant.")
);

// optional NGSD somatic header description if vcf contains NGSD somatic information
$column_desc_ngsd_som = array(
	array("NGSD_som_c", "Somatic variant count in the NGSD (tumor-normal)."),
	array("NGSD_som_p", "Project names containing this somatic variant in the NGSD (tumor-normal)."),
	array("NGSD_som_to_c", "Somatic variant count in the NGSD (tumor-only)."),
	array("NGSD_som_vicc_interpretation", "Somatic variant interpretation according VICC standard in the NGSD."),
	array("NGSD_som_vicc_comment", "Somatic VICC interpretation comment in the NGSD.")
);

// optional NGSD header description if vcf contains NGSD information
$column_desc_ngsd = array(
	array("NGSD_hom", "Homozygous variant count in NGSD."),
	array("NGSD_het", "Heterozygous variant count in NGSD."),
	array("NGSD_mosaic", "Mosaic variant count in NGSD."),
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

//optional IN_SHORTREAD_SAMPLE
$column_desc_short_read_overlap_annotation = array(
	array("in_short-read", "Variant was also found in corresponding short-read WGS sample.")
);

$column_desc_custom = [];
foreach($custom_columns as $key => $tmp)
{
	list($vcf, $col, $desc) = explode(";", $tmp, 3);
	$column_desc_custom[] = [$key, $desc];
}

if ($genotype_mode=="single")
{
	// add phasing info for long-reads
	if ($longread) array_unshift($column_desc, array("genotype_phased", "Phasing information of variant in sample. Containing phased genotype and id of the phased block."));

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
$handle = gzopen2($in, "r");
$handle_out = fopen2($out, "w");
fwrite($handle_out, "##GENOME_BUILD=GRCh38\n");
$skip_ngsd = true; // true as long as no NGSD header is found
$skip_ngsd_som = true; // true as long as no NGSD somatic header is found
$skip_cosmic_cmc = true; //true as long as no COSMIC Cancer Mutation Census (CMC) header is found.
$skip_cancerhotspots = true; //true as long as no CANCERHOTSPOTS header is found.
$skip_short_read_overlap_annotation = true; //true as long as no IN_SHORT_READ header is found.
$missing_domains = [];
$multisample_vcf = false; //determines if the input VCF is a standard multi-sample VCF or single-sample VCF/megSAP multi-sample VCF (no extra columns for each sample)
$annotate_refseq_consequences = false;

//write date (of input file)
fwrite($handle_out, "##CREATION_DATE=".date("Y-m-d", filemtime($in))."\n");

$clair_info = ["", ""];

while(!gzeof($handle))
{
	$line = nl_trim(gzgets($handle));
	if ($line=="" || trim($line)=="") continue;
	
	//write header
	if ($line[0]=="#") 
	{
		//variant calling date
		if (starts_with($line, "##fileDate="))
		{
			$date = trim(substr($line, 11));
			$date = substr_replace($date, "-", 6, 0);
			$date = substr_replace($date, "-", 4, 0);
			
			fwrite($handle_out, "##CALLING_DATE={$date}\n");
		}
		if (starts_with($line, "##DRAGENCommandLine=<ID=dragen,"))
		{
			$date = "";
			$parts = explode(",", $line);
			foreach($parts as $entry)
			{
				if (starts_with($entry, "Date="))
				{
					list(, $tmp) = explode('"', $entry);
					$tmp = explode(" ", strtr($tmp, ":", " "));
					$tmp = $tmp[1]." ".$tmp[2]." ".$tmp[7];
					$date = date("Y-m-d", strtotime($tmp));
				}
			}
			if ($date!="")
			{
				fwrite($handle_out, "##CALLING_DATE={$date}\n");
			}
		}
		
		//variant caller
		if (starts_with($line, "##source=freeBayes"))
		{
			fwrite($handle_out, "##SOURCE=freebayes ".trim(substr($line, 18))."\n");
		}
		if (starts_with($line, "##DRAGENCommandLine=<ID=dragen,"))
		{
			$dragen_ver = "";
			$parts = explode(",", strtr($line, "\"", ","));
			foreach($parts as $part)
			{
				$part = trim($part);
				if (starts_with($part, "SW:"))
				{
					$parts2 = explode(".", $part);
					$dragen_ver = implode(".", array_slice($parts2, -3));
				}
			}		
			fwrite($handle_out, "##SOURCE=DRAGEN {$dragen_ver}\n");
		}
		if (starts_with($line, "##source=Clair3"))
		{
			$clair_info[0] = "Clair3";
			if ($clair_info[0]!="" && $clair_info[1]!="") //write out clair info when complete - no matter which header is first
			{
				fwrite($handle_out, "##SOURCE=".$clair_info[0]." ".$clair_info[1]."\n");
			}
		}
		if (starts_with($line, "##clair3_version="))
		{
			$clair_info[1] = trim(substr($line, 17));
			if ($clair_info[0]!="" && $clair_info[1]!="") //write out clair info when complete - no matter which header is first
			{
				fwrite($handle_out, "##SOURCE=".$clair_info[0]." ".$clair_info[1]."\n");
			}
		}
		if (starts_with($line, "##source=DeepVariant"))
		{
			fwrite($handle_out, "##SOURCE=".trim(substr($line,9))."\n");
		}
		if (starts_with($line, "##GLnexusConfigName=")) //after gVCF merging of DeepVariant we need to reconstruct the DeepVariant version (GLnexus removes it...)
		{
			fwrite($handle_out, "##SOURCE=DeepVariant ".get_path("container_deepvariant")."\n");
		}
		if (starts_with($line, "##source=strelka2"))
		{
			fwrite($handle_out, "##SOURCE=".trim(substr($line,9))."\n");
		}
		if (starts_with($line, "##source=VarScan2"))
		{
			fwrite($handle_out, "##SOURCE=".trim(substr($line,9))."\n");
		}
		
		//filters
		if (starts_with($line, "##FILTER=<ID="))
		{
			$parts = explode(",Description=\"", substr(trim($line), 13, -2));
			fwrite($handle_out, "##FILTER=".$parts[0]."=".$parts[1]."\n");
		}
		
		//samples
		if (starts_with($line, "##SAMPLE="))
		{
			$line = trim($line);
			list($name) = explode(",", substr($line, 13, -1));

			if ($genotype_mode=="single")
			{
				// replace sample header with NGSD enry:
				if (db_is_enabled("NGSD") && !$test) $line = gsvar_sample_header($name, array("DiseaseStatus"=>"Affected"), "##", "");
				if ($column_desc[0][0]!="genotype")
				{
					trigger_error("Several sample header lines in 'single' mode!", E_USER_ERROR);
				}
				$column_desc[0][0] = $name;
			}
			else if ($genotype_mode=="multi")
			{
				$multi_cols[] = $name;
				if ($longread)
				{
					// add 2 columns per sample (genotype + phasing info)
					array_splice($column_desc, (2*count($multi_cols))-2, 0, array(array($name, "genotype of sample $name")));
					array_splice($column_desc, (2*count($multi_cols))-1, 0, array(array($name."_phased", "phasing information of sample $name")));
				}
				else
				{
					array_splice($column_desc, count($multi_cols)-1, 0, array(array($name, "genotype of sample $name")));
				}	
			}
			fwrite($handle_out, $line."\n");
		}
		
		//analysis type
		if (starts_with($line, "##ANALYSISTYPE="))
		{
			fwrite($handle_out, trim($line)."\n");
		}
		
		//pipeline
		if (starts_with($line, "##PIPELINE="))
		{
			fwrite($handle_out, trim($line)."\n");
		}
		
		//get annotation indices in CSQ field from VEP
		if (starts_with($line, "##INFO=<ID=CSQ,"))
		{
			$cols = explode("|", substr($line, 0, -2));
			$i_consequence = index_of($cols, "Consequence", "CSQ");
			$i_feature = index_of($cols, "Feature", "CSQ");
			$i_featuretype = index_of($cols, "Feature_type", "CSQ");
			$i_biotype = index_of($cols, "BIOTYPE", "CSQ");
			$i_domains = index_of($cols, "DOMAINS", "CSQ");
			$i_pubmed = index_of($cols, "PUBMED", "CSQ"); 
		}

		//get annotation indices in CSQ field from VcfAnnotateConsequence (also used for CSQ_REFSEQ)
		if (starts_with($line, "##INFO=<ID=CSQ2,"))
		{
			$cols = explode("|", substr($line, 0, -2));
			$i_vac_consequence = index_of($cols, "Consequence", "CSQ2");
			$i_vac_impact = index_of($cols, "IMPACT", "CSQ2");
			$i_vac_symbol = index_of($cols, "SYMBOL", "CSQ2");
			$i_vac_hgnc_id = index_of($cols, "HGNC_ID", "CSQ2");
			$i_vac_feature = index_of($cols, "Feature", "CSQ2");
			$i_vac_exon = index_of($cols, "EXON", "CSQ2");
			$i_vac_intron = index_of($cols, "INTRON", "CSQ2");
			$i_vac_hgvsc = index_of($cols, "HGVSc", "CSQ2");
			$i_vac_hgvsp = index_of($cols, "HGVSp", "CSQ2");		
		}
		//determine if RefSeq annotation is present
		if (starts_with($line, "##INFO=<ID=CSQ_REFSEQ,"))
		{
			$annotate_refseq_consequences = true;
			$column_desc[] = ["coding_and_splicing_refseq", "Variant consequence based on RefSeq transcripts (Gene, ENST number, type, impact, exon/intron number, HGVS.c, HGVS.p)."];
		}
		
		//Targeted info (Dragen targeted caller was used) - is written into filter colum, change header accordingly
		if (starts_with($line, "##INFO=<ID=TARGETED,"))
		{
			fwrite($handle_out, "##FILTER=targeted=Variant called by targeted caller of Illumina DRAGEN\n");	
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
			$i_cosmic_cmc_gene_name = index_of($cols, "GENE_NAME", "COSMIC_CMC");
			$i_cosmic_cmc_mut_id = index_of($cols, "GENOMIC_MUTATION_ID", "COSMIC_CMC");
			$i_cosmic_cmc_disease = index_of($cols, "DISEASE", "COSMIC_CMC");
			$i_cosmic_cmc_dnds_disease = index_of($cols, "DNDS_DISEASE_QVAL_SIG", "COSMIC_CMC");
			$i_cosmic_cmc_mut_sign_tier = index_of($cols, "MUTATION_SIGNIFICANCE_TIER", "COSMIC_CMC");
		}
		
		//Cancerhotspots.org header line
		if( starts_with($line, "##INFO=<ID=CANCERHOTSPOTS,") )
		{
			$skip_cancerhotspots = false;
			$cols = explode("|", substr($line, 0,-2) );
			$cols[0] = "GENE_SYMBOL";			
			$i_cancerhotspots_gene_symbol = index_of($cols, "GENE_SYMBOL", "CANCERHOTSPOTS");
			$i_cancerhotspots_transcript_id = index_of($cols, "ENSEMBL_TRANSCRIPT", "CANCERHOTSPOTS");
			$i_cancerhotspots_aa_pos = index_of($cols, "AA_POS", "CANCERHOTSPOTS");
			$i_cancerhotspots_aa_ref = index_of($cols, "AA_REF", "CANCERHOTSPOTS");
			$i_cancerhotspots_aa_alt = index_of($cols, "AA_ALT", "CANCERHOTSPOTS");
			$i_cancerhotspots_total_count = index_of($cols, "TOTAL_COUNT", "CANCERHOTSPOTS");
			$i_cancerhotspots_alt_count = index_of($cols, "ALT_COUNT", "CANCERHOTSPOTS");
		}

		//Annotation of short-read variants overlap
		if (starts_with($line, "##INFO=<ID=IN_SHORTREAD_SAMPLE"))
		{
			$skip_short_read_overlap_annotation = false;
			//extract process sample name 
			$sr_vcf_file_name = basename2(explode("'", $line)[1]);
			$sr_ps_name = explode("_", $sr_vcf_file_name)[0]."_".explode("_", $sr_vcf_file_name)[1];
		}
		
		//check VCF header
		if (starts_with($line, "#CHROM\t"))
		{
			if ($genotype_mode=="multi")
			{
				//determine multi-sample format according to VCF header
				$cols = explode("\t", trim($line, " \n\r\0\x0B"));
				$n_cols = count($cols);
				$multisample_vcf = ($n_cols > 10);
				if ($multisample_vcf)
				{
					if($n_cols-9 != count($multi_cols))
					{
						trigger_error("VCF column count doesn't match sample count in header!", E_USER_ERROR);
					}
					//verify that VCF header entries match sample entries in the comment section
					$vcf_sample_names =  array_slice($cols, 9);
					if ($vcf_sample_names != $multi_cols) trigger_error("VCF header entries doesn't match sample entries in the comment section!", E_USER_ERROR);
				}

			}
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
	if ($multisample_vcf)
	{
		//each sample has its own FORMAT column
		list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format) = $cols;
		$sample_format_values = array_slice($cols, 9, count($multi_cols));
	}
	else
	{
		list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $sample) = $cols;
	}
	$tag = "{$chr}:{$pos} {$ref}>{$alt}";
	if ($filter=="" || $filter=="." || $filter=="PASS")
	{
		$filter = array();
	}
	else
	{
		$filter = explode(";", $filter);
	}

	//parse variant data from VCF
	if(chr_check($chr, 22, false) === FALSE) continue; //skip bad chromosomes
	$start = $pos;
	$end = $pos;
	$ref = strtoupper($ref);
	$alt = strtoupper($alt);
	if(strlen($ref)>1 || strlen($alt)>1) //correct indels
	{
		list($start, $end, $ref, $alt) = correct_indel($start, $ref, $alt);
	}
	
	//skip too long variants (unique constraint in NGSD fails otherwise)
	if (strlen($ref)>500 || strlen($alt)>500) continue;
	
	//parse info from VCF
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
		
	//special handling for DRAGEN calling
	if (isset($info["TARGETED"]))
	{
		$filter[] = "targeted";
	}
	if ($chr!="chrMT" && isset($info["MOSAIC"]))
	{
		$filter[] = "mosaic";
	}
	
	//convert genotype information to TSV format
	if(!$multisample_vcf)
	{
		$sample = array_combine(explode(":", $format), explode(":", $sample));
	
		//rename Dragen format values for targeted calls.
		if (! isset($sample["DP"]) && isset($sample["JDP"])) $sample["DP"] = $sample["JDP"];
		if (! isset($sample["AF"]) && isset($sample["JAF"])) $sample["AF"] = $sample["JAF"];
		if (! isset($sample["AD"]) && isset($sample["JAD"])) $sample["AD"] = $sample["JAD"];
		if (! isset($sample["PL"]) && isset($sample["JPL"])) $sample["PL"] = $sample["JPL"];
	}
	
	if ($genotype_mode=="multi")
	{
		if($multisample_vcf)
		{
			$sample = array();
			//extract format info and arrange format values based on key
			foreach ($sample_format_values as $format_values) 
			{
				$tmp = array_combine(explode(":", $format), explode(":", $format_values));
				foreach ($tmp as $key => $value) 
				{
					$sample[$key][] = $value;	
				}
			}
			
			//rename Dragen format values for targeted calls.
			if (! isset($sample["DP"]) && isset($sample["JDP"])) $sample["DP"] = $sample["JDP"];
			if (! isset($sample["AF"]) && isset($sample["JAF"])) $sample["AF"] = $sample["JAF"];
			if (! isset($sample["AD"]) && isset($sample["JAD"])) $sample["AD"] = $sample["JAD"];
			if (! isset($sample["PL"]) && isset($sample["JPL"])) $sample["PL"] = $sample["JPL"];
			
			//determine human-readable genotype
			$genotypes = array();
			for ($i=0; $i < count($sample["GT"]); $i++) 
			{ 
				$gt = vcfgeno2human($sample["GT"][$i]);
				if ($sample["DP"][$i]<3) $gt = "n/a";
				$genotypes[] = $gt;
				if ($longread)
				{
					$phasing_info = "";
					if (strpos($sample["GT"][$i], "|") !== false) 
					{
						$phasing_info = $sample["GT"][$i]." (".$sample["PS"][$i].")";
					}
				}
				$genotypes[] = $phasing_info;
			}

			//combine certain values
			$genotype = "\t".implode("\t", $genotypes);
			$sample["DP"] = implode(",", $sample["DP"]);
			$sample["AF"] = implode(",", $sample["AF"]);
		}
		else
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
	}
	else if ($genotype_mode=="single")
	{
		if (!isset($sample["GT"])) 
		{
			trigger_error("VCF sample column does not contain GT value!", E_USER_ERROR);
		}

		// get phasing info for long-reads
		if ($longread)
		{
			$phasing_info = "";
			if (strpos($sample["GT"], "|") !== false) 
			{
				$phasing_info = $sample["GT"]." (".$sample["PS"].")";
			}
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
	
	//QD
	if (isset($sample["DP"])) //freebayes
	{
		if (is_numeric($sample["DP"]) && floatval($sample["DP"]) != 0)
		{
			$quality[] = "QD=".number_format($qual/$sample["DP"], 2);
		}
	}
	if (isset($sample["QD"])) //Dragen
	{
		$quality[] = "QD=".$sample["QD"];
	}
	
	//AF
	if (isset($sample["AO"]) && isset($sample["DP"])) //freebayes
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
	if (isset($sample["AF"])) //Dragen
	{
		$quality[] = "AF=".$sample["AF"];
	}
	
	//MQM
	if (isset($info["MQM"])) //freebayes
	{
		$quality[] = "MQM=".intval($info["MQM"]);
	}
	if (isset($info["MQ"])) //Dragen
	{
		$quality[] = "MQM=".intval($info["MQ"]);
	}
	
	//other quality fields specific to freebayes
	if (isset($info["SAP"])) 
	{
		$quality[] = "SAP=".intval($info["SAP"]);
	}
	if (isset($info["ABP"])) 
	{
		$quality[] = "ABP=".intval($info["ABP"]);
	}
	if (isset($info["SAR"]))
	{
		$quality[] = "SAR=".intval($info["SAR"]);
	}
	if (isset($info["SAF"]))
	{
		$quality[] = "SAF=".intval($info["SAF"]);
	}
	
	$phylop = array();
	if (isset($info["PHYLOP"])) 
	{
		$phylop[] = $info["PHYLOP"];
	}

	$revel = array();
	if (isset($info["REVEL"])) 
	{
		$revel = explode("&", $info["REVEL"]);
	}
	
	$alphamissense = [];
	if (isset($info["AM_MAIN"])) 
	{
		$alphamissense = explode("&", $info["AM_MAIN"]);
	}
	if (isset($info["AM_ISO"])) 
	{
		$alphamissense = explode("&", $info["AM_ISO"]);
	}
	
	//variant details
	$dbsnp = array();
	$genes = array();
	$variant_details = array();
	$coding_and_splicing_details = array();
	$coding_and_splicing_refseq = array();
	$af_gnomad_genome = array();
	$af_gnomad_afr = array();
	$af_gnomad_amr = array();
	$af_gnomad_eas = array();
	$af_gnomad_nfe = array();
	$af_gnomad_sas = array();
	$hom_gnomad = array();
	$hemi_gnomad = array();
	$wt_gnomad = array();
	$het_gnomad = array();
	$clinvar = array();
	$hgmd = array();
	$maxentscan = array();
	$regulatory = array();
	$pubmed = array();
	$custom_column_data = [];
	
	//variant details based on Ensembl (up/down-stream)
	$variant_details_updown = array();
	$genes_updown = array();
	$coding_and_splicing_details_updown = array(); 
	if (isset($info["CSQ"]) && isset($info["CSQ2"]))
	{
		//VEP - used for regulatory features, PubMed, Domains
		$vep = []; //transcript name without version > domain
		foreach(explode(",", $info["CSQ"]) as $entry)
		{			
			$parts = explode("|", $entry);
			
			//######################### general information (not transcript-specific) #########################
			
			//PubMed ids
			if ($i_pubmed!==FALSE)
			{
				$pubmed = array_merge($pubmed, explode("&", $parts[$i_pubmed]));
			}
			
			//######################### transcript-specific information #########################
			$feature_type = trim($parts[$i_featuretype]);
			if ($feature_type=="Transcript")
			{
				$transcript_id = trim($parts[$i_feature]);	

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
				if ($domain != "")
				{
					$domain_description = "";
					
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
						$domain_description .= $pfam_description[$domain];
					}

					// throw error if Pfam id is neither found in replacement data nor in description data
					if ($domain_description == "")
					{
						$missing_domains[$domain] = true;
					}

					// combine decription and id
					$domain = "$domain [$domain_description]";
				}
				
				$transcript_id_no_ver = explode('.', $transcript_id)[0];
				
				$vep[$transcript_id_no_ver] = $domain;			
			}
			else if ($feature_type=="RegulatoryFeature")
			{
				$regulatory[] = $parts[$i_consequence].":".$parts[$i_biotype];
			}
			else if ($feature_type=="MotifFeature")
			{
				$regulatory[] = $parts[$i_consequence];
			}
			else if ($feature_type!="") //feature type is empty for intergenic variants
			{				
				trigger_error("Unknown VEP feature type '{$feature_type}' for variant {$chr}:{$pos} {$ref}>{$alt}!", E_USER_ERROR);
			}
		}
		
		//VcfAnnotateConsequence (Ensembl)
		foreach(explode(",", $info["CSQ2"]) as $entry)
		{			
			$entry = trim($entry);
			if ($entry=="") continue;
			
			$parts = explode("|", $entry);
			
			$transcript_id = trim($parts[$i_vac_feature]);
			
			$consequence = $parts[$i_vac_consequence];
			$consequence = strtr($consequence, ["&NMD_transcript_variant"=>"", "splice_donor_variant&intron_variant"=>"splice_donor_variant", "splice_acceptor_variant&intron_variant"=>"splice_acceptor_variant"]);
			
			//extract variant type
			$variant_type = strtr($consequence, array("_variant"=>"", "_prime_"=>"'"));
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
			$gene = trim($parts[$i_vac_symbol]);
			if ($gene!="")
			{
				$hgnc_id = $parts[$i_vac_hgnc_id];
				$hgnc_id = trim(strtr($hgnc_id, array("HGNC:"=>"")));
				if (isset($hgnc[$hgnc_id]))
				{
					$hgnc_gene = $hgnc[$hgnc_id];
					if ($gene!=$hgnc_gene)
					{
						$gene = $hgnc_gene;
					}
				}
				else if ($hgnc_id!="")
				{
					@$hgnc_messages["ID '$hgnc_id' not valid or withdrawn for gene '$gene'"] += 1;
				}
				
				if (!$is_updown)
				{
					 $genes[] = $gene;
				}
				else
				{
					$genes_updown[] = $gene;
				}
			}
			
			//exon
			$exon = trim($parts[$i_vac_exon]);
			if ($exon!="") $exon = "exon".$exon;
			$intron = trim($parts[$i_vac_intron]);
			if ($intron!="") $intron = "intron".$intron;
			
			//hgvs
			$hgvs_c = trim($parts[$i_vac_hgvsc]);
			$hgvs_p = trim($parts[$i_vac_hgvsp]);
			$hgvs_p = str_replace("%3D", "=", $hgvs_p);
			
			//domain from VEP
			$domain = "";
			$transcript_id_no_ver = explode('.', $transcript_id)[0];
			if (isset($vep[$transcript_id_no_ver]))
			{
				$domain = $vep[$transcript_id_no_ver];
			}

			//add transcript information
			$transcript_entry = "{$gene}:{$transcript_id}:".$consequence.":".$parts[$i_vac_impact].":{$exon}{$intron}:{$hgvs_c}:{$hgvs_p}:{$domain}";
			if (!$is_updown)
			{
				$coding_and_splicing_details[] = $transcript_entry;
			}
			else
			{
				$coding_and_splicing_details_updown[] = $transcript_entry;
			}
		}	
	}

	if (isset($info["CSQ_REFSEQ"]))
	{
		//VcfAnnotateConsequence (RefSeq)
		foreach(explode(",", $info["CSQ_REFSEQ"]) as $entry)
		{
			$entry = trim($entry);
			if ($entry=="") continue;
			
			$parts = explode("|", $entry);
			
			$transcript_id = trim($parts[$i_vac_feature]);
			
			$consequence = $parts[$i_vac_consequence];
			$consequence = strtr($consequence, ["&NMD_transcript_variant"=>"", "splice_donor_variant&intron_variant"=>"splice_donor_variant", "splice_acceptor_variant&intron_variant"=>"splice_acceptor_variant"]);
			
			//determine gene name (update if neccessary)
			$gene = trim($parts[$i_vac_symbol]);
			if ($gene!="")
			{
				$hgnc_id = $parts[$i_vac_hgnc_id];
				$hgnc_id = trim(strtr($hgnc_id, array("HGNC:"=>"")));
				if (isset($hgnc[$hgnc_id]))
				{
					$hgnc_gene = $hgnc[$hgnc_id];
					if ($gene!=$hgnc_gene)
					{
						$gene = $hgnc_gene;
					}
				}
			}
			
			//exon
			$exon = trim($parts[$i_vac_exon]);
			if ($exon!="") $exon = "exon".$exon;
			$intron = trim($parts[$i_vac_intron]);
			if ($intron!="") $intron = "intron".$intron;
			
			//hgvs
			$hgvs_c = trim($parts[$i_vac_hgvsc]);
			$hgvs_p = trim($parts[$i_vac_hgvsp]);
			$hgvs_p = str_replace("%3D", "=", $hgvs_p);
			
			//add transcript information
			$coding_and_splicing_refseq[] = "{$gene}:{$transcript_id}:".$consequence.":".$parts[$i_vac_impact].":{$exon}{$intron}:{$hgvs_c}:{$hgvs_p}";
		}	
	}
	

	//dbSNP
	$dbsnp = [];
	if (isset($info["RS"]))
	{
		$rs = trim($info["RS"]);
		if ($rs!="") $dbsnp[] = $rs;
	}
	
	//MaxEntScan
	$mes_by_trans = [];
	if (isset($info["MES"])) //parse MES scores for native splice site
	{
		foreach(explode("|", $info["MES"]) as $mes_entry)
		{
			$parts = explode("&", trim($mes_entry));
			if (count($parts)!=3) continue;
			
			list($mes_ref, $mes_alt, $mes_trans) = $parts;
			$mes_by_trans[$mes_trans] = [$mes_ref, $mes_alt, "", "", "", ""];
		}
	}
	if (isset($info["MES_SWA"])) //parse MES scores generated by SWA
	{
		foreach(explode("|", $info["MES_SWA"]) as $mes_entry)
		{
			$parts = explode("&", trim($mes_entry));
			if (count($parts)!=7) continue;
			
			list($d_mes_ref, $d_mes_alt, $d_mes_comp, $a_mes_ref, $a_mes_alt, $a_mes_comp, $mes_trans) = $parts;
			
			$mes_ref = "";
			$mes_alt = "";
			if (isset($mes_by_trans[$mes_trans]))
			{
				list($mes_ref, $mes_alt) = $mes_by_trans[$mes_trans];
			}
			$mes_by_trans[$mes_trans] = [$mes_ref, $mes_alt, $a_mes_ref, $a_mes_alt, $d_mes_ref, $d_mes_alt];
		}
	}
	foreach($mes_by_trans as $mes_trans => list($mes_ref, $mes_alt, $a_mes_ref, $a_mes_alt, $d_mes_ref, $d_mes_alt))
	{
		$parts = ["", "", ""];
		if($mes_ref!="" && $mes_alt!="")
		{
			$parts[0] = $mes_ref.">".$mes_alt;
		}
		if($a_mes_ref!="" && $a_mes_alt!="")
		{
			$parts[1] = $a_mes_ref.">".$a_mes_alt;
		}
		if($d_mes_ref!="" && $d_mes_alt!="")
		{
			$parts[2] = $d_mes_ref.">".$d_mes_alt;
		}
		$maxentscan[] = implode("/", $parts);
	}

	//gnomAD AF
	if (isset($info["gnomADg_AF"])) $af_gnomad_genome = explode("&", $info["gnomADg_AF"]);
	if (isset($info["gnomADm_AF_hom"])) $af_gnomad_genome = explode("&", $info["gnomADm_AF_hom"]);
	
	//gnomAD hom/hemi
	if (isset($info["gnomADg_Hom"])) $hom_gnomad = explode("&", $info["gnomADg_Hom"]);
	if (isset($info["gnomADg_Hemi"])) $hemi_gnomad = explode("&", $info["gnomADg_Hemi"]);
	if (isset($info["gnomADg_Het"])) $het_gnomad = explode("&", $info["gnomADg_Het"]);
	if (isset($info["gnomADg_Wt"])) $wt_gnomad = explode("&", $info["gnomADg_Wt"]);

	//gnomAD sub-populations
	if (isset($info["gnomADg_AFR_AF"])) $af_gnomad_afr = explode("&", $info["gnomADg_AFR_AF"]);
	if (isset($info["gnomADg_AMR_AF"])) $af_gnomad_amr = explode("&", $info["gnomADg_AMR_AF"]);
	if (isset($info["gnomADg_EAS_AF"])) $af_gnomad_eas = explode("&", $info["gnomADg_EAS_AF"]);
	if (isset($info["gnomADg_NFE_AF"])) $af_gnomad_nfe = explode("&", $info["gnomADg_NFE_AF"]);
	if (isset($info["gnomADg_SAS_AF"])) $af_gnomad_sas = explode("&", $info["gnomADg_SAS_AF"]);

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
		
		$hgmd[] = trim($hgmd_id." [CLASS=".$hgmd_class." MUT=".$hgmd_mut." PHEN=".strtr($hgmd_phen, "_", " ")." GENE=".$hgmd_gene."];");
	}

	//AFs
	$dbsnp = implode(",", collapse($tag, "dbSNP", $dbsnp, "unique"));	
	$gnomad = collapse($tag, "gnomAD genome", $af_gnomad_genome, "max", 5);
	$gnomad_hom_hemi = collapse($tag, "gnomAD Hom", $hom_gnomad, "max").",".collapse($tag, "gnomAD Hemi", $hemi_gnomad, "max");
	if ($gnomad_hom_hemi==",") $gnomad_hom_hemi = "";
	$gnomad_sub = collapse($tag, "gnomAD AFR", $af_gnomad_afr, "max", 5).",".collapse($tag, "gnomAD AMR", $af_gnomad_amr, "max", 5).",".collapse($tag, "gnomAD EAS", $af_gnomad_eas, "max", 5).",".collapse($tag, "gnomAD NFE", $af_gnomad_nfe, "max", 5).",".collapse($tag, "gnomAD SAS", $af_gnomad_sas, "max", 5);
	if (str_replace(",", "", $gnomad_sub)=="") $gnomad_sub = "";
	$gnomad_het = collapse($tag, "gnomAD Het", $het_gnomad, "max");
	$gnomad_wt = collapse($tag, "gnomAD Wt", $wt_gnomad, "max");

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
			
			if (isset($info["NGSD_SOM_TO_C"]))
			{
				$ngsd_som_counts_to = intval(trim($info["NGSD_SOM_TO_C"]));
			}
			else
			{
				$ngsd_som_counts_to = "0";
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
		if (isset($info["NGSD_HAF"]) || $gnomad >= 0.05)
		{
			$ngsd_hom = "n/a (AF>5%)";
			$ngsd_het = "n/a (AF>5%)";
			$ngsd_mosaic = "n/a (AF>5%)";
			$ngsd_group = "n/a (AF>5%)";
		}
		elseif(isset($info["NGSD_COUNTS"]))
		{
			$ngsd_counts = explode(",", trim($info["NGSD_COUNTS"]));
			$ngsd_hom = $ngsd_counts[0];
			$ngsd_het = $ngsd_counts[1];
			$ngsd_mosaic = (count($ngsd_counts)<3 ? "n/a" : $ngsd_counts[2]);
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
			$ngsd_mosaic = "0";
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
		$tmp = [];
		$spliceai_info = trim($info["SpliceAI"]);
		$spliceai_values = array();

		$entries = explode(",", strtr($spliceai_info, "&", ",")); //both & and , are used as separator, depending on the sources of the SpliceAI annotation (pre-calcualted or calculated on the fly)
		foreach($entries as $entry)
		{
			$delta_scores = explode("|", $entry);
			if(count($delta_scores) == 10)
			{
				$tmp[] = implode("|", array_slice($delta_scores, 1));
			}
			else
			{
				trigger_error("Wrong SpliceAI annotation in line: ${line} in SpliceAI annotation: ${spliceai_info}! Delimiter for several genes must be ','.", E_USER_WARNING);
			}
		}
		$tmp = array_unique($tmp);
		sort($tmp);

		if(count($tmp)>0)
		{
			$spliceai = implode(",", $tmp);
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
	$revel = empty($revel) ? "" : collapse($tag, "REVEL", $revel, "max", 3);
	$alphamissense = empty($alphamissense) ? "" : collapse($tag, "AlphaMissense", $alphamissense, "max", 2);
	
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
	
	//custom columns
	foreach($custom_columns as $key => $tmp)
	{
		list($vcf, $col, $desc) = explode(";", $tmp, 3);
		$custom_column_data[$key] = isset($info["CUSTOM_".$col]) ? $info["CUSTOM_".$col] : "";
	}
	
	//skip common MODIFIER variants in WGS mode
	if ($wgs && skip_in_wgs_mode($chr, $coding_and_splicing_details, $gnomad, $clinvar, $hgmd, $ngsd_clas))
	{
		++$c_skipped_wgs;
		continue;
	}

	//sr overlap annotation
	if (!$skip_short_read_overlap_annotation)
	{
		$in_short_read_sample = "";
		if (isset($info["IN_SHORTREAD_SAMPLE"])) $in_short_read_sample = "1";
	}
	
	//write data
	++$c_written;
	$genes = array_unique($genes);
	fwrite($handle_out, "$chr\t$start\t$end\t$ref\t{$alt}{$genotype}");
	if($longread && ($genotype_mode != "multi")) //for multi-sample the phasing info is added before
	{
		fwrite($handle_out,"\t".$phasing_info);
	}
	fwrite($handle_out,"\t".implode(";", $filter)."\t".implode(";", $quality)."\t".implode(",", $genes)."\t$variant_details\t$coding_and_splicing_details\t$regulatory\t$omim\t$clinvar\t$hgmd\t$repeatmasker\t$dbsnp\t$gnomad\t$gnomad_sub\t$gnomad_hom_hemi\t$gnomad_het\t$gnomad_wt\t$phylop\t$cadd\t$revel\t$alphamissense\t$maxentscan\t$spliceai\t$pubmed");
	if ($annotate_refseq_consequences)
	{
		fwrite($handle_out, "\t".implode(",", $coding_and_splicing_refseq));
	}
	if (!$skip_ngsd_som)
	{
		fwrite($handle_out, "\t$ngsd_som_counts\t$ngsd_som_projects\t$ngsd_som_counts_to\t$ngsd_som_vicc\t$ngsd_som_vicc_comment");
	}
	if (!$skip_ngsd)
	{
		fwrite($handle_out, "\t$ngsd_hom\t$ngsd_het\t$ngsd_mosaic\t$ngsd_group\t$ngsd_clas\t$ngsd_clas_com\t$ngsd_val\t$ngsd_com\t$ngsd_gene_info");
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

	if (!$skip_short_read_overlap_annotation)
	{
		fwrite($handle_out, "\t".$in_short_read_sample);
	}
	
	foreach($custom_columns as $key => $tmp)
	{
		fwrite($handle_out, "\t".$custom_column_data[$key]);
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

gzclose($handle);
fclose($handle_out);

if (count($missing_domains)>0)
{
	trigger_error("No description found for the folling domains: ".implode(", ", array_keys($missing_domains)), E_USER_WARNING);
}

//print debug output
print "Variants written: {$c_written}\n";
if ($wgs)
{
	print "Variants skipped because WGS mode is enabled: {$c_skipped_wgs}\n";
}


?>
