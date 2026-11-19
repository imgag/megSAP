<?php 
/** 
	@page an_aidiva_hpo
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("an_aidiva_hpo", "Annotates a VCF file with aiDIVA score (independent of pheontype).");
$parser->addInfile("in",  "Input file in VCF or VCF.GZ format that is annotated with 'an_aidiva'.", false);
$parser->addOutfile("out", "Output file in VCF format.", false);
$parser->addOutfile("hpo_list", "Comma-speparated list of HPO terms: HP:XXXXXX,HP:XXXXXX", false);
$parser->addInteger("threads", "Number of threads to use.", true, 1); //TODO use
extract($parser->parse($argv));

//TODO convert to function
//- AIDIVA_SCORE via flag: if true, check that annotated with an_aidiva.php



//convert VCF to intput TSV
$tsv_in = $parser->tempFile(".tsv");
$ho = fopen($tsv_in, "w");
fputs($ho, "#CHROM\tPOS\tREF\tALT\tFILTER\tIS_INDEL\tAllele\tConsequence\tIMPACT\tSYMBOL\tALPHA_MISSENSE_SCORE\tCADD_PHRED\tREVEL\tphyloP_vertebrate\toe_lof\tHIGH_IMPACT\tGT\tMAX_AF\tSpliceAI\n"); //TODO GT suffix with PS name: GT_DX110151_01

$h = gzopen2($in, "r");
while(!gzeof($h))
{
	$line = nl_trim(gzgets($h));
	if ($line=="" || $line[0]=="#") continue;
	
	list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $sample) = explode("\t", $line);
	$info = explode(";", $info);
	
	//use only low_conf_region filter
	if (contains($filter, "low_conf_region")) $filter = "low_conf_region";
	else $filter = "";
	
	//get GT
	$format = explode(":", $format);
	if ($format[0]!="GT") continue;
	$gt = explode(":", $sample)[0];
	
	//normalize GT
	$gt = strtr($gt, ["."=>"0", "|"=>"/"]);
	if ($gt=="0/0") continue;
	
	//parse INFO
	$allele = "";
	$consequence = "";
	$impact = "";
	$symbol = "";
	$high_impact = "0";
	$alphamissense = "";
	$cadd = "";
	$revel = "";
	$phylop = "";
	$oe_lof = "";
	$maxaf = 0;
	$spliceai = "";
	foreach($info as $key)
	{
		if (contains($key, "=")) list($key, $value) = explode("=", $key, 2);
		
		if ($key=="CADD_SNV") $cadd = $value;
		if ($key=="CADD_INDEL") $cadd = $value;
		
		if ($key=="gnomADg_AFR_AF" || $key=="gnomADg_AMR_AF" || $key=="gnomADg_EAS_AF" || $key=="gnomADg_NFE_AF" || $key=="gnomADg_SAS_AF") $maxaf = max($maxaf, $value);
		
		if ($key=="REVEL") $revel = $value;
		
		if ($key=="PHYLOP") $phylop = $value;
		
		//SpliceAI=G|PTCH1_ENST00000331920|0.01|0.00|0.00|0.00|-72|-14|255|-128,G|PTCH1_ENST00000437951|0.01|0.00|0.00|0.00|-72|-14|255|-128
		if ($key=="SpliceAI")
		{
			$spliceai = 0;
			foreach(explode(",", $value) as $entry)
			{
				$parts = explode("|", $entry);
				for ($i=2; $i<6; ++$i)
				{
					$spliceai = max($spliceai, $parts[$i]);
				}
			}
			
		}
		
		//NGSD_GENE_INFO=USHBP1%20(oe_syn%3D0.89%20oe_mis%3D1.02%20oe_lof%3D0.88%20pli%3D0.00)&BABAM1%20(oe_syn%3D1.02%20oe_mis%3D0.91%20oe_lof%3D0.41%20pli%3D0.34)
		if ($key=="NGSD_GENE_INFO")
		{
			$oe_lof = 0;
			foreach(explode("&", $value) as $entry)
			{
				$parts = explode(" ", strtr($entry, ["%20"=>" ", "%3D"=>"="]));
				foreach($parts as $part)
				{
					if (starts_with($part, "oe_lof="))
					{
						$oe_lof = max($oe_lof, explode("=", $part)[1]);
					}
				}
			}
		}
		
		if ($key=="AM_MAIN")
		{
			$alphamissense = $value;
		}
		
		//CSQ=A|intron_variant|MODIFIER|CHEK2|HGNC:16627|ENST00000928683.1|Transcript||1/14|c.-6-99dup|,A|intron_variant|MODIFIER|CHEK2|HGNC:16627|ENST00000928687.1|Transcript||1/10|c.-6-99dup|,A|upstream_gene_variant|MODIFIER|CHEK2|HGNC:16627|ENST00000402731.6|Transcript||||,A|intron_variant|MODIFIER|CHEK2|HGNC:16627|ENST00000928694.1|Transcript||3/16|c.-6-99dup|,A|intron_variant|MODIFIER|CHEK2|HGNC:16627|ENST00000928693.1|Transcript||1/14|c.-6-99dup|,A|intron_variant|MODIFIER|CHEK2|HGNC:16627|ENST00000425190.7|Transcript||1/13|c.-345+6943dup|,A|intron_variant|MODIFIER|CHEK2|HGNC:16627|ENST00000928690.1|Transcript||2/15|c.-6-99dup|,A|intron_variant|MODIFIER|CHEK2|HGNC:16627|ENST00000928688.1|Transcript||1/12|c.-6-99dup|,A|intron_variant|MODIFIER|CHEK2|HGNC:16627|ENST00000928689.1|Transcript||1/9|c.-6-99dup|,A|intron_variant|MODIFIER|CHEK2|HGNC:16627|ENST00000404276.6|Transcript||1/14|c.-6-99dup|,A|intron_variant|MODIFIER|CHEK2|HGNC:16627|ENST00000928686.1|Transcript||1/8|c.-6-99dup|,A|intron_variant|MODIFIER|CHEK2|HGNC:16627|ENST00000348295.7|Transcript||1/13|c.-6-99dup|,A|intron_variant|MODIFIER|CHEK2|HGNC:16627|ENST00000928684.1|Transcript||1/13|c.-6-99dup|,A|5_prime_UTR_variant|MODIFIER|CHEK2|HGNC:16627|ENST00000899937.1|Transcript|1/14||c.-105dup|,A|intron_variant|MODIFIER|CHEK2|HGNC:16627|ENST00000928691.1|Transcript||2/15|c.-6-99dup|,A|intron_variant|MODIFIER|CHEK2|HGNC:16627|ENST00000954618.1|Transcript||2/15|c.-6-99dup|,A|intron_variant|MODIFIER|CHEK2|HGNC:16627|ENST00000928695.1|Transcript||1/13|c.-6-99dup|,A|5_prime_UTR_variant|MODIFIER|CHEK2|HGNC:16627|ENST00000928692.1|Transcript|2/14||c.-105dup|,A|intron_variant|MODIFIER|CHEK2|HGNC:16627|ENST00000928685.1|Transcript||1/14|c.-6-99dup|,A|intron_variant|MODIFIER|CHEK2|HGNC:16627|ENST00000928696.1|Transcript||1/12|c.-6-99dup|,A|intron_variant|MODIFIER|CHEK2|HGNC:16627|ENST00000649563.1|Transcript||1/11|c.-72+6943dup|,A|intron_variant|MODIFIER|CHEK2|HGNC:16627|ENST00000405598.5|Transcript||2/15|c.-6-99dup|,A|intron_variant|MODIFIER|CHEK2|HGNC:16627|ENST00000382580.6|Transcript||1/15|c.-6-99dup|,A|intron_variant|MODIFIER|CHEK2|HGNC:16627|ENST00000650281.1|Transcript||2/15|c.-6-99dup|,A|upstream_gene_variant|MODIFIER|CHEK2|HGNC:16627|ENST00000403642.5|Transcript||||
		if ($key=="CSQ")
		{
			$max_im = "MODIFIER";
			$max = [];
			foreach(explode(",", $value) as $entry)
			{
				$parts = explode("|", $entry);
				if ($parts[2]=="HIGH")
				{
					$max = $parts;
					$max_impact = "HIGH";
				}
				if ($parts[2]=="MODERATE" && $max_im!="HIGH")
				{
					$max = $parts;
					$max_impact = "MODERATE";
				}
				if ($parts[2]=="LOW" && $max_im!="HIGH" && $max_im!="MODERATE")
				{
					$max = $parts;
					$max_impact = "MODERATE";
				}
				if ($parts[2]=="MODIFIER" && $max_im=="MODIFIER")
				{
					$max = $parts;
					$max_impact = "MODERATE";
				}
			}
			list($allele, $consequence, $impact, $symbol) = $max;
			if ($impact=="HIGH") $high_impact = "1";
		}
	}
	
	//is indel?
	$is_indel = strlen($ref)>1 || strlen($alt)>1 ? "1" : "0";
	
	//score only AF<=1%
	if ($maxaf>0.01) continue;
	
	fputs($ho, implode("\t", [$chr, $pos, $ref, $alt, $filter, $is_indel, $allele, $consequence, $impact, $symbol, $alphamissense, $cadd, $revel, $phylop, $oe_lof, $high_impact, $gt, $maxaf, $spliceai])."\n");
}
gzclose($h);
fclose($ho);

//create config file based on template
$model = get_path("data_folder")."/dbs/aiDIVA/aiDIVA_megSAP_model.pkl.gz"; 
$content = file_get_contents(repository_basedir()."/data/misc/aidiva/config_score_hpo.yaml");
$content = str_replace("[MODEL]", $model, $content);
$config = $parser->tempFile("aidiva_config.yaml");
file_put_contents($config, $content);

//run aiDIVA
$out_prefix = $parser->tempFolder()."/out";
$parser->execApptainer("aidiva", "python3", "/opt/aidiva_local/aidiva/run_aidiva-rf.py --in_data $tsv_in --config $config --out_prefix $out_prefix", [$model]);

//get aiDIVA scores form output TSV
$scores = [];
foreach(file("{$out_prefix}_result_aidiva-score.tsv") as $line)
{
	if ($line=="" || $line[0]=="#") continue;
	$parts = explode("\t", $line);
	list($chr, $pos, $ref, $alt) = $parts;
	$score = trim($parts[20]);
	if (!is_numeric($score)) continue;
	
	$scores["{$chr}_{$pos}_{$ref}_{$alt}"] = number_format($score, 3);
}

//annotate aiDIVA scores to VCF
$ho = fopen($out, "w");
$h = gzopen2($in, "r");
while(!gzeof($h))
{
	$line = nl_trim(gzgets($h));
	if ($line=="") continue;
	
	if($line[0]=="#")
	{
		if (!starts_with($line, "##"))
		{
			fputs($ho, "##INFO=<ID=AIDIVA,Number=1,Type=Float,Description=\"aiDIVA pathogenicity score (independent of phenotype).\">\n");
		}
		fputs($ho, $line."\n");
		continue;
	}
	
	$parts = explode("\t", $line);
	list($chr, $pos, $id, $ref, $alt) = $parts;
	$tag = "{$chr}_{$pos}_{$ref}_{$alt}";
	if (isset($scores[$tag]))
	{
		$parts[7] .= ";AIDIVA=".$scores[$tag];
	}
	fputs($ho, implode("\t", $parts)."\n");
}
gzclose($h);
fclose($ho);

//TODO:
//- read pre-calculated HPO data from cache
//- use transcript information from Ensembl BioMart from cache

?>