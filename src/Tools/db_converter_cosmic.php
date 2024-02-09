<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once("/mnt/storage2/users/ahott1a1/main/megSAP/src/Common/all.php");

// parse command line arguments
$parser = new ToolBase("db_converter_cosmic", "Converts COSMIC Cancer Mutation Census (CMC) database file (TSV) to tabix-indexed VCF.GZ.");
// optional
$parser->addInfile("in_cmc",  "Input file in TSV-format with COSMIC CMC scores. [or '-' for stdin]", false, false);
$parser->addInfile("in_genome_vcf",  "Input file in VCF.GZ-format with COSMIC genome mutations.", false, false);
$parser->addInfile("in_non_coding_vcf",  "Input file in VCF.GZ-format with COSMIC non-coding mutations.", false, false);
$parser->addInfile("in_target_screens_vcf",  "Input file in VCF.GZ-format with COSMIC non-coding mutations.", false, false);
$parser->addOutfile("out",  "Output file in vcf.gz-format with CMC data.", false, false);
$parser->addEnum("build", "The genome build to use.", true, ["GRCh37", "GRCh38"], "GRCh38");
extract($parser->parse($argv));

function parse_cmc_file($in_cmc, $build)
{
	if ($in_cmc == "-") $in_cmc = "php://stdin";
	$in_fp = fopen2($in_cmc, "r");
	
	//get indices (header/first line of input file)
	$indices = array();
	$parts = explode( "\t",trim(fgets($in_fp)) );
	for($i=0; $i<count($parts); ++$i)
	{
		$indices[$parts[$i]] = $i;
	}
	//columns to be included in INFO field (without coordinates/reference/alteration fields)
	$ann_cols = array_values(
		array_diff(
			array_keys($indices), 
			["GENOMIC_WT_ALLELE_SEQ","GENOMIC_MUT_ALLELE_SEQ", "Mutation genome position GRCh37", "Mutation genome position GRCh38", "EXAC_AF" , "EXAC_AFR_AF" , "EXAC_AMR_AF" , "EXAC_EAS_AF" , "EXAC_FIN_AF" , "EXAC_NFE_AF" , "EXAC_SAS_AF" , "GNOMAD_EXOMES_AF" , "GNOMAD_EXOMES_AFR_AF" , "GNOMAD_EXOMES_AMR_AF" , "GNOMAD_EXOMES_ASJ_AF" , "GNOMAD_EXOMES_EAS_AF" , "GNOMAD_EXOMES_FIN_AF" , "GNOMAD_EXOMES_NFE_AF" , "GNOMAD_EXOMES_SAS_AF" , "GNOMAD_EXOMES_OTH_AF" , "GNOMAD_GENOMES_AF" , "GNOMAD_GENOMES_AFR_AF" , "GNOMAD_GENOMES_AMR_AF" , "GNOMAD_GENOMES_ASJ_AF" , "GNOMAD_GENOMES_EAS_AF" , "GNOMAD_GENOMES_FIN_AF" , "GNOMAD_GENOMES_NFE_AF" , "GNOMAD_GENOMES_OTH_AF"]
		) 
	);
	
	$cosmic_id_idx = array_search("GENOMIC_MUTATION_ID", $parts);
	
	if ($cosmic_id_idx === false) trigger_error("Column 'GENOMIC_MUTATION_ID' not found in given cmc file.", E_USER_ERROR);
	
	$vars_id_anno = array();
	
	// parse input file 
	while(!feof($in_fp))
	{
		$line = trim( fgets($in_fp) );

		//skip comments and empty lines
		if( starts_with($line, "#") ) continue;
		if( empty($line) ) continue;
		
		
		$parts = explode( "\t" , $line );
		
		if(count($parts) != count($indices) )
		{
			trigger_error("Invalid number of columns (" . count($parts) . "). File header has " . count($indices) . " columns."  , E_USER_ERROR);
		}
		
		//variant has no positions for the given build -> skip it!
		$coords = $parts[ $indices["Mutation genome position $build"] ];
		if ($coords == "") continue;
		
		
		$cosmic_id = $parts[$cosmic_id_idx];
		
		$cmc_data = array();
		foreach($ann_cols as $col_name)
		{
			$temp = $parts[$indices[$col_name]];
			$temp = str_replace(",", ";", $temp); //replace "," by ";", otherwise there is no in-field separator for later processing left.
			$cmc_data[] = vcf_encode_url_string($temp);
		}
		
		$vars_id_anno[$cosmic_id] = "COSMIC_CMC=". implode("|", $cmc_data);
	}

	return [$vars_id_anno, $ann_cols];
}

function write_variants_from_file(&$vars_id_anno, &$out_fp, $cosmic_vcf)
{
	//gather positions and ref/alt allels from the cosmic genome mutation vcf file!
	$cos_vcf_fp = gzopen($cosmic_vcf, "r");
	$count = 0;

	while($cos_vcf_fp && !gzeof($cos_vcf_fp))
	{
		$line = gzgets($cos_vcf_fp);
		
		if (empty($line) || $line[0] == '#') continue;
		
		list($chr, $pos, $cos_id, $ref, $alt, $qual, $filter, $info) = explode("\t", $line);
		
		if (array_key_exists($cos_id, $vars_id_anno))
		{
			if ($count % 500000 == 0 && $count != 0) echo "Found var_ids in comsic vcf file: $count\n";
			
			//convert chromsome names
			$chr = "chr" . $chr;
			$chr = str_replace("chr23", "chrX", $chr);
			$chr = str_replace("chr24", "chrY", $chr);
			$chr = str_replace("chr25", "chrMT", $chr);
			
			fwrite($out_fp, "$chr\t$pos\t.\t$ref\t$alt\t.\t.\t".$vars_id_anno[$cos_id]."\n");
			$count++;
			//unset variants that were already written:
			unset($vars_id_anno[$cos_id]);
		}
	}
	gzclose($cos_vcf_fp);
}



//init
$ngsbits = get_path("ngs-bits");
$genome_fa = genome_fasta($build, false);


list($vars_id_anno, $ann_cols) = parse_cmc_file($in_cmc, $build);
echo "Count of var_ids from cmc file: ".count($vars_id_anno)."\n";

// open output file
$temp_file = temp_file(".vcf", "cosmic_cmc_not_normalized");
$out_fp = fopen2($temp_file, "w");

// write VCF header
fwrite($out_fp, "##fileformat=VCFv4.2\n");
fwrite($out_fp, "##fileDate=".date("Ymd")."\n");
fwrite($out_fp, "##source={$in_cmc}\n");
fwrite($out_fp, "##reference={$genome_fa}\n");
fwrite($out_fp, "##INFO=<ID=COSMIC_CMC,Number=1,Type=String,Description=\"COSMIC Cancer Mutation Census (CMC) Data. Format: ". implode("|",$ann_cols)    ."\">\n");
fwrite($out_fp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

write_variants_from_file($vars_id_anno, $out_fp, $in_genome_vcf);
echo "Variants left after genome screen variants: ".count($vars_id_anno)."\n";

write_variants_from_file($vars_id_anno, $out_fp, $in_non_coding_vcf);
echo "Variants left after non-coding variants: ".count($vars_id_anno)."\n";

write_variants_from_file($vars_id_anno, $out_fp, $in_target_screens_vcf);
echo "Variants left after target_screens variants: ".count($vars_id_anno)."\n";

fclose($out_fp);

if(count($vars_id_anno) > 0)
{
	echo "Variants whose id was not found in the vcf file:\n";
	var_dump($vars_id_anno);
	trigger_error("Not all variants were found in the comsic vcf file! Not found: ".count($vars_id_anno), E_USER_WARNING);
}

//POST PROCESSING:

//Normalize file and sort 
$temp_file2 = temp_file(".vcf", "cosmic_cmc_normalized");
$parser->exec("{$ngsbits}/VcfLeftNormalize", "-stream -in $temp_file -out $temp_file2 -ref {$genome_fa}");

$temp_file3 = temp_file(".vcf", "cosmic_cmc_normalized_sorted");
$parser->exec("{$ngsbits}/VcfSort", " -in $temp_file2 -out $temp_file3");

//remove columns after INFO columns and compress
$parser->execPipeline([["cut -f1-8","$temp_file3"], ["bgzip" , " -c > $out"]], "trim and compress");
$parser->exec("tabix" , "-p vcf $out");

//Check converted VCF file
$parser->exec("{$ngsbits}/VcfCheck", "-in $out -ref {$genome_fa} -lines 0", true);

?>
