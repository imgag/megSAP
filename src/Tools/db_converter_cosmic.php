<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("db_converter_cosmic", "Converts COSMIC Cancer Mutation Census (CMC) database file (TSV) to tabix-indexed VCF.GZ.");
// optional
$parser->addInfile("in",  "Input file in TSV-format with COSMIC CMC scores. ('-' for STDIN)", false, false);
$parser->addOutfile("out",  "Output file in vcf.gz-format with CMC data.", false, false);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
extract($parser->parse($argv));

//check
if($build != "GRCh37" && $build != "GRCh38")
{
	trigger_error("COSMIC only contains annotation coordinates for GRCh37 and GRCh38.", E_USER_ERROR);
}

//init
$ngsbits = get_path("ngs-bits");
$genome_fa = genome_fasta($build, false);


// open input file
if ($in == "-") $in = "php://stdin";
$in_fp = fopen2($in, "r");


// open output file
$temp_file = temp_file(".vcf", "cosmic_cmc_not_normalized");
$out_fp = fopen2($temp_file, "w");

// write VCF header
fwrite($out_fp, "##fileformat=VCFv4.2\n");
fwrite($out_fp, "##fileDate=".date("Ymd")."\n");
fwrite($out_fp, "##source={$in}\n");
fwrite($out_fp, "##reference={$genome_fa}\n");


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

fwrite($out_fp, "##INFO=<ID=COSMIC_CMC,Number=1,Type=String,Description=\"COSMIC Cancer Mutation Census (CMC) Data. Format: ". implode("|",$ann_cols)    ."\">\n");
fwrite($out_fp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

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
	
	
	//coordinates in format 1:124-214321
	$coords = $parts[ $indices["Mutation genome position $build"] ];
	list($chr, $start, $end) = preg_split("/[:-]/", $coords);
	
	
	if( empty($chr) || empty($start) || empty($end) ) continue; //Skip entries without genomic coordinates
	
	//convert chromsome names
	$chr = "chr" . $chr;
	$chr = str_replace("chr23", "chrX", $chr);
	$chr = str_replace("chr24", "chrY", $chr);
	$chr = str_replace("chr25", "chrMT", $chr);
	
	
	$pos = (int) $start;
	$ref = trim( $parts[$indices["GENOMIC_WT_ALLELE_SEQ"]] );
	$alt = trim( $parts[$indices["GENOMIC_MUT_ALLELE_SEQ"]] );
	
	$cmc_data = array();
	foreach($ann_cols as $col_name)
	{
		$temp = $parts[$indices[$col_name]];
		$temp = str_replace(",", ";", $temp); //replace "," by ";", otherwise there is no in-field separator for later processing left.
		$cmc_data[] = vcf_encode_url_string($temp);
	}
	
	fwrite($out_fp, "$chr\t$pos\t.\t$ref\t$alt\t.\t.\tCOSMIC_CMC=". implode("|", $cmc_data) ."\n");
}

fclose($in_fp);
fclose($out_fp);

//Normalize file and sort 
$temp_file2 = temp_file(".vcf", "cosmic_cmc_normalized");
$parser->exec("{$ngsbits}/VcfLeftNormalize", " -in $temp_file -out $temp_file2 -ref {$genome_fa}");

$temp_file3 = temp_file(".vcf", "cosmic_cmc_normalized_sorted");
$parser->exec("{$ngsbits}/VcfSort", " -in $temp_file2 -out $temp_file3");

//remove columns after INFO columns and compress
$parser->execPipeline([["cut -f1-8","$temp_file3"], ["bgzip" , " -c > $out"]], "trim and compress");
$parser->exec("tabix" , "-p vcf $out");

//Check converted VCF file
$parser->exec("{$ngsbits}/VcfCheck", "-in $out -ref {$genome_fa}", true);

?>
