<?php 
/** 
	@page an_somatic_cancerhotspots
	
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("an_somatic_cancerhotspots", "Variant annotation with somatc cancer hot spots from cancerhotspots.org.");
$parser->addInfile("in",  "Input file in VCF format (VEP annotation mandatory).", false);
$parser->addOutfile("out", "Output file in VCF format.", false);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
extract($parser->parse($argv));

//get index of columnn in QSC header.
function index_of($cols, $name)
{
	$index = array_search($name, $cols);
	if ($index===FALSE)
	{
		trigger_error("Could not find column '$name' in QSC2 annotation. Valid column names are: ".implode(", ", array_values($cols)), E_USER_ERROR);
	}
	return $index;
}

$i_vep_symbol = -1;
$i_vep_hgvsp = -1;
$i_vep_feature = -1;

$in_db = get_path("data_folder") . "/dbs/cancerhotspots/cancerhotspots_snv.tsv";
if(!file_exists($in_db))
{
	trigger_error("Could not find cancerhotspots annotation file {$in_db}." , E_USER_ERROR);
}

$handle_in = fopen2($in, "r");
$handle_out = fopen2($out, "w");

while(!feof($handle_in))
{
	$line_in = trim( fgets($handle_in) );
	
	
	if(empty($line_in))
	{
		continue;
	}
	
	//Parse header
	if(starts_with($line_in, "#"))
	{
		//get annotation indices in CSQ field
		if (starts_with($line_in, "##INFO=<ID=CSQ2,"))
		{
			$cols = explode("|", substr($line_in, 0, -2));
			$i_vep_symbol = index_of($cols, "SYMBOL");
			$i_vep_hgvsp = index_of($cols, "HGVSp");
			$i_vep_feature = index_of($cols, "Feature");
			
			//write INFO field for CANCERHOTSPOTS just before INFO field of CSQ annotation
			fwrite($handle_out,"##INFO=<ID=CANCERHOTSPOTS,Number=.,Type=String,Description=\"Cancerhotspots.org data. Format: GENE_SYMBOL|ENSEMBL_TRANSCRIPT|AA_POS|AA_REF|AA_ALT|TOTAL_COUNT|ALT_COUNT\">\n");
		}
		//write info field of CSQ
		fwrite($handle_out, "{$line_in}\n");
		continue;
	}
	
	//parse variants	
	$parts = explode("\t", $line_in);
	if(count($parts) <9)
	{
		trigger_error("Data line with less than 9 fields found: {$line_in}.", E_USER_ERROR);
	}
	
	//Make list of amino acid changes predicted from VEP for this variant
	$vep_alterations = array();
	$info_parts = explode(";", $parts[7]);
	foreach($info_parts as $info_part)
	{
		if(!starts_with($info_part, "CSQ2=")) continue;
		
		$csq_parts = explode(",", $info_part);
		foreach($csq_parts as $csq_part)
		{
			$vep_parts = explode("|", str_replace("CSQ2=","", $csq_part) );
			
			
			$vep_gene = $vep_parts [$i_vep_symbol];
			$vep_trans = explode(".", $vep_parts [$i_vep_feature].".")[0]; # remove transcript version if annotated
			
			//Parse VEP amino acid change
			$vep_hgvsp_parts = explode(":", vcf_decode_url_string( $vep_parts [$i_vep_hgvsp] ) );
			if(count($vep_hgvsp_parts) < 2 ) continue;
			$vep_raw_aa_change = $vep_hgvsp_parts[1];
			
			//Skip protein changes other than p.AAA21321AAA (i.e. substitutions)
			if(!preg_match("/^p.[a-zA-Z]{3}\d+[a-zA-Z]{3}$/", $vep_raw_aa_change)  ) continue;
			
			$vep_aa_ref = substr( str_replace("p.","", $vep_raw_aa_change), 0,3 );
			$vep_aa_alt = substr( $vep_raw_aa_change, -3 );
			$vep_aa_pos = preg_replace( "/[.a-zA-Z]/", "", $vep_raw_aa_change ); //removes alphabetic chars and "." from unparsed string, only aa position remains
			
			
			$vep_alterations[] = array("gene" => $vep_gene, "transcript" => $vep_trans, "aa_ref" => aa3_to_aa1($vep_aa_ref), "aa_pos" => $vep_aa_pos, "aa_alt" => aa3_to_aa1($vep_aa_alt) );
		}
	}
	
	//query database file
	$matches = array();
	foreach($vep_alterations as $vep_alt)
	{
		list($res,,$exit_code) = exec2("grep " . $vep_alt["transcript"] ." $in_db", false);
		
		if($exit_code != 0) continue; //grep returns 1 if entry not found
		
		foreach($res as $res2)
		{
			list($gene,$entrez_id,$trans_id,$pos,$ref,$alt,$total_count, $alt_count) = explode("\t", $res2);
			
			if($pos != $vep_alt["aa_pos"] || $ref != $vep_alt["aa_ref"] || $alt != $vep_alt["aa_alt"]) continue;
			
			$matches[] = implode("|",[$gene,$trans_id,$pos,$ref,$alt,$total_count,$alt_count]);
		}
	}
	
	//add matches from cancerhotspots input file to INFO field
	if(!empty($matches)) $parts[7] .= ";CANCERHOTSPOTS=" . implode("," , $matches);
	fwrite( $handle_out, implode("\t", $parts) ."\n" );
}

//check vcf file
$parser->exec(get_path("ngs-bits")."VcfCheck", "-in $out -ref ".genome_fasta($build), true);
?>
