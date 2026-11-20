<?php 
/** 
	@page an_aidiva
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");
require_once(__DIR__."/an_aidiva_common.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("an_aidiva", "Annotates a VCF file with aiDIVA score (independent of pheontype). INFO tag in output: AIDIVA");
$parser->addInfile("in",  "Input file in VCF or VCF.GZ format.", false);
$parser->addOutfile("out", "Output file in VCF format.", false);
extract($parser->parse($argv));

//convert VCF to intput TSV
$tsv_in = $parser->tempFile(".tsv");
$c_vars = aidiva_vcf2tsv($in, $tsv_in);

//run aiDIVA only if there are variants to score - no variants leads to missing aiDIVA output file
$scores = [];
if ($c_vars>0)
{
	//create config file based on template
	$model = get_path("data_folder")."/dbs/aiDIVA/aiDIVA_megSAP_model.pkl.gz"; 
	$config = aidiva_config(repository_basedir()."/data/misc/aidiva/config_score.yaml", $model);

	//run aiDIVA
	$out_prefix = $parser->tempFolder()."/out";
	$parser->execApptainer("aidiva", "python3", "/opt/aidiva_local/aidiva/run_aidiva-score.py --in_data $tsv_in --config $config --out_prefix $out_prefix", [$model]);

	//get aiDIVA scores from output TSV
	$scores = aidiva_results("{$out_prefix}_result_aidiva-score.tsv", ['AIDIVA_SCORE']);
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
		$parts[7] .= ";AIDIVA=".number_format($scores[$tag]['AIDIVA_SCORE'], 3);
	}
	fputs($ho, implode("\t", $parts)."\n");
}
gzclose($h);
fclose($ho);

?>