<?php 
/** 
	@page an_aidiva_hpo
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");
require_once(__DIR__."/an_aidiva_common.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("an_aidiva_hpo", "Ranks variants in a VCF using the phenotype-specific aiDIVA score. INFO tags in output: FINAL_AIDIVA_SCORE and AIDIVA_RANK");
$parser->addInfile("in",  "Input file in VCF or VCF.GZ format that is annotated with 'an_aidiva'.", false);
$parser->addOutfile("out", "Output file in VCF format.", false);
$parser->addString("hpo_list", "Comma-speparated list of HPO terms: HP:XXXXXX,HP:XXXXXX", false);
$parser->addInt("threads", "Number of threads to use.", true, 1);
extract($parser->parse($argv));

//convert VCF to intput TSV
$tsv_in = $parser->tempFile(".tsv");
$c_vars = aidiva_vcf2tsv($in, $tsv_in, true);

//run aiDIVA only if there are variants to score - no variants leads to missing aiDIVA output file
$infos = [];
if ($c_vars>0)
{
	//create config file based on template
	$model = get_path("data_folder")."/dbs/aiDIVA/aiDIVA_megSAP_model.pkl.gz";
	$config = aidiva_config(repository_basedir()."/data/misc/aidiva/config_score_hpo.yaml", $model);

	//run aiDIVA
	$out_prefix = $parser->tempFolder()."/out";
	$parser->execApptainer("aidiva", "python3", "/opt/aidiva_local/aidiva/run_aidiva-rf.py --in_data $tsv_in --config $config --out_prefix $out_prefix --hpo_list $hpo_list --threads $threads", [$model]);

	//get aiDIVA infos from output TSV
	$infos = aidiva_results("{$out_prefix}_result_filtered_aidiva-rf.tsv", ['FINAL_AIDIVA_SCORE', 'AIDIVA_RANK']);
}

//annotate AIDIVA infos to VCF
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
			fputs($ho, "##INFO=<ID=FINAL_AIDIVA_SCORE,Number=1,Type=Float,Description=\"aiDIVA pathogenicity score (phenotype-specific).\">\n");
			fputs($ho, "##INFO=<ID=AIDIVA_RANK,Number=1,Type=Int,Description=\"aiDIVA rank.\">\n");
		}
		fputs($ho, $line."\n");
		continue;
	}
	
	$parts = explode("\t", $line);
	list($chr, $pos, $id, $ref, $alt) = $parts;
	$tag = "{$chr}_{$pos}_{$ref}_{$alt}";
	if (isset($infos[$tag]))
	{
		$parts[7] .= ";FINAL_AIDIVA_SCORE=".number_format($infos[$tag]['FINAL_AIDIVA_SCORE'], 3).";AIDIVA_RANK=".(int)($infos[$tag]['AIDIVA_RANK']);
	}
	fputs($ho, implode("\t", $parts)."\n");
}
gzclose($h);
fclose($ho);

?>
