<?php

/** 
	@page hla_type
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("hla_genotyper", "Predicts four-digit HLA class I and II genotypes for a sample.");
$parser->addInfile("bam",  "Input BAM file.", false);
$parser->addOutfile("out",  "Output file with HLA haplotypes. If not specified, writes to stdout", true);
//optional params
$parser->addString("genes", "comma-separated list of HLA genes to be calculated (A,B,C,DRB1,DQA1,DQB1)", true, "A,B,C");
$all_genes = array("A", "B", "C", "DRB1", "DQA1", "DQB1");
$parser->addEnum("type", "Type of the analysis (exome,genome,rnaseq)", true, array("exome", "genome", "rnaseq"), "exome");
$parser->addEnum("ethnicity", "Ethnicity of sample.", true, array("EUR","AFA","HIS","API","UNK","UNK"), "UNK");
$parser->addString("name", "Sample name for output file.", true, "SAMPLE");
$parser->addEnum("build", "The genome build to use.", true, array("GRCh37", "GRCh38"), "GRCh38") ;
$parser->addOutfile("out_dose", "Output file with per allele dose.", true);
extract($parser->parse($argv));

//Check HLA genes
$genes = explode(",", $genes);
foreach($genes as $gene)
{
	if(!in_array($gene, $all_genes))
	{
		trigger_error("This tool only supports the class I and II HLA-genes " . implode(", ", $all_genes), E_USER_ERROR);
	}
}


//Execute hla-genotyper
$outdir = temp_folder();
$args = array(
	"-o", $outdir,
	"-e", $ethnicity,
	"-r", ($build == "GRCh38" ? "38" : "37"),
	"-s", $name,
	"-g", implode(",", $genes),
	"--{$type}"
);
$parser->exec(get_path("python3")." ". get_path("hla_genotyper"), implode(" ", $args) . " $bam");

//Print output
$res_prefix = "{$outdir}/hla.{$ethnicity}.{$name}.{$type}";
if(isset($out))
{
	$parser->copyFile("{$res_prefix}.txt", $out);
}
else
{
	$file_data = file("{$res_prefix}.txt");
	$out_stream = fopen2("php://stdout", "w");
	
	foreach($file_data as $line)
	{
		fwrite($out_stream, $line);
	}
	fclose($out_stream);
}

if(isset($out_dose))
{
	$parser->copyFile("{$res_prefix}.dose" , $out_dose);
}

?>

