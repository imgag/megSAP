<?php
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("vc_cfdna", "cfDNA mutation caller.");
$parser->addInfile("bam", "Input BAM file, with deduplicated alignments and DP tags.", false);
$parser->addInfile("target", "Target region BED file.", false);
$parser->addString("build", "Reference genome build.", false);
$parser->addOutfile("vcf", "Variant call output as VCF file.", false);
//$parser->addOutfile("tsv", "Variant call output as TSV file.", true);
//$parser->addOutfile("mrd", "MRD probability output file.", true);
$parser->addInfile("model", "Error model parameters.", true);
//$parser->addInfile("model_region", "Target region to use for parameter estimation.", true);

extract($parser->parse($argv));

$genome = genome_fasta($build);

// $workdir = $parser->tempFolder("cfdna_tmp");
// $workdir = dirname($vcf);
// $vcf_name = basename($vcf);
// $temp_vcf = "{$workdir}/{$vcf_name}";
$args = [
    "--bam", $bam,
    "--ref", $genome,
    "--bed", $target,
    "--out_file", $vcf
];
if (isset($model))
{
    $args[] = "--param";
    $args[] = $model;
}
else
{
    //check that target BED file has sufficient bases for error modeling
    $ret = $parser->exec(get_path("ngs-bits")."BedInfo", "-in {$target}");
    $n_bases = intval(explode(":", $ret[0][1])[1]);
    if ($n_bases < 1000)
    {
        trigger_error("Target region has to be at least 1,000 bases for reliable error model parameter estimation!", E_USER_WARNING);
    }
}

$parser->exec(get_path("cfdna_caller"), implode(" ", $args));

?>
