<?php
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("vc_cfdna", "cfDNA mutation caller.");
$parser->addInfile("bam", "Input BAM file, with deduplicated alignments and DP tags.", false);
$parser->addInfile("target", "Target region BED file.", false);
$parser->addString("build", "Reference genome build.", false);
$parser->addString("folder", "Analysis data folder with cfDNA variant calling.", false);
$parser->addInfile("model", "Error model parameters.", true);
$parser->addInfile("monitoring_vcf", "VCF containing monitoring SNVs and sample identifier.", true);

extract($parser->parse($argv));

$genome = genome_fasta($build);
$tempdir = $parser->tempFolder("umiVar");


// remove previous analysis
create_directory($folder);
$folder = realpath($folder);
$parser->exec("rm", "-rf $folder");

$args = [
    "--tbam", realpath($bam),
    "--ref", realpath($genome),
    "--bed", realpath($target),
    "--out_folder", $folder,
    "--temp_dir", $tempdir];
if (isset($monitoring_vcf))
{
    $args[] = "--monitoring ${monitoring_vcf}";
}
if (isset($model))
{
    $args[] = "--param {$model}";
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

//set environment variables
putenv("umiVar_python_binary=\"".get_path("python3")."\"");
putenv("umiVar_R_binary=\"".get_path("rscript")."\"");
putenv("umiVar_samtools_binary=\"".get_path("samtools")."\"");

// call umiVar2 in virtual environment
$umiVar2 = get_path("umiVar2");
$parser->exec(get_path("python3"), $umiVar2."/umiVar.py ".implode(" ", $args));

// sort VCF file(s)
$vcf = $folder."/".basename2($bam).".vcf";
if (file_exists($vcf))
{
    $parser->exec(get_path("ngs-bits")."VcfSort","-in $vcf -out $vcf", true);
}
$vcf_hq = $folder."/".basename2($bam)."_hq.vcf";
if (file_exists($vcf_hq))
{
    $parser->exec(get_path("ngs-bits")."VcfSort","-in ${vcf_hq} -out ${vcf_hq}", true);
}


?>
