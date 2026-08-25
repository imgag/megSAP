<?php 
/** 
	@page vc_himito
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("vc_himito", "Mito variant calling with Himito.");
$parser->addInfile("bam",  "Indexed BAM/CRAM file to call variants on.", false);
$parser->addOutfile("out", "Output file in VCF.GZ format.", false);
$parser->addEnum("data_type", "Data type of the input data", false, ["pacbio", "ont-r9", "ont-r10"]);
//optional
$parser->addInt("threads", "The maximum number of threads used.", true, 1);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
$parser->addFlag("keep_additional_data", "Keep all calling files in separate folder (same directory as '-out')");
$parser->addString("name", "Sample name used for output VCF. If unset, BAM base name is used.", true, "");
extract($parser->parse($argv));

//init
$genome = genome_fasta($build);
$gz_output = ends_with(strtolower($out), ".vcf.gz");
$tmp_folder = $parser->tempFolder("himito");
if ($name=="") $name = basename2($bam);

//extract chrMT
$genome_mito = $parser->tempFile("chrMT.fa");
$parser->execApptainer("samtools", "samtools faidx", "-o {$genome_mito} {$genome} chrMT", [$genome]);
$parser->execApptainer("samtools", "samtools faidx", "{$genome_mito}");

//run himito
$args = [];
$args[] = "--input-bam {$bam}";
$args[] = "--threads {$threads}";
$args[] = "--chromo chrMT";
$args[] = "--reference-path {$genome_mito}";
$args[] = "--data-type {$data_type}";
$args[] = "--output-prefix {$tmp_folder}/{$name}";
$args[] = "--sample-id {$name}";
$parser->execApptainer("Himito", "Himito quick-start", implode(" ", $args), [$bam]);

//post-processing
$pipeline = [];
$pipeline[] = ["cat", "{$tmp_folder}/{$name}.vcf"];
$pipeline[] = ["sed", "'s/^chrM/chrMT/1'"];
$pipeline[] = ["", $parser->execApptainer("vcflib", "vcfallelicprimitives", "-kg", [], [], true)];
$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfBreakMulti", "", [], [], true)];
$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfLeftNormalize", "-stream -ref $genome", [$genome], [], true)];
if ($gz_output) 
{
	$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfStreamSort", "", [], [], true)];
	$pipeline[] = ["", $parser->execApptainer("htslib", "bgzip", "-c > $out", [], [dirname($out)], true)];
}
else
{
	$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfStreamSort", "-out {$out}", [], [dirname($out)], true)];
}
$parser->execPipeline($pipeline, "Himito post processing");

if ($gz_output) $parser->execApptainer("htslib", "tabix", "-f -p vcf $out", [], [dirname($out)]);

if ($keep_additional_data)
{
	//copy additional data
	$out_folder = dirname($out)."/himito/";
	if (!file_exists($out_folder)) mkdir($out_folder);
	$parser->exec("cp", "{$tmp_folder}/* {$out_folder}/");
}

?>
