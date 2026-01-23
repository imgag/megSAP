<?php 
/** 
	@page normalize_small_variants
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("normalize_small_variants", "Variant normalization pipeline.");
$parser->addInfile("in",  "Input VCF/VCF.GZ file .", false);
$parser->addOutfile("out", "Output file in VCF.", false);

//optional
$parser->addString("build", "The genome build to use.", true, "GRCh38");
$parser->addString("mode", "Set mode for vcf_fix.php (deepvariant, mosaic, clair3, dragen, mixed)", true, "freebayes");
$parser->addFlag("primitives", "Enables vcfallelicprimitives");
$parser->addFlag("fix", "Enables vcf_fix.php step");

extract($parser->parse($argv));

//init
$genome = genome_fasta($build);

//filter variants according to variant quality>5
$pipeline[] = array(ends_with(strtolower($in), ".vcf") ? "cat" : "zcat", $in);

//Annotate source variants
$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfStoreSourceVariant", "", [], [], true)];

//split complex variants to primitives
//this step has to be performed before VcfBreakMulti - otherwise mulitallelic variants that contain both 'hom' and 'het' genotypes fail - see NA12878 amplicon test chr2:215632236-215632276
if ($primitives) $pipeline[] = ["", $parser->execApptainer("vcflib", "vcfallelicprimitives", "-kg", [], [], true)];

//split multi-allelic variants - -no_errors flag can be removed, when vcfallelicprimitives is replaced
$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfBreakMulti", "-no_errors", [], [], true)];
//normalize all variants and align INDELs to the left
$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfLeftNormalize", "-stream -ref $genome", [$genome], [], true)];

//sort variants by genomic position
$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfStreamSort", "", [], [], true)];

//fix errors and merge variants
$mode_arg = "";
if ($mode != "freebayes") $mode_arg = "--{$mode}_mode";
if ($mode == "mixed") $mode_arg = "--dragen_mode";
if ($fix) $pipeline[] = ["php ".repository_basedir()."/src/Tools/vcf_fix.php", "{$mode_arg}", false];
if ($fix && $mode == "mixed") $pipeline[] = ["php ".repository_basedir()."/src/Tools/vcf_fix.php", "--clair3_mode", false];

//Prune unchanged source variant annotation
$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfPruneSourceVariant", "-ref $genome -out $out", [$genome], [], true)];

//execute pipeline
$parser->execPipeline($pipeline, "Variant normalization pipeline");

?>
