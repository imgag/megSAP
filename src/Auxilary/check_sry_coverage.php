<?php
/** 
	@page check_sry_coverage
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("check_sry_coverage", "Checks SRY coverage of tumor for female samples, can be a hint for contamination with male DNA.");
$parser->addInfile("t_bam", "Tumor sample to check SRY coverage for.", false);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
extract($parser->parse($argv));

$t_id = basename2($t_bam);
$genome = genome_fasta($build);

if( db_is_enabled("NGSD") )
{
	$db = DB::getInstance("NGSD");
	$tinfo = get_processed_sample_info($db, $t_id, false);
	
	if(!is_null($tinfo) && $tinfo["gender"] == "female")
	{
        $out = $parser->execApptainer("ngs-bits", "SampleGender", "-in $t_bam -build ".ngsbits_build($build)." -method sry -ref {$genome}", [$t_bam, $genome]);
        list(,,$cov_sry) = explode("\t", $out[0][1]);

        if(is_numeric($cov_sry) && (float)$cov_sry >= 30)
        {
            trigger_error("Detected contamination of female tumor sample {$t_id} with male genomic DNA on SRY. SRY coverage is at {$cov_sry}x.", E_USER_ERROR);
        }
	}
}
else
{
    trigger_error("No NGSD access, skipping check of female tumor sample  $t_bam for contamination with male genomic DNA.", E_USER_WARNING);
}

?>