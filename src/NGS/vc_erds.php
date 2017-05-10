<?php 
/** 
	@page erds
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("erds", "Run ERDS. Duplication / Deletion / Breakpoint detection. Only WGS.");
$parser->addInfile("bam_file",  "Bam file.", false);
$parser->addInfile("vcf_file",  "VCF file.", false);
$parser->addString("out", "Output folder.", false);
//optional
$parser->addString("build", "The genome build to use.", true, "GRCh37");
extract($parser->parse($argv));

//TODO check processing system => only WG
//TODO write to temp folder and copy only event files
$processing_id = basename($bam_file, ".bam");
if(!is_dir($out))	trigger_error("Could not find output folder ".$out, E_USER_ERROR);
$out .= "/".$processing_id;
if(!is_dir($out))	mkdir($out);
chdir($out);	//move to out_folder
$parser->setLogFile(realpath($parser->getLogFile()));				//redirect log file
$parser->exec(get_path("ERDS"), "-b $bam_file -v $vcf_file -out $out -r ".get_path("local_data")."/{$build}.fa", true);

//annotate genes, use sample name from vcf file
list($sid,) = explode("_",$processing_id);
$event_file = $out."/ut/".$sid.".events";
$parser->exec(get_path("ngs-bits")."BedAnnotateGenes", "-in $event_file -out $event_file",true);

//add headers
$events = Matrix::fromTSV($event_file);
$events->setHeaders(array("chr","start","end","length","type","summed_score","precise_boundary","reference_cn","inferred_cn","genes"));
$events->toTSV($event_file);
?>
