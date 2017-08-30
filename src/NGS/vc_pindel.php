<?php 
/** 
	@page pindel
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_pindel", "Run pindel. No Haloplex.");
$parser->addInfileArray("bam_files",  "Bam files for pindel. Please note that MappingQC output file must reside next to bam file to extract correct insert size.", false);
$parser->addString("out", "Output folder (will contain files prefixed with pindel).", false);
//optional
$parser->addString("build", "The genome build to use.", true, "GRCh37");
extract($parser->parse($argv));

//check bam file and qcml file
if(!is_dir($out))	trigger_error("Output folder '".$out."' is not valid!", E_USER_ERROR);
foreach($bam_files as $bam_file)
{
	if(!is_file($bam_file))	trigger_error("Could not find bam file '".$bam_file."'!", E_USER_ERROR);
	$qcml_file = dirname($bam_file)."/".basename($bam_file, ".bam")."_stats_map.qcML";
	if(!is_file($qcml_file))	trigger_error("Could not find mapping qcml file '".$qcml_file."'!", E_USER_ERROR);
}

//create config file
foreach($bam_files as $bam_file)
{
	//get insert size
	$insert_size_name = "insert size";
	$insert_size_id = "QC:2000023";
	$qcml_file = dirname($bam_file)."/".basename($bam_file, ".bam")."_stats_map.qcML";
	$is = get_qc_from_qcml($qcml_file, $insert_size_id, $insert_size_name);
	if($is===FALSE)
	{
		trigger_error("QC value with ID ".$insert_size_id." was not found within qcml file ".$qcml_file.". Skipping this sample.",E_USER_WARNING);
		continue;
	}
	
	//run pindel for this sample, nb unable to run pindel without a config file and calculated insert size
	//TODO move _CONF file to temp folder
	$p_id = basename($bam_file, ".bam");
	$conf_file = $out."/".$p_id."_CONF"; //$parser->tempFile("_CONF");
	$conf = new Matrix();
	$conf->addRow(array($bam_file,$is,$p_id));
	$conf->toTSV($conf_file);
	$parser->exec(get_path("pindel"), "-f ".get_path("local_data")."/{$build}.fa -i ".$conf_file." -c ALL -o ".$out."/".$p_id,true);
}
$conf->toTSV($conf_file);

//combine all results
$combined_results = new Matrix();
foreach($bam_files as $bam_file)
{
	$p_id = basename($bam_file, ".bam");

	//different SV types
	//_D = deletions, _INV = inversions, _SI = small indels, TD = tandem duplications, LI = long indels, BP = unassigned breakpoints
	$s_variants = array("_D","_INV","_SI","_TD","_LI","_BP");
	foreach($s_variants as $s_variant)
	{
		$handle = fopen($out."/".$p_id.$s_variant, "r");
		while (!feof($handle)) 
		{
			$row = fgets($handle);
			if($row[0]=="#")
			{
				$row = trim(fgets($handle));
				$columns = explode("\t",$row);
				$chr = substr($columns[3],5);	//ChrID 
				$start = substr($columns[4],2);	//BP 
				$end = $columns[5];
				$type = trim($s_variant,'_');
				$support = substr($columns[8],8);
				$tmp = $s_variant;
				if($s_variant=="_SI")	$tmp = "I";
				$size = substr($columns[1],strlen($tmp));	//TD 126, currently only _SI variants found / maybe update need for other types
				$id = $p_id;
				$combined_results->addRow(array($chr,$start,$end,"",$type,trim($support),$size,$id));
			}
		}
		fclose($handle);
	}
}
$result_file = $out."/_combined_results.tsv";
$combined_results->setHeaders(array("chr","start","end","gene","type","supporting_reads","region_size","sample_id"));
$combined_results->toTSV($result_file);

//annotate genes
$parser->exec(get_path("ngs-bits")."BedAnnotateGenes", " -in $result_file -out $result_file", true);
$parser->exec(get_path("ngs-bits")."BedSort", " -in $result_file -out $result_file", true);
//add headers
$combined_results = Matrix::fromTSV($result_file);
$combined_results->toTSV($result_file);
