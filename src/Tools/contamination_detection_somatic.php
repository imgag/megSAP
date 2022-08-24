<?php
/** 
	@page contamination_detection_somatic
	
	Test data:
	- DX206781_01 and DX206802_01 are contaminated by the same DNA of unknown origin
	- DNA2206018A1_01/DNA2205431A1_01 is contaminated by DNA2206062A1_01
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("contamination_detection_somatic", "Detects which sample could have caused the contamination of a tumor sample normal pair.");
$parser->addString("gsvar", "GSvar file of tumor/normal analysis.", false);
$parser->addFloat("max_af", "Maximum variant allele frequency of variant in input tumor sample.", true, 0.3);
$parser->addFloat("max_gnomad", "Maximum gnomAD allele frequency of variant in input tumor sample.", true, 0.05);
$parser->addInt("max_ngsd", "Maximum occurances in NGSD for variants to be considered.", true, 1000);
$parser->addFloat("min_overlap_percentage", "Minimum variant overlap percentage to display potential contaminant sample.", true, 8.0);
extract($parser->parse($argv));

//init
if(!db_is_enabled("NGSD"))
{
	trigger_error("NGSD is needed for somatic contamination detection.", E_USER_ERROR);
}
$db = DB::getInstance("NGSD");

$i_gnomad = -1; //annotation index of gnomAD column
$c_vars_used = 0; //count of variants used for comparison
$var_refound_by_ps_id = array(); //count of variants that have been found in another sample in NGSD

$handle = fopen2($gsvar,"r");
while(!feof($handle))
{
	$line = fgets($handle);
	if(nl_trim($line)=="") continue;
	
	$parts = explode("\t", trim($line));
	
	//header
	if(starts_with($line, "#"))
	{
		//determine index of gnomAD column
		if(starts_with($line, "#chr"))
		{
			for($i=0;$i<count($parts);++$i)
			{
				if($parts[$i] == "gnomAD") $i_gnomad  = $i;
			}
			if ($i_gnomad==-1) trigger_error("Column 'gnomAD' not found in tumor/normal GSvar file!");
		}
	
		continue;
	}
	
	//content
	list($chr, $start, $end, $ref, $obs, $tumor_af, $tumor_dp) = $parts;
	
	//filter for tumor AF
	if($tumor_af > $max_af) continue;
	
	//filter for gnomAD AF
	$gnomad_af = trim($parts[$i_gnomad]);
	if ($gnomad_af=="") continue; //we want known germline variants
	if($gnomad_af > $max_gnomad) continue; //we want rare variants
	
	//get NGSD variant id
	$var_id = $db->getValue("SELECT id FROM variant WHERE chr='{$chr}' AND start='{$start}' AND end='{$end}' AND ref='{$ref}' AND obs='{$obs}'", -1);
	if ($var_id==-1) continue;	
	
	//determine germline samples that contain the same variant
	$res = $db->getValues("SELECT processed_sample_id FROM detected_variant WHERE variant_id={$var_id} LIMIT ".($max_ngsd+1));
	if(count($res) > $max_ngsd) continue; //Skip variants that occur too often in NGSD (artefacts)
	foreach($res as $ps_id)
	{
		if(!isset($var_refound_by_ps_id[$ps_id])) $var_refound_by_ps_id[$ps_id] = 0;
		++$var_refound_by_ps_id[$ps_id];
	}
	
	//determine somatic samples that contain the same variant
	$res = $db->getValues("SELECT processed_sample_id_tumor FROM detected_somatic_variant WHERE variant_id={$var_id} LIMIT ".($max_ngsd+1));
	foreach($res as $ps_id)
	{
		if(!isset($var_refound_by_ps_id[$ps_id])) $var_refound_by_ps_id[$ps_id] = 0;
		++$var_refound_by_ps_id[$ps_id];
	}
	
	++$c_vars_used;
}
fclose($handle);

arsort($var_refound_by_ps_id);

//print output
print "##{$c_vars_used} variants are considered for detection of contaminant sample.\n";
print "#sample\toverlap_perc\toverlap_count\n";
foreach($var_refound_by_ps_id as $ps_id => $vars_refound)
{	
	$overlap_perc = 100.0*$vars_refound / $c_vars_used;
	if($overlap_perc < $min_overlap_percentage) continue;
	
	$ps_name = $db->getValue("SELECT CONCAT(s.name,'_',LPAD(ps.process_id,2,'0')) FROM processed_sample as ps, sample as s WHERE ps.sample_id = s.id AND ps.id={$ps_id}");
	if (contains($gsvar, $ps_name)) continue; //skip tumor and normal sample
		
	print "{$ps_name}\t ". number_format($overlap_perc, 2) ."\t{$vars_refound}\n";
}

?>
