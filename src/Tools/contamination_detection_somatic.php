<?php
/** 
	@page contamination_detection_somatic
	
	Test data: DX206781_01 and DX206802_01 are contaminated by the same DNA of unknown origin
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("contamination_detection_somatic", "Detects which sample could have caused the contamination of a tumor-normal pair.");
$parser->addString("t_id", "Processed sample id of contaminated tumor.", false);
$parser->addString("n_id", "Processed sample id of contaminated normal.", false);
$parser->addFloat("max_af", "Maximum variant allele frequency of variant in input tumor sample.", true, 0.2);
$parser->addFloat("min_af", "Minimum variant allele frequency of variant in input tumor sample.", true, 0.00);
$parser->addFloat("max_gnomad", "Maximum gnomAD allele frequency of variant in input tumor sample.", true, 0.01);
$parser->addInt("min_hits", "Minimum number of variant hits for reporting a potential contaminant sample.", true, 3);
$parser->addFloat("min_overlap_percentage", "Minimum variant overlap percentage with potential contaminant sample.", true, 0.05);
$parser->addInt("max_ngsd", "Maximum count inNGSD of variants to be considered.", true, 1000);
extract($parser->parse($argv));

if(!db_is_enabled("NGSD"))
{
	trigger_error("NGSD is needed for somatic contamination detection.", E_USER_ERROR);
}

//init
$db = DB::getInstance("NGSD");

//get required somatic input GSVar file
$tinfo = get_processed_sample_info($db, $t_id);
$gsvar_input = $tinfo["project_folder"] . "/Somatic_{$t_id}-{$n_id}/{$t_id}-{$n_id}.GSvar";
if(!file_exists($gsvar_input))
{
	trigger_error("Could not find $gsvar_input of contaminated pair {$t_id}-{$n_id}.", E_USER_ERROR);
}


$handle = fopen($gsvar_input,"r");
$i_gnomad = -1; //annotation index of gnomAD column

$c_vars_filtered = 0; //variant count of variants that pass gnomAD filtering of input variant list
$var_refound_count = array(); //count of variants that have been found in another sample in NGSD

while(!feof($handle))
{
	$line = fgets($handle);
	if( empty($line) ) continue;
	
	$parts = explode("\t", trim($line));
	
	//Determine index of gnomAD column
	if(starts_with($line, "#chr") )
	{
		for($i=0;$i<count($parts);++$i)
		{
			if($parts[$i] == "gnomAD") 
			{
				$i_gnomad  = $i;
				break;
			}
		}
	}
	
	if(starts_with($line, "#")) continue;
	
	list($chr,$start,$end,$ref,$obs,$tumor_af,$tumor_dp) = $parts;
	$gnomad_af = $parts[$i_gnomad];
	
	if($tumor_af > $max_af) continue;
	if(!empty($gnomad_af) && $gnomad_af > $max_gnomad) continue;
	
	
	//Get NGSD variant id
	$var_id = $db->getValue("SELECT id FROM variant WHERE chr='{$chr}' AND start='{$start}' AND end='{$end}' AND ref='{$ref}' AND obs='{$obs}'");
	if (is_null($var_id)) continue;	
	
	//determine germline samples that contain the same variant
	$res = $db->getValues("SELECT processed_sample_id FROM detected_variant WHERE variant_id={$var_id}");
	//Skip variants that occur too often in NGSD
	if(count($res) > $max_ngsd) continue;
	
	foreach($res as $ps_id)
	{
		if(array_key_exists($ps_id, $var_refound_count)) ++$var_refound_count[$ps_id];
		else $var_refound_count[$ps_id] = 1;
	}
	
	
	//determine somatic samples that contain the same variant
	$res2 = $db->getValues("SELECT processed_sample_id_tumor FROM detected_somatic_variant WHERE variant_id={$var_id}");
	foreach($res2 as $ps_id)
	{
		if(array_key_exists($ps_id, $var_refound_count)) ++$var_refound_count[$ps_id];
		else $var_refound_count[$ps_id] = 1;
	}
	
	++$c_vars_filtered;
}

fclose($handle);

arsort($var_refound_count);

//print output
print "##{$c_vars_filtered} variants are considered for detection after filtering and were found in NGSD.\n";
print "#sample\toverlap_perc\toverlap_count\n";
foreach($var_refound_count as $ps_id => $var_count)
{
	if ($var_count<$min_hits) continue;
	
	$overlap_perc = $var_count / $c_vars_filtered;
	
	if($overlap_perc < $min_overlap_percentage) continue;
	
	$ps_name = $db->getValue("SELECT CONCAT(s.name,'_',LPAD(ps.process_id,2,'0')) FROM processed_sample as ps, sample as s WHERE ps.sample_id = s.id AND ps.id={$ps_id}");
	
	//skip same sample
	if (explode("_", $t_id)[0]==explode("_", $ps_name)[0]) continue;
	
	print "{$ps_name}\t ". sprintf("%.3f", round($overlap_perc,3)) ."\t{$var_count}\n";
}

?>
