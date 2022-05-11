<?php
/** 
	@page find_same_sample 
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("find_same_sample", "Calculates the rare variant overlap for samples.");
$parser->addString("sample", "Processed sample name.", false);
$parser->addStringArray("test", "Processed sample names to samples to test.", false);
$parser->addFloat("max_af", "Maximum AF of variants to use for comparison in ExAC, gnomAD", true, 0.02);
extract($parser->parse($argv));

function load_rare_variants($sample, $max_af)
{
	//get infos
	$db_conn = DB::getInstance("NGSD");
	$info = get_processed_sample_info($db_conn, $sample, false);
	if (is_null($info)) return "No sample info in NGSD";
		
	//load variant list
	$gsvar_file = $info['ps_folder']."/".$info['ps_name'].".GSvar";
	if (!file_exists($gsvar_file)) return "GSvar file missing";
	
	$file = Matrix::fromTSV($gsvar_file);
	$vars_rare = array();
	$i_exa = $file->getColumnIndex("ExAC", false, false);
	if ($i_exa==-1) return "ExAC column missing";
	$i_gno = $file->getColumnIndex("gnomAD", false, false);
	if ($i_gno==-1) return "gnomAD column missing";
	$i_fil = $file->getColumnIndex("filter", false, false);
	if ($i_fil==-1) return "filter column missing";
	
	for($i=0; $i<$file->rows(); ++$i)
	{
		$row = $file->getRow($i);
		if (contains($row[$i_fil], "off-target")) continue;
		if ($row[$i_exa]==0 || $row[$i_exa]>$max_af) continue;
		if ($row[$i_gno]==0 || $row[$i_gno]>$max_af) continue;
		
		$vars_rare[] = implode(" ", array_slice($row, 0, 5));
	}
	
	return array($vars_rare, $file->rows(), $info['sys_name']);
}


//load reference sample variant list
list($vars1, $all_var_count1, $sys1) = load_rare_variants($sample, $max_af);
print "##$sample ($sys1): Select ".count($vars1)." rare, known variant from overall {$all_var_count1} variants.\n";
print "#sample\tprocessing_system\toverlap_abs\toverlap_perc\tcomments\n";
foreach($test as $sample2)
{
	if (starts_with($sample2, "#")) continue;
	
	$output = load_rare_variants($sample2, $max_af);
	if (is_array($output))
	{
		list($vars2, $all_var_count2, $sys2) = $output;
		$ol_abs = count(array_intersect($vars1, $vars2));
		$ol_perc = ($ol_abs==0) ? 0 :  number_format(100.0*$ol_abs/min(count($vars1), count($vars2)), 2);
		print "{$sample2}\t{$sys2}\t{$ol_abs}\t{$ol_perc}\t\n";
	}
	else
	{
		print "{$sample2}\t-\t0\t0.00\t{$output}\n";
	}
}

?>