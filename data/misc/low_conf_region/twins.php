<?php

include("/mnt/storage2/megSAP/pipeline/src/Common/all.php");

$db = DB::getInstance("NGSD");

print "#ps1\tps2\tsystem\tfolder\n";

$twin_samples = [];
$tmp = $db->executeQuery("SELECT `sample1_id`,`sample2_id`  FROM `sample_relations` WHERE `relation` = 'twins (monozygotic)'");
foreach($tmp as $row)
{
	$s1 = $db->getValue("SELECT name FROM sample WHERE id=".$row['sample1_id']);
	$s2 = $db->getValue("SELECT name FROM sample WHERE id=".$row['sample2_id']);
	$twin_samples[] = [$s1, $s2];
}

$gsvar_files = $db->getValues("SELECT gsvar_file FROM secondary_analysis WHERE type='multi sample'");
foreach($gsvar_files as $gsvar)
{
	if (!file_exists($gsvar))
	{
		print "##$gsvar skipped: missing file\n";
		continue;
	}
	if (!contains($gsvar, "/diagnostic/"))
	{
		print "##$gsvar skipped: not diagnostic\n";
		continue;
	}
	
	$parts = explode("/", $gsvar);
	$ps_concat = substr($parts[7], 6);
	$contains_samples = false;
	foreach($twin_samples as list($s1, $s2))
	{
		if (contains($ps_concat, $s1) && contains($ps_concat, $s2)) $contains_samples = true;
	}
	if (!$contains_samples) continue;
	
	$parts2 = explode("_", $ps_concat);
	if (count($parts2)!=4)
	{
		print "##$gsvar skipped: contains more than two samples\n";
		continue;
	}
	
	$ps1 = $parts2[0]."_".$parts2[1];
	$info1 = get_processed_sample_info($db, $ps1);
	$system = $info1['sys_name_short']; 
	$ps2 = $parts2[2]."_".$parts2[3];
	$info2 = get_processed_sample_info($db, $ps2);
	if ($info2['sys_name_short']!=$system)
	{
		print "##$gsvar skipped: processing system mismatch\n";
		continue;
	}
	
	print implode("\t",[$ps1, $ps2, $system, dirname($gsvar)])."\n";
	
}

?>