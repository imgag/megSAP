<?php

include("/mnt/storage2/megSAP/pipeline/src/Common/all.php");

$db = DB::getInstance("NGSD");

function infos(&$db, $ps)
{
	$info = get_processed_sample_info($db, $ps);

	$quality = $info['ps_quality'];
	$system = $info['sys_name_short'];
	$bam = $info['ps_bam'];

	return [$quality, $system, $bam];
}


print "#gsvar\tsystem\tps_c\tps_f\tps_m\n";
$gsvar_files = $db->getValues("SELECT gsvar_file FROM secondary_analysis WHERE type='trio'");
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
	$tmp = explode("_", $parts[count($parts)-2]);
	$ps_c = $tmp[1]."_".$tmp[2];
	$ps_f = $tmp[3]."_".$tmp[4];
	$ps_m = $tmp[5]."_".$tmp[6];
	list($q_c, $sys_c, $bam_c) = infos($db, $ps_c);
	list($q_f, $sys_f) = infos($db, $ps_f);
	list($q_m, $sys_m) = infos($db, $ps_m);
	if ($q_c=="bad" || $q_f=="bad" || $q_m=="bad")
	{
		print "##$gsvar skipped: bad quality ($q_c $q_f $q_m)\n";
		continue;
	}
	if ($sys_c!=$sys_f || $sys_c!=$sys_m)
	{
		print "##$gsvar skipped: differing systems ($sys_c $sys_f $sys_m)\n";
		continue;
	}
	
	if (!file_exists($bam_c))
	{
		print "##$gsvar skipped: index BAM missing\n";
		continue;
	}
	list($stdout) = exec2("SampleGender -in $bam_c -method sry -sry_cov 10");
	if (!contains($stdout[1], "female"))
	{
		print "##$gsvar skipped: male index (".$stdout[1].")\n";
		continue;
	}
	
	print $gsvar."\t".$sys_c."\t".$ps_c."\t".$ps_f."\t".$ps_m."\n";
}

?>