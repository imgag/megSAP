<?php

include("/mnt/storage2/megSAP/pipeline/src/Common/all.php");

print "#gsvar\tsystem\tinh_errors\tinh_error_perc\tcomments\n";
$db = DB::getInstance("NGSD");
$sys_filter = $argv[2];
$max_error_perc = $argv[3];
$out_folder = "trio_errors_{$sys_filter}";
exec2("mkdir -p $out_folder");

$file = file($argv[1]);
foreach($file as $line)
{
	$line = trim($line);
	if ($line=="" || $line[0]=="#") continue;
	
	list($gsvar, $sys, $ps_c, $ps_f, $ps_m) = explode("\t", $line);
	if ($sys!=$sys_filter) continue;
		
	$c_all = 0;
	$c_inh_error = 0;
	$output = [];
	$comments = [];
	$h = fopen2($gsvar, "r");
	while(!feof($h))
	{
		$line =  trim(fgets($h));
		if ($line=="") continue;
		
		$parts = explode("\t", $line);
		if ($line[0]=="#")
		{
			$c_idx = array_search($ps_c, $parts);
			$f_idx = array_search($ps_f, $parts);
			$m_idx = array_search($ps_m, $parts);
			continue;
		}
		++$c_all;
		
		$chr = $parts[0];
		$c = $parts[$c_idx];
		$f = $parts[$f_idx];
		$m = $parts[$m_idx];
		
		//skip all chromosomes that are not diploid - chrX is diploid as we look at female index cases only
		if($chr=="chrMT" || $chr=="chrY") continue; 
		
		$inh_error = false;
		if ($c=="wt" && ($f=="hom" || $m=="hom")) $inh_error = true;
		if ($c=="het" && (($f=="hom" && $m=="hom") || ($f=="wt" && $m=="wt"))) $inh_error = true;
		if ($c=="hom" && ($f=="wt" || $m=="wt")) $inh_error = true;
		if ($inh_error)
		{
			++$c_inh_error;
			$output[] = $line."\n";
		}
		
	}
	fclose($h);
	$error_perc = number_format(100.0 * $c_inh_error / $c_all, 2);
	if ($error_perc<=$max_error_perc)
	{
		file_put_contents("{$out_folder}/{$ps_c}_{$ps_f}_{$ps_m}.GSvar", $output);
	}
	else
	{
		$comments[] = "too many errors > skipped";
	}
	print $gsvar."\t".$sys."\t".$c_inh_error."\t".$error_perc."\t".implode("; ", $comments)."\n";
}

?>