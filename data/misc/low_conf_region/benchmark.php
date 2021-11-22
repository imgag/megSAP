<?php

include("/mnt/users/bioinf/megSAP/src/Common/all.php");

//parse CLI arguments
$vcfs = [];
$regs = [];
for ($i=1; $i<$argc; ++$i)
{
	$file = trim($argv[$i]);
	if ($file=="") continue;
	
	if (ends_with($file, ".vcf"))
	{
		$vcfs[] = $file;
	}
	else
	{
		$regs[] = $file;
	}
}

//header
print "#reg	MB_removed_WGS	MB_removed_WES";
foreach($vcfs as $vcf)
{
	print "\t".basename($vcf, ".vcf");
}
print "\n";

list($stdout) = exec2("BedInfo -in /mnt/share/data/enrichment/WGS_hg38.bed");
$wgs_size = trim(explode(":", $stdout[1])[1]);


list($stdout) = exec2("BedInfo -in gene_exons.bed");
$wes_size = trim(explode(":", $stdout[1])[1]);

//content
foreach($regs as $reg)
{
	list($stdout) = exec2("BedInfo -in $reg");
	$reg_size = trim(explode(":", $stdout[1])[1]);
	print $reg."\t".number_format(($wgs_size - $reg_size)/1024/1024, 3);
	
	list($stdout) = exec2("BedIntersect -in $reg -in2 gene_exons.bed | BedInfo");
	$wes_size2 = trim(explode(":", $stdout[1])[1]);
	print "\t".number_format(($wes_size - $wes_size2)/1024/1024, 3);

	foreach($vcfs as $vcf)
	{
		list($vars) = exec2("egrep '^chr' {$vcf} | wc -l");
		$vars = trim(implode("", $vars));
		list($vars_filter) = exec2("VcfFilter -in {$vcf} -reg {$reg} | egrep '^chr' |wc -l");
		$vars_filter = trim(implode("", $vars_filter));
		print "\t".number_format(100.0*$vars_filter/$vars, 2);
	}
		
	print "\n";		
}

?>