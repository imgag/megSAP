<?php

include("../../../src/Common/all.php");

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

$wgs_size = bed_size("/mnt/storage2/megSAP/data/enrichment/WGS_hg38.bed");
$wes_size = bed_size("gene_exons.bed");

//content
foreach($regs as $reg)
{
	print $reg;
	
	list($stdout) = exec2("BedIntersect -in $reg -in2 /mnt/storage2/megSAP/data/enrichment/WGS_hg38.bed | BedInfo");
	$weg_removed = trim(explode(":", $stdout[1])[1]);
	print "\t".number_format($weg_removed/1000000, 3, '.', '');
	
	list($stdout) = exec2("BedIntersect -in $reg -in2 gene_exons.bed | BedInfo");
	$wes_removed = trim(explode(":", $stdout[1])[1]);
	print "\t".number_format($wes_removed/1000000, 3);

	foreach($vcfs as $vcf)
	{
		list($vars) = exec2("egrep '^chr' {$vcf} | wc -l");
		$vars = trim(implode("", $vars));
		list($vars_filter) = exec2("VcfFilter -in {$vcf} -reg {$reg} | egrep '^chr' |wc -l", false);
		$vars_filter = trim(implode("", $vars_filter));
		print "\t".number_format(100.0*$vars_filter/$vars, 2);
	}
		
	print "\n";		
}

?>