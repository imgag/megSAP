<?php

include("/mnt/storage2/megSAP/pipeline/src/Common/all.php");

$type = $argv[1];

$db = DB::getInstance("NGSD");

print "##fileformat=VCFv4.1\n";
print "##INFO=<ID=BLA,Number=.,Type=String,Description=\"bla\">\n";
print "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";

if ($type=="causal")
{
	$cond = "rcv.causal='1'";
}
else
{
	$cond = "rcv.exclude_artefact='1'";
}
$res = $db->executeQuery("SELECT v.chr, v.start, v.end, v.ref, v.obs FROM variant v, report_configuration_variant rcv WHERE v.id=rcv.variant_id AND {$cond}");
foreach($res as $row)
{
	$chr = $row['chr'];
	$pos = $row['start'];
	$ref = $row['ref'];
	$obs = $row['obs'];
	if ($ref=='-') //insertion
	{
		$base = get_ref_seq("GRCh38", $chr, $pos, $pos);
		$ref = $base;
		$obs = $base.$obs;
	}
	else if ($obs=='-') //deletion
	{
		$pos -= 1;
		$base = get_ref_seq("GRCh38", $chr, $pos, $pos);
		$obs = $base;
		$ref = $base.$ref;
	}
	print "$chr\t$pos\t.\t$ref\t$obs\t30\tPASS\tBLA=bla\n";
}
?>