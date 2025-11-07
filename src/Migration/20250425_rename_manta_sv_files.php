<?php 

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("20250425_rename_manta_sv_files", "Rename old MANTA SV files to new name.");
$parser->addString("ps", "Processed sample name.", false);
$parser->addString("folder", "Processed sample folder.", false);
extract($parser->parse($argv));


//check if old files exist
$old_bedpe = "{$folder}/{$ps}_manta_var_structural.bedpe";
$old_vcf = "{$folder}/{$ps}_manta_var_structural.vcf.gz";
$old_tbi = "{$folder}/{$ps}_manta_var_structural.vcf.gz.tbi";
$old = [];
if (file_exists($old_bedpe)) $old[] = "bedpe";
if (file_exists($old_vcf)) $old[] = "vcf.gz";
if (file_exists($old_tbi)) $old[] = "vcf.gz.tbi";

//no old data > nothing to do
if (count($old)==0) exit(0);

//not all old files present > abort
if (count($old)!=3)
{
	print "{$ps} skipped - only part of old files present: ".implode("/", $old)."\n";
	exit(0);
}

//check if new files exist
$new_bedpe = "{$folder}/{$ps}_var_structural_variants.bedpe";
$new_vcf = "{$folder}/{$ps}_var_structural_variants.vcf.gz";
$new_tbi = "{$folder}/{$ps}_var_structural_variants.vcf.gz.tbi";
$new = [];
if (file_exists($new_bedpe)) $new[] = "bedpe";
if (file_exists($new_vcf)) $new[] = "vcf.gz";
if (file_exists($new_tbi)) $new[] = "vcf.gz.tbi";

//no new data available > move it
if (count($new)==0)
{
	print "{$ps}: renaming old to new\n";
	exec2("mv $old_bedpe $new_bedpe");
	exec2("mv $old_vcf $new_vcf");
	exec2("mv $old_tbi $new_tbi");
	exit(0);
}

//not all new files present > abort
if (count($new)!=3)
{
	print "{$ps} skipped - only part of new files present: ".implode("/", $old)." new: ".implode("/", $new)."\n";
	exit(0);
}

//all new files present > delete old
print "{$ps}: deleting old\n";
exec2("rm -rf {$folder}/{$ps}_manta_var_structural.*");

?>