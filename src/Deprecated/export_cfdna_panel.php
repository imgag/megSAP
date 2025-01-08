<?php
/** 
	@page export_cfdna_panel
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("export_cfdna_panel", "Exports the cfDNA panel VCF/BED file for a given processed sample.");
$parser->addString("name", "Processed sample name of the cfDNA sample.", false);
$parser->addOutfile("out", "Output file (in VCF or BED format).", false);
$formats = array("VCF", "BED");
$parser->addString("format", "Output format. Possible formats: ".implode(", ", $formats), true, "VCF");
$parser->addEnum("db",  "Database to connect to.", true, db_names(), "NGSD");
extract($parser->parse($argv));

//database
if (db_is_enabled("NGSD"))
{
	$db = DB::getInstance("NGSD", false);
}
else
{
	trigger_error("NGSD connection is required to determine tumor id!", E_USER_ERROR);
}

if (!isset($db)) trigger_error("NGSD connection is required to determine tumor id!", E_USER_ERROR);

//related samples
list($sample_name, $ps_num) = explode("_", $name);
$res_samples = $db->executeQuery("SELECT id, name FROM sample WHERE name=:name", ["name" => $sample_name]);
if (count($res_samples) !== 1)
{
	trigger_error("Could not find sample for processed sample {$name}!", E_USER_WARNING);
}
$sample_id = $res_samples[0]['id'];
$system = null;
$sys = load_system($system, $name);
$res = $db->executeQuery("SELECT * FROM sample_relations WHERE relation='tumor-cfDNA' AND (sample1_id=:sid OR sample2_id=:sid)", array("sid" => $sample_id));

$psamples = [];
$ps_ids = [];
foreach ($res as $row)
{
	$sample_id_annotation = $row['sample1_id'] != $sample_id ? $row['sample1_id'] : $row['sample2_id'];
	$res = $db->executeQuery("SELECT ps.sample_id, ps.process_id, ps.processing_system_id, ps.quality, sys.id, sys.type, CONCAT(s.name, '_', LPAD(ps.process_id, 2, '0')) as psample, ps.id FROM processed_sample as ps, processing_system as sys, sample as s WHERE ps.sample_id=:sid AND ps.processing_system_id=sys.id AND ps.sample_id=s.id AND sys.type!='cfDNA (patient-specific)' AND (NOT ps.quality='bad') ORDER BY ps.process_id ASC", array("sid" => $sample_id_annotation));
	$psamples = array_merge($psamples, array_column($res, 'psample'));
}
if (count($psamples) > 1)
{
	trigger_error("Found more than one referenced tumor, using first one: " . implode(" ", $psamples), E_USER_WARNING);
}
if (count($psamples) === 0)
{
	trigger_error("Could not find any related tumor processed sample!", E_USER_ERROR);
}

$tumor_id = $psamples[0];
trigger_error("Tumor id extracted from the NGSD: ${tumor_id}", E_USER_NOTICE);

// get patient-specific target region/SNVs from NGSD

$tumor_ps_id = get_processed_sample_id($db, $tumor_id);
$sys_id = $db->getValue("SELECT id FROM processing_system WHERE name_short='".$sys["name_short"]."'");

if ($format == "VCF")
{
	$vcf_content = $db->getValue("SELECT vcf FROM cfdna_panels WHERE tumor_id=${tumor_ps_id} AND processing_system_id=${sys_id}");
	file_put_contents($out, $vcf_content);
}
else if ($format == "BED")
{
	$bed_content = $db->getValue("SELECT bed FROM cfdna_panels WHERE tumor_id=${tumor_ps_id} AND processing_system_id=${sys_id}");
	file_put_contents($out, $bed_content);
}
else
{
	trigger_error("Invalid output format! Possible formats: ".implode(", ", $formats), E_USER_ERROR);
}


?>
