<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");
error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("genlab_dna_rna", "Export DNA-RNA relations from GENLAB.");
extract($parser->parse($argv));

//check DB is enabled
if(!GenLabDB::isEnabled())
{
	trigger_error("GenLab database is not enabled - please add credentials to settings file!", E_USER_ERROR);
}

//get gender
$db = GenLabDB::getInstance();
$res = $db->executeQuery("SELECT * FROM v_ngs_dnarna");
foreach($res as $row)
{
	// var_dump($row);
	$dna_str = $row["LABORNUMMER"];
	if (!preg_match('/DX-?([0-9]{6})(_[0-9]{2})?/', $dna_str, $match_dna))
	// if (!preg_match('/DX([0-9]{6})/', $dna_str, $match_dna))
	{
		trigger_error("Could not extract DNA sample ID from string '{$dna_str}' - skipping!", E_USER_WARNING);
		// var_dump($row);
		continue;
	}	

	$rna_str = $row["T_UNTERSUCHUNG_1_MATERIALINFO"];
	if (!preg_match('/RNA-([0-9]{6})/', $rna_str, $match_rna))
	{
		trigger_error("Could not extract RNA sample ID from string '{$rna_str}' - skipping!", E_USER_WARNING);
		// var_dump($row);
		continue;
	}
	
	$dna_sample = "DX" . $match_dna[1];
	$rna_sample = "RX" . $match_rna[1];
	$line = implode("\t", [$dna_sample, "same sample", $rna_sample]);
	print($line . "\n");
}

?>