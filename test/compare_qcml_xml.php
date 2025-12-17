<?php

/*
 * Compare two qcML files
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/all.php");

// Rewrites quality parameters from qcML file into tabular format
function rewrite_qcml($qcml, &$tmpfile)
{
	// load qcML with SimpleXML
	$sxml = simplexml_load_file($qcml);
	if ($sxml===false) trigger_error("Could not load XML file: $qcml", E_USER_ERROR);

	// array to hold entries as accession => [name, value]
	$values = [];
	// output file
	$tmpfile = temp_file("", "rewrite_qcml");
	// read values into array
	foreach ($sxml->runQuality->qualityParameter as $qp)
	{
		$values[(string) $qp["accession"]] = [
			(string) $qp["name"],
			(string) $qp["value"]
		];
	}
	// sort array by key
	ksort($values);
	// write values as tab-separated (accession, name, value)
	$handle = fopen2($tmpfile, "w");
	fwrite($handle, implode("\t", [ "#id", "name", "value" ]) . "\n");
	foreach ($values as $key => $value)
	{
		fwrite($handle, implode("\t", [ $key, $value[0], $value[1] ]) . "\n");
	}
	fclose($handle);
}

rewrite_qcml($argv[1], $tmp1);
rewrite_qcml($argv[2], $tmp2);

//diff tsv outputs
$args = [];
$args[] =  "-in1 {$tmp1}";
$args[] =  "-in2 {$tmp2}";
if (isset($argv[3]))
{
	$args[] = "-diff_abs value=".$argv[3];
}
list($stdout, $stderr, $code) = execApptainer("ngs-bits", "TsvDiff", implode(" ", $args), [$tmp1, $tmp2], [], false, false, false);
foreach ($stdout as $line)
{
	if ($line=="") continue;
	print $line."\n";
}
foreach ($stderr as $line)
{
	if ($line=="") continue;
	print $line."\n";
}

exit($code);

?>