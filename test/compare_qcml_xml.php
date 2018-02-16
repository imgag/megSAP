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
	// array to hold entries as accession => [name, value]
	$values = [];
	// output file
	$tmpfile = temp_file();
	// read values into array
	foreach ($sxml->runQuality->qualityParameter as $qp)
	{
		$values[(string) $qp["accession"]] = [
			(string) $qp["name"],
			(string) $qp["value"]
		];
	}
	// sorty array by key
	ksort($values);
	// write values as tab-separated (accession, name, value)
	$handle = fopen($tmpfile, "w");
	foreach ($values as $key => $value)
	{
		fwrite($handle, implode("\t", [ $key, $value[0], $value[1] ]) . "\n");
	}
	fclose($handle);
}

rewrite_qcml($argv[1], $tmp1);
rewrite_qcml($argv[2], $tmp2);

// diff tsv outputs
list($stdout, $stderr, $code) = exec2("diff -b $tmp1 $tmp2", false);
foreach ($stdout as $line)
{
	print $line."\n";
}
foreach ($stderr as $line)
{
	print $line."\n";
}

exit($code);

?>