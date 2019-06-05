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
	foreach ($values as $key => $value)
	{
		fwrite($handle, implode("\t", [ $key, $value[0], $value[1] ]) . "\n");
	}
	fclose($handle);
}

rewrite_qcml($argv[1], $tmp1);
rewrite_qcml($argv[2], $tmp2);

//diff tsv outputs
if (isset($argv[3]))
{
	$diff_command = "numdiff --absolute-tolerance ".$argv[3];
}
else
{
	$diff_command = "diff -b";
}
list($stdout, $stderr, $code) = exec2("$diff_command $tmp1 $tmp2", false);
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