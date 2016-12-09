<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/all.php");

/* example:
NA12878_stats_fastq.qcML:    <qualityParameter ID="qp0001" name="read count" description="Total number of reads (forward and reverse reads of paired-end sequencing count as two reads)." value="1059810" cvRef="QC" accession="QC:2000005"/>
NA12878_stats_map.qcML:    <qualityParameter ID="qp0019" name="error estimation SNV percentage" description="SNV error percentage determined on special target region after mapping." value="0.1081" cvRef="QC" accession="QC:2000035"/>
NA12878_stats_vc.qcML:    <qualityParameter ID="qp0006" name="transition/transversion ratio" description="Transition/transversion ratio of SNV variants." value="3.32" cvRef="QC" accession="QC:2000018"/>
*/

function extract_values($filename)
{
	$output = array();
	$file = file($filename);
	foreach($file as $line)
	{
		list(, $line) = explode("<", $line);		
		$line = trim(substr($line, 16, -2));
		$parts =  explode("\"", $line);
		$output[] = $parts[11]."\t".$parts[3]."\t".$parts[7]."\n";
	}
	
	return $output;
}

//extract data
$data1 = extract_values($argv[1]);
$data2 = extract_values($argv[2]);

//store data to tmp files
$file1 = temp_file();
file_put_contents($file1, $data1);
$file2 = temp_file();
file_put_contents($file2, $data2);

//diff
list($stdout, $stderr, $code) = exec2("diff -b $file1 $file2", false);

//filter diff
$stdout2 = array();
foreach($stdout as $line)
{
	if ($line=="" || ($line[0]!=">" && $line[0]!="<")) continue;
	$stdout2[] = explode("\t", strtr($line, array("> "=>">\t", "< "=>"<\t")));
}

//sort diff
function custom_sort($line1, $line2)
{
	if ($line1[1]!=$line2[1]) return $line1[1]<$line2[1];
	return $line1[0]>$line2[0];
}
usort($stdout2, 'custom_sort');

//output
foreach($stdout2 as $line)
{
	print implode(" ", $line)."\n";
}

exit($code);
?>