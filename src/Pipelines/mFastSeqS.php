<?php

/**
	@page mFastSeqS
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("mFastSeqS", "mFastSeqS pipeline.");
$parser->addString("folder", "Sample data folder.", false);
$parser->addString("name", "Base file name, typically the processed sample ID (e.g. 'GS120001_01').", false);
$steps_all = array("ma", "rc", "zc");
$parser->addString("steps", "Comma-separated list of steps to perform:\nma=mapping, rc=read counting, zc=z-score calculation.", true, implode(",", $steps_all));
extract($parser->parse($argv));

//check steps
$steps = explode(",", $steps);
foreach($steps as $step)
{
	if (!in_array($step, $steps_all)) trigger_error("Unknown processing step '$step'!", E_USER_ERROR);
}

//arms to skip because normalized read count is too low (<0.2%)
$skip = array("chr13p", "chr14p", "chr15p", "chr21p", "chr22p", "chrYp", "chrYq");	

//chromosome arm coordinates
$arms = array(
	"chr1p" => "chr1:0-125000000",
	"chr1q" => "chr1:125000000-249250621",
	"chr2p" => "chr2:0-93300000",
	"chr2q" => "chr2:93300000-243199373",
	"chr3p" => "chr3:0-91000000",
	"chr3q" => "chr3:91000000-198022430",
	"chr4p" => "chr4:0-50400000",
	"chr4q" => "chr4:50400000-191154276",
	"chr5p" => "chr5:0-48400000",
	"chr5q" => "chr5:48400000-180915260",
	"chr6p" => "chr6:0-61000000",
	"chr6q" => "chr6:61000000-171115067",
	"chr7p" => "chr7:0-59900000",
	"chr7q" => "chr7:59900000-159138663",
	"chr8p" => "chr8:0-45600000",
	"chr8q" => "chr8:45600000-146364022",
	"chr9p" => "chr9:0-49000000",
	"chr9q" => "chr9:49000000-141213431",
	"chr10p" => "chr10:0-40200000",
	"chr10q" => "chr10:40200000-135534747",
	"chr11p" => "chr11:0-53700000",
	"chr11q" => "chr11:53700000-135006516",
	"chr12p" => "chr12:0-35800000",
	"chr12q" => "chr12:35800000-133851895",
	"chr13p" => "chr13:0-17900000",
	"chr13q" => "chr13:17900000-115169878",
	"chr14p" => "chr14:0-17600000",
	"chr14q" => "chr14:17600000-107349540",
	"chr15p" => "chr15:0-19000000",
	"chr15q" => "chr15:19000000-102531392",
	"chr16p" => "chr16:0-36600000",
	"chr16q" => "chr16:36600000-90354753",
	"chr17p" => "chr17:0-24000000",
	"chr17q" => "chr17:24000000-81195210",
	"chr18p" => "chr18:0-17200000",
	"chr18q" => "chr18:17200000-78077248",
	"chr19p" => "chr19:0-26500000",
	"chr19q" => "chr19:26500000-59128983",
	"chr20p" => "chr20:0-27500000",
	"chr20q" => "chr20:27500000-63025520",
	"chr21p" => "chr21:0-13200000",
	"chr21q" => "chr21:13200000-48129895",
	"chr22p" => "chr22:0-14700000",
	"chr22q" => "chr22:14700000-51304566",
	"chrXp" => "chrX:0-60600000",
	"chrXq" => "chrX:60600000-155270560",
	"chrYp" => "chrY:0-12500000",
	"chrYq" => "chrY:12500000-59373566"
	);

/* code for creating stats
	$file = file("ref_data.tsv");
	$square_sums = array();
	print "\$stats = array(\n";
	foreach($file as $line)
	{
		$line = trim($line);
		if ($line=="" || $line[0]=="#") continue;
		
		$parts = explode("\t", $line);
		$arm = $parts[0];
		$values = array_slice($parts, 1);
		$mean = mean($values);
		$stdev = stdev($values, $mean);
		if (!in_array($arm, $skip))
		{
			for($i=0; $i<count($values); ++$i)
			{
				$z = ($values[$i]-$mean)/$stdev;
				@$square_sums[$i] += $z*$z;
			}
		}
		print "\t\"$arm\" => array(".number_format($mean,4).", ".number_format($stdev,4)."),\n";
	}
	$mean = mean($square_sums);
	$stdev = stdev($square_sums, $mean);
	print_r($square_sums);
	print "\t\"genome-wide\" => array(".number_format($mean,4).", ".number_format($stdev,4).")\n";
	print ");\n";
	die;
*/
$stats = array(
        "chr1p" => array(3.2349, 0.0645),
        "chr1q" => array(3.5428, 0.0899),
        "chr2p" => array(2.9624, 0.0683),
        "chr2q" => array(5.6349, 0.1317),
        "chr3p" => array(3.0651, 0.0591),
        "chr3q" => array(4.5280, 0.0990),
        "chr4p" => array(1.8703, 0.0761),
        "chr4q" => array(6.3310, 0.2354),
        "chr5p" => array(2.0458, 0.0996),
        "chr5q" => array(5.2965, 0.1297),
        "chr6p" => array(1.9309, 0.0555),
        "chr6q" => array(4.7302, 0.1074),
        "chr7p" => array(1.8760, 0.0461),
        "chr7q" => array(3.6107, 0.0970),
        "chr8p" => array(1.4167, 0.0333),
        "chr8q" => array(4.0481, 0.1386),
        "chr9p" => array(1.6605, 0.0425),
        "chr9q" => array(2.0224, 0.0515),
        "chr10p" => array(1.4337, 0.0373),
        "chr10q" => array(2.9959, 0.0574),
        "chr11p" => array(1.9585, 0.0661),
        "chr11q" => array(2.5985, 0.1546),
        "chr12p" => array(1.4189, 0.1200),
        "chr12q" => array(3.5250, 0.0711),
        "chr13p" => array(0.0000, 0.0000),
        "chr13q" => array(3.7513, 0.0762),
        "chr14p" => array(0.0000, 0.0000),
        "chr14q" => array(2.7927, 0.0977),
        "chr15p" => array(0.0000, 0.0000),
        "chr15q" => array(2.1658, 0.0595),
        "chr16p" => array(0.9361, 0.0354),
        "chr16q" => array(1.1298, 0.0868),
        "chr17p" => array(0.6110, 0.0501),
        "chr17q" => array(1.1839, 0.0562),
        "chr18p" => array(0.4683, 0.0189),
        "chr18q" => array(2.1008, 0.1184),
        "chr19p" => array(0.3091, 0.0457),
        "chr19q" => array(0.6930, 0.0673),
        "chr20p" => array(0.8775, 0.0357),
        "chr20q" => array(0.5769, 0.0377),
        "chr21p" => array(0.0797, 0.0153),
        "chr21q" => array(1.4428, 0.0489),
        "chr22p" => array(0.0000, 0.0000),
        "chr22q" => array(0.5262, 0.0268),
        "chrXp" => array(1.9091, 0.5845),
        "chrXq" => array(4.5122, 1.4394),
        "chrYp" => array(0.0744, 0.0654),
        "chrYq" => array(0.1215, 0.1074),
		"genome-wide" => array(41.0000, 14.2488)
);


//(1) mapping 
if (in_array("ma", $steps))
{
	$parser->execTool("Pipelines/analyze.php", "-folder {$folder} -name {$name} -steps ma -no_abra");
}

//(2) read counting
$reads = "{$folder}/{$name}_counts.tsv";
if (in_array("rc", $steps))
{	
	//count reads per arm
	$counts = array();
	foreach($arms as $arm => $coords)
	{	
		list($stdout) = $parser->exec(get_path("samtools"), "view -q 15 {$folder}/{$name}.bam $coords | wc -l", false);
		$counts[$arm] = trim($stdout[0]);
	}
	$sum = array_sum(array_values($counts));

	//write count percentages
	$h = fopen($reads, "w");
	fwrite($h, "#arm\tread_count_percentage\n");
	foreach($arms as $arm => $coords)
	{
		fwrite($h, $arm."\t".number_format(100*$counts[$arm]/$sum, 3)."\n");
	}
	fwrite($h, "##reads_overall: {$sum}\n");
	fclose($h);
}

//(3) z-scores calculation
$zscores = "{$folder}/{$name}_zscores.tsv";
if (in_array("zc", $steps))
{	
	//read count percentages file
	$values = array();
	$file = file($reads);
	foreach($file as $line)
	{
		$line = trim($line);
		if ($line=="" || $line[0]=="#") continue;
		list($a, $v) = explode("\t", $line);
		$values[$a] = $v;
	}

	//write z-scores for arms (and genome-wide)
	$square_sum = 0;
	$h = fopen($zscores, "w");
	fwrite($h, "#arm\tz-score\n");
	foreach($arms as $arm => $coords)
	{
		if (in_array($arm, $skip))
		{
			$z = "n/a";
		}
		else
		{
			$z = number_format(($values[$arm]-$stats[$arm][0])/$stats[$arm][1], 2);
			$square_sum += $z*$z;
		}
		fwrite($h, "{$arm}\t{$z}\n");
	}
	fwrite($h, "genome-wide\t".number_format(($square_sum-$stats['genome-wide'][0])/$stats['genome-wide'][1], 2)."\n");
	fclose($h);
}

?>
