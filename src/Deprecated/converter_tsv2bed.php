<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("converter_tsv2bed", "Generate bed-files for ds-experiment.");
$parser->addInfile("tsv",  "Input file in tsv-format.", false);
$parser->addOutfile("out",  "Output file in bed-format.", false);
extract($parser->parse($argv));

//read tsv-file
$input = Matrix::fromTSV($tsv);
$output = array();
for($i = 0; $i < $input->rows(); ++$i)
{
	$row = $input->getRow($i);
	$output[] = $row[0]."\t".($row[1]-1)."\t".$row[2]."\n";
}
file_put_contents($out, $output);
