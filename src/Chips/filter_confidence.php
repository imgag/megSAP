<?php

/**
	@page filter_confidence
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("filter_confidence", "\$Rev: 2 $", "Filters an Affymetrics 6.0 chip export according to conficence.");
$parser->addInfile("in",  "SNP array genotype file.", false);
$parser->addFloat("thres", "Threshold over which the genotypes are set to NoCall.", false);
$parser->addOutfile("out", "Filtered genotype file.", false);

extract($parser->parse($argv));

$data = Matrix::fromTSV($in);
$rows = $data->rows();
$cols = $data->cols();
for ($r=1; $r<$rows; ++$r)
{
	for ($c=2; $c<$cols; $c+=2)
	{
		$confidence = $data->get($r, $c);
		if ($confidence<=1 && $confidence>$thres)
		{
			$data->set($r, $c-1, "NoCall");
		}
	}
}
$data->toTSV($out);

?>
