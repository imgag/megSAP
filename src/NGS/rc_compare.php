<?php
/**
 * @page rc_compare
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("rc_compare", "Compare normalized read counts produced by analyze_rna.");
$parser->addInfile("in1", "Counts file for the tumor/case sample", false, true);
$parser->addInfile("in2", "Counts file for the reference/control samples", false, true);
$parser->addOutfile("out", "Output file containing the comparison.", false);
// optional parameters
$parser->addEnum("method", "The comparison method.", true, array("log1p_fc", "log_fc", "fc"), "log1p_fc");
extract($parser->parse($argv));

// counts per million mapped reads
function expr_diff($expr_case, $expr_ctrl, $method = "log_fc")
{
	switch ($method) {
		case "log1p_fc":
			return log($expr_case + 1, 2) - log($expr_ctrl + 1, 2);
			break;
		case "log_fc":
			return log($expr_case, 2) - log($expr_ctrl, 2);
			break;
		case "fc":
			return $expr_case / $expr_ctrl;
			break;
		default:
			trigger_error("Method '{$method}' not available!", E_USER_ERROR);
			break;
	}
}
$rc_case = Matrix::fromTSV($in1);
$rc_ctrl = Matrix::fromTSV($in2);

$rc_diff = new Matrix();
$rc_diff->setHeaders($rc_case->getHeaders());

if ($rc_case->rows() != $rc_ctrl->rows())
{
	trigger_error("Case and control read count files have different number of rows!", E_USER_ERROR);
}

for ($i = 0; $i < $rc_case->rows(); $i ++)
{
	$gene_id = $rc_case->get($i, 0);
	$gene_id_ctrl = $rc_case->get($i, 0);
	
	if ($gene_id !== $gene_id_ctrl)
	{
		trigger_error("Case and control read count files have different order or list of genes!", E_USER_ERROR);
	}
	
	$gene_symbol = $rc_case->get($i, 2);
	
	$diff = expr_diff($rc_case->get($i, 1), $rc_ctrl->get($i, 1), $method);
	$rc_diff->addRow(array($gene_id, $diff, $gene_symbol));
}

$rc_diff->toTSV($out);
