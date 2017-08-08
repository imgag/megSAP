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
$parser->addString("column", "The value column to base the comparison on.", true, "fpkm");
extract($parser->parse($argv));

// calculate expression difference
function expr_diff($expr_case, $expr_ctrl, $method)
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

if ($rc_case->rows() != $rc_ctrl->rows())
{
	trigger_error("Case and control read count files have different number of rows!", E_USER_ERROR);
}

$rc_case_id_idx = $rc_case->getColumnIndex("gene_id");
$rc_ctrl_id_idx = $rc_ctrl->getColumnIndex("gene_id");

$rc_case_value_idx = $rc_case->getColumnIndex($column);
$rc_ctrl_value_idx = $rc_ctrl->getColumnIndex($column);

$rc_diff = new Matrix();
$rc_diff->setHeaders(array("gene_id", "case", "ctrl", "difference"));
$rc_diff->addComment("TUMOR=".$in1);
$rc_diff->addComment("NORMAL=".$in2);
$rc_diff->addComment("METHOD=".$method);
$rc_diff->addComment("VALUE=".$column);

for ($i = 0; $i < $rc_case->rows(); $i ++)
{
	$gene_id = $rc_case->get($i, $rc_case_id_idx);
	$gene_id_ctrl = $rc_case->get($i, $rc_ctrl_id_idx);
	
	if ($gene_id !== $gene_id_ctrl)
	{
		trigger_error("Case and control read count files have different order or list of genes!", E_USER_ERROR);
	}
	
	$expr_case = $rc_case->get($i, $rc_case_value_idx);
	$expr_ctrl = $rc_ctrl->get($i, $rc_ctrl_value_idx);
	$diff = expr_diff($expr_case, $expr_ctrl, $method);
	$rc_diff->addRow(array($gene_id, $expr_case, $expr_ctrl, $diff));
}

$rc_diff->toTSV($out);
