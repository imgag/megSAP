<?php
/** 
	@page converter_tsv2upd 
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

//parse command line arguments
$parser = new ToolBase("converter_tsv2upd", "");
$parser->addString("f",  "Analyzed folder (father).", false);
$parser->addString("m",  "Analyzed folder (mother).", false);
$parser->addString("c",  "Analyzed folder (child).", false);
$parser->addOutfile("out",  "Output file for UPDtool.", false, "tsv");
extract($parser->parse($argv));

//get files
$tmp_tsv = glob($f."/*.GSvar");
$tmp_bam = glob($f."/*.bam");
if(count($tmp_tsv) != 1 || count($tmp_bam) != 1) trigger_error ("Could not determine father data.", E_USER_ERROR);
$f_tsv = $tmp_tsv[0];
$f_bam = $tmp_bam[0];
$parser->log("Using '".$f_tsv."' & '".$f_bam."' as father.");

$tmp_tsv = glob($m."/*.GSvar");
$tmp_bam = glob($m."/*.bam");
if(count($tmp_tsv) != 1 || count($tmp_bam) != 1) trigger_error ("Could not determine mother data.", E_USER_ERROR);
$m_tsv = $tmp_tsv[0];
$m_bam = $tmp_bam[0];
$parser->log("Using '".$m_tsv."' & '".$m_bam."' as mother.");

$tmp_tsv = glob($c."/*.GSvar");
$tmp_bam = glob($c."/*.bam");
if(count($tmp_tsv) != 1 || count($tmp_bam) != 1) trigger_error ("Could not determine child data.", E_USER_ERROR);
$c_tsv = $tmp_tsv[0];
$c_bam = $tmp_bam[0];
$parser->log("Using '".$c_tsv."' & '".$c_bam."' as child.");

//SampleOverview => combine trio in one file
$so = $parser->tempFile("_annovar.tsv");
$parser->exec(get_path("ngs-bits")."SampleOverview", "-in $f_tsv $m_tsv $c_tsv -out $so", true);	

//VariantAnnotateFrequency => add depth and variant frequency
$vaf_options = " -depth -ref ".get_path("local_data")."/hg19.fa";
$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency", "-in $so -bam $f_bam -out $so -name father $vaf_options", true);	
$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency", "-in $so -bam $m_bam -out $so -name mother $vaf_options", true);	
$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency", "-in $so -bam $c_bam -out $so -name child $vaf_options", true);	

//filter all variants and generate output
$combined = Matrix::fromTSV($so);
$idx_mg = $combined->getColumnIndex(basename($m_tsv));
$idx_md = $combined->getColumnIndex("mother_depth");
$idx_fg = $combined->getColumnIndex(basename($f_tsv));
$idx_fd = $combined->getColumnIndex("father_depth");
$idx_cg = $combined->getColumnIndex(basename($c_tsv));
$idx_cd = $combined->getColumnIndex("child_depth");
$output = new Matrix();
for($i=0; $i<$combined->rows(); ++$i)
{
    $row = $combined->getRow($i);
    
    //skip indels
    if($row[5] == "INDEL") continue;
    
    //check that all are covered at 10x
    if($row[$idx_md] < 10 || $row[$idx_fd] <10 || $row[$idx_cd] < 10) continue;
    
    //mother
    $g_m = "AA";
    if($row[$idx_mg] == "yes (het)")
    {
	    $g_m = "AB";
    }
    if($row[$idx_mg] == "yes (hom)")
    {
	    $g_m = "BB";
    }
    
    //father
    $g_f = "AA";
    if($row[$idx_fg] == "yes (het)")
    {
	    $g_f = "AB";
    }
    if($row[$idx_fg] == "yes (hom)")
    {
	    $g_f = "BB";
    }
    
    //child
    $g_c = "AA";
    if($row[$idx_cg] == "yes (het)")
    {
	    $g_c = "AB";
    }
    if($row[$idx_cg] == "yes (hom)")
    {
	    $g_c = "BB";
    }

    $output->addRow(array($row[0], $row[1], $g_f, $g_m, $g_c));
}

//output
$output->toTSV($out);
?>
