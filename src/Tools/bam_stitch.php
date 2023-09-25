<?php
/** 
	@page bam_stitch
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("bam_stitch", "Stitches two BAM files together e.g to create demo data.");
$parser->addInfile("in1", "Main BAM file.", false);
$parser->addInfile("in2", "Secondary BAM file.", false);
$parser->addInfile("roi", "ROI to etract from 'in1'.", false);
$parser->addString("regs", "Comma-separeted list of regions in 'in1' to replace with data from 'in2'.", false);
$parser->addOutfile("out", "Output BAM file.", false);
extract($parser->parse($argv));

//store 'regs' as BED file
$regs = explode(",", $regs);
for($i=0; $i<count($regs); ++$i)
{
	$regs[$i] = strtr($regs[$i], array("-"=>"\t", ":"=>"\t", " "=>""));
}
$regs_file = $parser->tempFile(".bed");
file_put_contents($regs_file, implode("\n", $regs));

//extract sample names to correct readgroup info in SAM header
$ps_name1 = basename2($in1);
list($s_name1) = explode("_", $ps_name1."_");
$ps_name2 = basename2($in2);
$ps_name3 = basename2($out);

//extract data from 'in1'
print "Extracting data from 'in1'.\n";
$roi_sub = $parser->tempFile(".bed");
exec2(get_path("ngs-bits")."BedSubtract -in {$roi} -in2 {$regs_file} -out {$roi_sub}");
$bam_in1 = $parser->tempFile(".bam");
exec2(get_path("samtools")." view -h -L {$roi_sub} -M {$in1} | sed \"s/{$ps_name1}/{$ps_name3}/g\" | sed \"s/{$s_name1}/{$ps_name3}/g\" | ".get_path("samtools")." view -Sb > {$bam_in1}");

//extract data from 'in2'
print "Extracting data from 'in2'.\n";
$bam_in2 = $parser->tempFile(".bam");
exec2(get_path("samtools")." view -h -L {$regs_file} -M  {$in2} | sed \"s/{$ps_name2}/{$ps_name3}/g\" | ".get_path("samtools")." view -Sb > {$bam_in2}");

//merge bam files
print "Merging data.\n";
$bam_merge = $parser->tempFile(".bam");
exec2(get_path("samtools")." merge -f {$bam_merge} {$bam_in1} {$bam_in2}");

//sort and index output
print "Sorting and indexing 'out'.\n";
$parser->sortBam($bam_merge, $out, 4);
$parser->indexBam($out, 4);

?>

