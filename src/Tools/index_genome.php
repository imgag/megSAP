<?php
/** 
	@page index_genome
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("index_genome", "Indexes a genome FASTA file with BWA and samtools.");
$parser->addInfile("in", "FASTA file to index.", false);
extract($parser->parse($argv));

$use_bwa1 = get_path("use_bwa1", true);
if ($use_bwa1==true || $use_bwa1=="true")
{
	exec2(get_path("bwa")." index -a bwtsw $in");
}
else
{
	exec2(get_path("bwa-mem2")." index $in");
}

exec2(get_path("samtools")." faidx $in");
exec2("md5sum -b $in > {$in}.md5");

?>
