<?php
/**
	@page vc_strelka
	@todo use all.vcf from Strelka (currently only somatic variants are returned)
	@todo sparate SNV and INDEL vcfs, might be better option to pass vcf integrity test; vcf2gsvar - use snv and indel file for conversion
	@todo fix header (add tool name to filter name)
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// add parameter for command line ${input1.metadata.bam_index}
// parse command line arguments
$parser = new ToolBase("vc_strelka", "Call somatic variants with strelka. Creates an VCF file.");
$parser->addInfile("t_bam",  "Tumor BAM file. Remember: Place bai-file in the same folder as *.bam-file. Name it *.bam.bai.", false);
$parser->addInfile("n_bam",  "Normal BAM format. Remember: Place bai-file in the same folder as *.bam-file. Name it *.bam.bai.", false);
$parser->addOutfile("out", "Output file in VCF format (gzipped and tabix indexed).", false);
//optional
$parser->addFlag("k", "Keep all variants. Otherwise all variants that do not pass all filters will be removed.");
$parser->addFlag("amplicon",  "Enables amplicon mode.");
$parser->addString("build", "The genome build to use.", true, "GRCh37");
$parser->addInfile("config", "Config file for strelka.", true, true);
$parser->addInt("threads", "Number of threads (make jobs) to use.", true, 4);
extract($parser->parse($argv));

//get processed sample names
$t_ps = basename($t_bam, ".bam");
$n_ps = basename($n_bam, ".bam");
$t_bam = realpath($t_bam);
$n_bam = realpath($n_bam);

//get config file
if (!isset($config))
{
	if ($amplicon)
	{
		$config = get_path('strelka')."/../etc/strelka_config_bwa_amplicon.ini";
	}
	else
	{
		$config = get_path('strelka')."/../etc/strelka_config_bwa.ini";
	}
}

//run strelka
$strelka_folder = $parser->tempFolder()."/strelkaAnalysis";
$parser->exec(get_path('strelka')."/configureStrelkaWorkflow.pl",
	"--tumor {$t_bam} --normal {$n_bam} --ref ".get_path("local_data")."/{$build}.fa --output-dir={$strelka_folder} --config={$config}",
	true);
$parser->exec("make", "-C {$strelka_folder} -j {$threads}", true);

//merge vcf files
$vcf_combined = $parser->tempFile("_combined.vcf");	//$out."_test";
$file1 = Matrix::fromTSV("$strelka_folder/results/all.somatic.snvs.vcf");
$file2 = Matrix::fromTSV("$strelka_folder/results/all.somatic.indels.vcf");
$filec = new Matrix();

//set headers
$tmp_headers = $file1->getHeaders();
$t_index = -1;
$n_index = -1;
for($i=0;$i<count($tmp_headers);++$i)
{
	$value = $tmp_headers[$i];
	if($value == "TUMOR")
	{
		$tmp_headers[$i] = $t_ps;	//replace TUMOR header by tumor ID
		$idx_tumor = $i;
	}
	if($value == "NORMAL")
	{
		$tmp_headers[$i] = $n_ps;	//replace NORMAL header by normal ID
		$idx_normal = $i;
	}
}
$filec->setHeaders($tmp_headers);
$filec->setComments($file2->getComments());	//from indel file
foreach($file1->getComments() as $comment)	//from snv file
{
	//filter identical comments
    foreach($filec->getComments() as $c)
    {
		if($comment == $c) continue;
    }
    $filec->addComment($comment);
}
//set comments
$tmp_comments = array();
$tmp_filters = array();
foreach($filec->getComments() as $c)
{
	if(strpos($c,"#FILTER")===0)
	{
		$id = NULL;
		$desc = NULL;
		list(,$m) = explode("=",$c,2);
		$m = trim($m,'><');
		$fs = explode(",",$m);
		foreach($fs as $f)
		{
			list($n,$v) = explode("=",$f);
			if($n=="ID")	$id = $v;
			if($n=="Description")	$desc = trim($v,"\"");
		}

		if(is_null($id) || is_null($desc))	trigger_error("Could not identify filter in line '$c'!",E_USER_ERROR);
		if(isset($tmp_filters[$id]) && $tmp_filters[$id]!=$desc)
		{
			$desc = $tmp_filters[$id]." OR ".$desc;
			trigger_error("Duplicate filter ID '$id' with different descriptions '$desc' != '".$tmp_filters[$id]."'!",E_USER_WARNING);
		}
		$tmp_filters[$id] = $desc;
		continue;
	}
	//##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth for tier1 (used+filtered)">
	if(strpos($c,"#FORMAT")===0)
	{
		$id = NULL;
		$desc = NULL;
		list(,$m) = explode("=",$c,2);
		$m = trim($m,'><');
		$fs = explode(",",$m,4);
		foreach($fs as $f)
		{
			list($n,$v) = explode("=",$f);
			if($n=="ID")	$id = $v;
			if($n=="Type")	$desc = trim($v,"\"");
		}
		if($id=="DP")
		{
			if($desc=="Read depth for tier1 (used+filtered)")	$desc = "Read depth for tier1";
			$c = "#FORMAT=<ID=$id,Number=1,Type=Integer,Description=\"$desc\">";
		}
	}
	$tmp_comments[] = $c;
}
//if variant list is reduced add basic filter information. Example from freebayes: ##filter="QUAL > 5 & AO > 2"
if(!$k)	$tmp_comments[] = "#filter=\"".implode(" & ",array_keys($tmp_filters))."\"";
else
{
	foreach($tmp_filters as $id => $desc)
	{
		$tmp_comments[] = "#FILTER=<ID=s.$id,Description=\"$desc (strelka).\">";
	}
}
$tmp_comments[] = "#FILTER=<ID=special_chromosome,Description=\"Special chromosome.\">";
$tmp_comments[] = "#PEDIGREE=<Tumor=$t_ps,Normal=$n_ps>";	//add pedigree information for SNPeff
$tmp_comments = array_unique($tmp_comments);	//filter duplicate vcf comments
$tmp_comments = sort_vcf_comments($tmp_comments);	//sort vcf comments
$filec->setComments($tmp_comments);
//combine variants
for($i=0; $i<$file1->rows();++$i)
{
	$row = $file1->getRow($i);
	$filter = explode(";",$row[6]);
	foreach($filter as $j => $f)	//prefix filter with 's.'
	{
		if(empty($f) || $f == "." || $f == "PASS")	continue;
		$filter[$j] = "s.".$f;
	}
	$row[6] = implode(";",$filter);
	if(chr_check($row[0], 22, false) === FALSE)
	{
		if($row[6] == "PASS")	$row[6] = "";
		$row[6] .= ";s.special_chromosome"; //skip bad chromosomes
	}
	$row[6] = trim($row[6],';');
	$filec->addRow($row);
}
for($i=0; $i<$file2->rows();++$i)
{
	//remove overfluent '.' at end or beginning of indels that can be found with long indels - bug?
	$tmp = $file2->getRow($i);
	$tmp[3] = trim($tmp[3], ".");
	$tmp[4] = trim($tmp[4], ".");
	$filter = explode(";",$tmp[6]);
	foreach($filter as $j => $f)	//prefix filter with 's.'
	{
		if(empty($f) || $f == "." || $f == "PASS")	continue;
		$filter[$j] = "s.".$f;
	}
	$tmp[6] = implode(";",$filter);
	if(chr_check($tmp[0], 22, false) === FALSE)
	{
		if($tmp[6] == "PASS")	$row[6] = "";
		$tmp[6] .= ";special_chromosome"; //skip bad chromosomes
	}
	$tmp[6] = trim($tmp[6],';');
	$filec->addRow($tmp);
}
if(!$k)
{
	$tmp_filec = new Matrix();
	$tmp_filec->setHeaders($filec->getHeaders());
	$tmp_filec->setComments($filec->getComments());
	for($i=0;$i<$filec->rows();++$i)
	{
		$row = $filec->getRow($i);
		if($row[6]!="PASS")	continue;
		$tmp_filec->addRow($row);
	}
	$filec = $tmp_filec;
}
$filec->toTSV($vcf_combined);

//align INDELs to the left (this improves the ability of varFilter to remove overlapping indels)
$vcf_aligned = $parser->tempFile("_aligned.vcf");
$parser->exec(get_path("ngs-bits")."VcfLeftNormalize"," -in {$vcf_combined} -out {$vcf_aligned} -ref ".get_path("local_data")."/{$build}.fa", true);

//sort variants
$vcf_sorted = $parser->tempFile("_sorted.vcf");
//Nb: VcfStreamSort cannot be used since two vcf files (variants, indels) are concatenated and variants are not grouped by chromosome; total number of variants should be low (somatic).
$parser->exec(get_path("ngs-bits")."VcfSort","-in {$vcf_aligned} -out {$vcf_sorted}", true);

// filter high quality variants
$vcf_filtered1 = $strelka_folder."/strelka_filtered1.vcf";
$parser->execTool("NGS/filter_vcf.php", "-in $vcf_sorted -out $vcf_filtered1 -type somatic-lq -keep", true);

// flag off-target variants
$vcf_filtered2 = $vcf_filtered1;
if(!empty($target))
{
	$vcf_filtered2 = $strelka_folder."/strelka_filtered2.vcf";
	$parser->exec(get_path("ngs-bits")."VariantFilterRegions", "-in $vcf_filtered1 -mark off-target -reg $target -out $vcf_filtered2", true);
}

//zip and index output file
$parser->exec("bgzip", "-c $vcf_filtered2 > $out", true);
$parser->exec("tabix", "-p vcf $out", true);

?>