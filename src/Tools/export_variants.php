<?php
/** 
	@page export_variants 
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("export_variants", "Exports variants for a set of samples and genes from NGSD.");
$parser->addOutfile("out", "Output VCF file.", false);
$parser->addInfile("samples", "Processed sample names file (one name per line).", false);
$parser->addInfile("genes", "Gene names file (one name per line).", false);
$parser->addFlag("anno", "If set, produces an annotated version of the output VCF.");
$parser->addFlag("gsvar", "If set, produces the GSvar file from the annotated VCF.");
$parser->addEnum("db",  "Database to export from to.", true, db_names(), "NGSD");
$parser->addFlag("vcf_mode",  "The 'samples' file contains VCF.GZ file names instead of sample names.");
extract($parser->parse($argv));

//init
$ngsbits = get_path("ngs-bits");

//determine ROI from genes
print "determining ROI from genes...\n";
$bed = $parser->tempFile(".bed", "export_variants");
$pipeline = [
	[$ngsbits."GenesToBed", "-in {$genes} -source ensembl -mode gene"],
	[$ngsbits."BedExtend", "-n 5000"],
	[$ngsbits."BedMerge", "-out $bed"]
];
list($stdout, $stderr) = $parser->execPipeline($pipeline, "genes to bed");
foreach(array_merge($stdout, $stderr) as $line)
{
	$line = trim($line);
	if ($line!="") print "  Notice: $line\n";
}

//determine input VCFs from processed samples
$vcfs = array();
if ($vcf_mode)
{
	$tmp = file($samples);
	for($i=0; $i<count($tmp); ++$i)
	{
		$vcf = trim($tmp[$i]);
		if ($vcf=="") continue;
		
		if (!file_exists($vcf))
		{
			print "Notice: Skipping VCF file because it does not exist: {$vcf}\n";
			continue;
		}
		
		$vcfs[$i] = $vcf;
	}
}
else
{
	print "\n";
	print "determining input VCFs for samples...\n";
	$db = DB::getInstance($db);
	$samples = file($samples);
	foreach($samples as $ps)
	{
		$ps = trim($ps);
		if ($ps=="" || $ps[0]=="#") continue;
		
		$info = get_processed_sample_info($db, $ps, false);
		if (is_null($info))
		{
			print "Notice: Skipping processed sample {$ps} - not found in NGSD!\n";
			continue;
		}
		
		$vcf = $info['ps_folder']."/".$info['ps_name']."_var.vcf.gz";
		if (!file_exists($vcf))
		{
			print "Notice: Skipping processed sample {$ps} - VCF not found: {$vcf}\n";
			continue;
		}
		
		$vcfs[$ps] = $vcf;
	}
}


//merge variants from VCF files
print "\n";
print "merging variants from ".count($vcfs)." VCFs...\n";
$variants = array();
foreach($vcfs as $ps => $vcf)
{
	print "  $vcf\n";
	list($stdout) = $parser->exec("tabix", "{$vcf} -R {$bed}", false);
	foreach($stdout as $line)
	{
		$line = trim($line);
		if($line=="") continue;
		
		$parts = explode("\t", $line);
		$parts[2] = "."; //ID
		$parts[5] = "100"; //QUAL
		$parts[6] = "."; //FILTER
		$parts[7] = "."; //INFO
		$parts[8] = "GT:DP:AO"; //FORMAT
		$parts[9] = "0/1:100:50"; //SAMPLE
		
		$variant = implode("\t", $parts);
		
		$variants[$variant] = true;
	}
}

//write output
print "\n";
print "writing ".count($variants)." distinct variants to VCF file...\n";
$sample_name = basename($out, ".vcf");
$h = fopen2($out, "w");
fputs($h, "##fileformat=VCFv4.2\n");
fputs($h, "##ANALYSISTYPE=GERMLINE_SINGLESAMPLE\n");
fputs($h, "##SAMPLE=<ID={$sample_name},Gender=n/a,IsTumor=no,IsFFPE=no,DiseaseGroup=n/a,DiseaseStatus=affected>\n");
fputs($h, "##FILTER=<ID=off-target,Description=\"Variant marked as 'off-target'\">\n");
fputs($h, "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	{$sample_name}\n");
foreach($variants as $variant => $dummy)
{
	fputs($h, $variant."\n");
}
fclose($h);

//sort
$parser->exec($ngsbits."VcfSort", "-in $out -out $out");

//annotate
$vcf_anno = substr($out, 0, 4)."_anno.vcf";
if ($anno)
{
	print "\n";
	print "annotating VCF file...\n";
	$parser->execTool("NGS/an_vep.php", "-in {$out} -out {$vcf_anno} -threads 5");
}

//convert to GSvar
if ($anno && $gsvar)
{
	print "\n";
	print "generating GSvar file...\n";
	$gsvar = substr($out, 0, 4).".GSvar";
	$parser->execTool("NGS/vcf2gsvar.php", "-in {$vcf_anno} -out {$gsvar}");
}
?>