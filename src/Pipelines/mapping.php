<?php

/**
	@page mapping
	
	@todo Add pipeline test for MIPs/HaloPlex HS (only a few exons)
	
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("mapping", "Mapping pipeline.");
$parser->addInfileArray("in_for",  "Forward reads FASTQ file(s).", false);
$parser->addInfileArray("in_rev",  "Reverse reads FASTQ file(s).", false);
$parser->addString("out_folder", "Output folder.", false);
$parser->addString("out_name", "Output file base name (e.g. 'GS120001_01').", false);
//optional
$parser->addInfile("system",  "Processing system INI file (determined from 'out_name' by default).", true);
$parser->addInt("threads", "The maximum number of threads used.", true, 2);
$parser->addFlag("clip_overlap", "Soft-clip overlapping read pairs.", true);
$parser->addFlag("no_abra", "Skip realignment with ABRA.", true);
extract($parser->parse($argv));

//extract processing system information from DB
$sys = load_system($system, $out_name);

//detmerine sample sheet name
$sheet = dirname($in_for[0])."/SampleSheet.csv";

// determine output file base name
$basename = $out_folder."/".$out_name;
$out = $basename.".bam";

// check FASTQ quality encoding
$files = array_merge($in_for, $in_rev);
foreach($files as $file)
{
	list($stdout, $stderr) = $parser->exec(get_path("ngs-bits")."FastqFormat", "-in $file", true);
	if (!contains($stdout[2], "Sanger"))
	{
		trigger_error("Input file '$file' is not in Sanger/Illumina 1.8 format!", E_USER_ERROR);
	}
}

// perform adapter trimming, read QC and merging
if ($sys["adapter1_p5"]=="" && $sys["adapter2_p7"]=="")
{
	trigger_error("No forward and/or reverse adapter sequence given!\nForward: ".$sys["adapter1_p5"]."\nReverse: ".$sys["adapter2_p7"], E_USER_ERROR);
}
$trimmed1 = $parser->tempFile("_trimmed1.fastq.gz");
$trimmed2 = $parser->tempFile("_trimmed2.fastq.gz");
$stafile1 = $basename."_stats_fastq.qcML";
$parser->exec(get_path("ngs-bits")."SeqPurge", "-in1 ".implode(" ", $in_for)." -in2 ".implode(" ", $in_rev)." -out1 $trimmed1 -out2 $trimmed2 -a1 ".$sys["adapter1_p5"]." -a2 ".$sys["adapter2_p7"]." -qc $stafile1 -threads ".($threads>2 ? 2 : 1), true);

// MIPs: move molecular barcode to separate file
$index_file = $basename."_index.fastq.gz";
if($sys['type']=="Panel MIPs")
{
	$trimmed_mips2 = $parser->tempFile("_mips_indices.fastq.gz"); ///@todo rename to [sample]_MB_001.fastq.gz and store to sample folder
	$parser->exec(get_path("ngs-bits")."FastqExtractBarcode", "-in $trimmed2 -out_main $trimmed_mips2 -cut 8 -out_index $index_file",true);
	$trimmed2 = $trimmed_mips2;
}

//HaloPlex HS special handling
if($sys['type']=="Panel Haloplex HS")
{
	//cut extra C at beginning of read2 and end of read1
	$trimmed_hs1 = $parser->tempFile("_trimmed_hs1.fastq.gz");
	$trimmed_hs2 = $parser->tempFile("_trimmed_hs2.fastq.gz");
	$parser->exec(get_path("ngs-bits")."FastqTrim", "-in $trimmed1 -out $trimmed_hs1 -end 1",true);
	$parser->exec(get_path("ngs-bits")."FastqTrim", "-in $trimmed2 -out $trimmed_hs2 -start 1",true);
	$trimmed1 = $trimmed_hs1;
	$trimmed2 = $trimmed_hs2;
	
	//merge index files (in case sample was distributed over several lanes)
	$index_file = $basename."_index.fastq.gz";
	if (file_exists($index_file))
	{
		unlink($index_file);
	}
	$index_files = glob("$out_folder/*_index_*.fastq.gz");
	if (count($index_files)>0)
	{
		$parser->exec("cat",implode(" ",$index_files)." > $index_file",true);
	}
	
	if(!file_exists($index_file))
	{
		trigger_error("Index file for Haloplex HS enrichment $index_file was not found. The data is processed like a normal Haloplex panel without molecular barcode!", E_USER_WARNING);
	}
}

// mapping
$mapping_options = "";
if ($sys['shotgun']) $mapping_options .= " -dedup";
$parser->execTool("NGS/mapping_bwa.php", "-in1 $trimmed1 -in2 $trimmed2 -out $out $mapping_options -build ".$sys['build']." -threads ".$threads);

//perform indel realignment
if (!$no_abra && ($sys['target_file']!="" || $sys['type']=="WGS"))
{
	//execute
	$args = array();
	$args[] = "-in $out";
	$tmp_bam = $parser->tempFile("_indel_realign.bam");
	$args[] = "-out $tmp_bam";
	$args[] = "-build ".$sys['build'];
	$args[] = "-threads ".$threads;
	if ($sys['target_file']!="")
	{
		$args[] = "-roi ".$sys['target_file'];
	}
	else if ($sys['type']=="WGS") //for WGS use exome target region
	{
		$args[] = "-roi ".get_path("data_folder")."/enrichment/ssHAEv6_2017_01_05.bed";
	}
	
	$parser->execTool("NGS/indel_realign_abra.php", implode(" ", $args));
	
	//copy/index output
	copy2($tmp_bam, $out);
	$parser->exec(get_path("samtools")." index", " ".$out, true);
}

//remove too short reads from amplicon data
if (contains($sys['type'],"Haloplex")) //matches both HaloPlex and Haloplex HS
{
	$bam_clean = $parser->tempFile("_clipped1.bam");
	$parser->exec(get_path("ngs-bits")."BamCleanHaloplex -min_match 30", " -in $out -out $bam_clean", true);
	copy2($bam_clean, $out);
	$parser->exec(get_path("samtools")." index", " $out", true);
}

//clip overlapping reads
if($clip_overlap)
{
	$bam_clip1 = $parser->tempFile("_clipped1.bam");
	$parser->exec(get_path("ngs-bits")."BamClipOverlap", " -in $out -out $bam_clip1", true);
	$tmp1 = $parser->tempFile();
	$parser->exec(get_path("samtools"),"sort -T $tmp1 -o $out $bam_clip1", true);
	$parser->exec(get_path("samtools")." index", " $out", true);
}

//MIPs: remove duplicates by molecular barcode and cut extension/ligation arms
if($sys['type']=="Panel MIPs")
{
	$mip_file = isset($sys['mip_file']) ? $sys['mip_file'] : "/mnt/share/data/mipfiles/".$sys["name_short"].".txt";
	$bam_dedup1 = $parser->tempFile("_dedup1.bam");
	$parser->exec(get_path("ngs-bits")."BamDeduplicateByBarcode", " -bam $out -index $index_file -mip_file $mip_file -out $bam_dedup1 -del_amb -stats ".$basename."_bar_stats.tsv -dist 1", true);
	$tmp1 = $parser->tempFile();
	$parser->exec(get_path("samtools"),"sort -T $bam_dedup1 -o $out $bam_dedup1", true);
	$parser->exec(get_path("samtools")." index", " $out", true);
}

//HaloPlex HS: remove duplicates by molecular barcode
if($sys['type']=="Panel Haloplex HS" && file_exists($index_file))
{
	$bam_dedup1 = $parser->tempFile("_dedup1.bam");
	$min_group = isset($sys['min_group']) ? $sys['min_group'] : 1;
	$dist = isset($sys['dist']) ? $sys['dist'] : 1;
	$amplicon_file = substr($sys['target_file'], 0, -4)."_amplicons.bed";
	
	$parser->exec(get_path("ngs-bits")."BamDeduplicateByBarcode", " -bam $out -index $index_file -out $bam_dedup1 -min_group $min_group -stats ".$basename."_bar_stats.tsv -dist $dist -hs_file $amplicon_file", true);
	$tmp1 = $parser->tempFile();
	$parser->exec(get_path("samtools"),"sort -T $bam_dedup1 -o $out $bam_dedup1", true);
	$parser->exec(get_path("samtools")." index", " $out", true);
}

//run mapping QC
$stafile2 = $basename."_stats_map.qcML";
$params = array();
if ($sys['target_file']=="" || $sys['type']=="WGS")
{
	$params[] = "-wgs";
}
else
{
	$params[] = "-roi ".$sys['target_file'];
}
if ($sys['build']=="hg19" || $sys['build']=="GRCh37")
{
	$params[] = "-3exons";
}

$parser->exec(get_path("ngs-bits")."MappingQC", "-in $out -out $stafile2 ".implode(" ", $params), true);

?>
