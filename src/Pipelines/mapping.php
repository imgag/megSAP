<?php

/**
	@page mapping
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
		trigger_error("Input file '$in_for' is not in Sanger/Illumina 1.8 format!", E_USER_ERROR);
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
if($sys['type']=="Panel MIPs")
{
	$trimmed_mips2 = $parser->tempFile("_mips_indices.fastq.gz"); ///@todo rename to [sample]_MB_001.fastq.gz and store to sample folder
	$mips_indices = $basename."_index.fastq.gz";
	$parser->exec(get_path("ngs-bits")."FastqExtractBarcode", "-in $trimmed2 -out_main $trimmed_mips2 -cut 8 -out_index $mips_indices",true);
	$trimmed2 = $trimmed_mips2;
}

/*
//HaloPlex special handling
if($sys['type']=="Panel Haloplex")
{
	//cut enzyme footprint at beginning of read2 and end of read1
	$trimmed_hs1 = $parser->tempFile("_trimmed_hs1.fastq.gz");
	$trimmed_hs2 = $parser->tempFile("_trimmed_hs2.fastq.gz");
	$parser->exec(get_path("ngs-bits")."FastqTrim", "-in $trimmed1 -out $trimmed_hs1 -start 1 -end 1", true);
	$parser->exec(get_path("ngs-bits")."FastqTrim", "-in $trimmed2 -out $trimmed_hs2 -start 1 -end 1", true);
	$trimmed1 = $trimmed_hs1;
	$trimmed2 = $trimmed_hs2;
}
*/
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
	
	//merge index files
	$index_hs = $basename."_index.fastq.gz";
	$index_files = glob("$out_folder/*_index_*.fastq.gz");
	$index_files = array_diff(glob("$out_folder/*_index_*.fastq.gz"),array($index_hs));
	$parser->exec("cat",implode(" ",$index_files)." > $index_hs",true);
}

// mapping
$mapping_options = "";
if (file_exists($sheet)) $mapping_options .= " -sheet $sheet";
if ($sys['shotgun']) $mapping_options .= " -dedup";
$parser->execTool("NGS/mapping_bwa.php", "-in1 $trimmed1 -in2 $trimmed2 -out $out $mapping_options -build ".$sys['build']." -threads ".$threads);

//perform indel realignment
if (!$no_abra)
{
	$abra_args = "-in $out -out $out -build ".$sys['build']." -mer 0.1 -threads ".$threads;
	if ($sys['target_file']!="" && $sys['type']!="WGS") //for panels on complete target region
	{
		$parser->execTool("NGS/indel_realign_abra.php", "$abra_args -roi ".$sys['target_file']);
	}
	else if ($sys['build']=="hg19" && $sys['type']=="WGS") //for WGS on exome target region
	{
		$parser->execTool("NGS/indel_realign_abra.php", "$abra_args -roi ".get_path("data_folder")."/enrichment/ssHAEv6_2016_09_01.bed");
	}
}

//remove too short reads from amplicon data
if (contains($sys['type'],"Haloplex")) //matches both HaloPlex and Haloplex HS
{
	$bam_clean = $parser->tempFile("_clipped1.bam");
	$parser->exec(get_path("ngs-bits")."BamCleanHaloplex -min_match 30", " -in $out -out $bam_clean", true);
	copy2($bam_clean, $out);
	$parser->exec(get_path("ngs-bits")."BamIndex", "-in $out", true);
}

//clip overlapping reads
if($clip_overlap)
{
	$bam_clip1 = $parser->tempFile("_clipped1.bam");
	$parser->exec(get_path("ngs-bits")."BamClipOverlap", " -in $out -out $bam_clip1", true);
	$tmp1 = $parser->tempFile();
	$parser->exec(get_path("samtools"),"sort -T $tmp1 -o $out $bam_clip1", true);
	$parser->exec(get_path("ngs-bits")."BamIndex", "-in $out", true);
}

//MIPs: remove duplicates by molecular barcode and cut extension/ligation arms
if($sys['type']=="Panel MIPs") //@todo Add pipeline test for MIPs/HaloPlex HS (only a few exons)
{
	$bam_dedup1 = $parser->tempFile("_dedup1.bam");
	$parser->exec(get_path("ngs-bits")."BamDeduplicateByBarcode", " -bam $out -index ".$basename."_index.fastq.gz -mip_file /mnt/share/data/mipfiles/".$sys["name_short"].".txt -out $bam_dedup1 -stats ".$basename."_bar_stats.tsv -dist 1", true);
	$tmp1 = $parser->tempFile();
	$parser->exec(get_path("samtools"),"sort -T $tmp1 -o $out $bam_dedup1", true);
	$parser->exec(get_path("ngs-bits")."BamIndex", "-in $out", true);
}

//HaloPlex HS: remove duplicates by molecular barcode
if($sys['type']=="Panel Haloplex HS")
{
	$bam_dedup1 = $parser->tempFile("_dedup1.bam");
	$min_group = 1;
	if(isset($sys['min_group']))	$min_group = $sys['min_group'];
	$dist = 1;
	if(isset($sys['dist']))	$dist = $sys['dist'];
	
	$extension_pos = strrpos($sys['target_file'], '.'); // find position of the last dot, so where the extension starts
	$amplicon_file=substr($sys['target_file'], 0, $extension_pos) . '_amplicons' . substr($sys['target_file'], $extension_pos);//add '_amplicons' just before the dot (i.e. at end of file base name)
	
	$parser->exec(get_path("ngs-bits")."BamDeduplicateByBarcode", " -bam $out -index $index_hs -out $bam_dedup1 -min_group $min_group -stats ".$basename."_bar_stats.tsv -dist $dist -hs_file $amplicon_file", true);
	$tmp1 = $parser->tempFile();
	$parser->exec(get_path("samtools"),"sort -T $tmp1 -o $out $bam_dedup1", true);
	$parser->exec(get_path("ngs-bits")."BamIndex", "-in $out", true);
}

// run mapping QC
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
if ($sys['build']=="hg19")
{
	$params[] = "-3exons";
}

$parser->exec(get_path("ngs-bits")."MappingQC", "-in $out -out $stafile2 ".implode(" ", $params), true);

?>
