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

$stafile1 = $basename."_stats_fastq.qcML";
$trimmed1 = $parser->tempFile("_trimmed1.fastq.gz");
$trimmed2 = $parser->tempFile("_trimmed2.fastq.gz");
$args = array();
if($sys['umi_type']!="n/a") $args[] = "-qcut 0 -ncut 0";
$parser->exec(get_path("ngs-bits")."SeqPurge", "-in1 ".implode(" ", $in_for)." -in2 ".implode(" ", $in_rev)." -out1 $trimmed1 -out2 $trimmed2 -a1 ".$sys["adapter1_p5"]." -a2 ".$sys["adapter2_p7"]." -qc $stafile1 -threads ".($threads>2 ? 2 : 1)." ".implode(" ",$args), true);

//extract barcodes and add to fastq header
$barcode_correction = false;
$index_file = $parser->tempFile("_index.fastq.gz");
if($sys['umi_type']!="n/a")
{
	//1. barcode within insert - Adapter trimming needed, MIPs (single UMI, 8 bp at the beginning of read1), ThruPlex (dual UMI, 6 bp both at beginning and end of pair) => extract to index file
	if($sys['umi_type']=="MIPs")
	{
		// move trimmed reads to temp
		$trimmed_mips2 = $parser->tempFile("_MB_001.fastq.gz");
		$parser->exec(get_path("ngs-bits")."FastqExtractBarcode", "-in $trimmed2 -out_main $trimmed_mips2 -cut 8 -out_index $index_file",true);
		$trimmed2 = $trimmed_mips2;		
	}
	//2. barcode within separate index read -  SureSelect HS (single UMI, index read 10 bp), HaloPlex HS (single UMI, index read 10 bp) => combine to single index file
	else if($sys['umi_type']=="HaloPlex HS" || $sys['umi_type']=="SureSelect HS" )
	{
		//HaloPlex HS mismatch bases at the end and beginning
		if($sys['umi_type']=="HaloPlex HS")
		{
			//cut extra C at beginning of read2 and end of read1
			$trimmed_hs1 = $parser->tempFile("_trimmed_hs1.fastq.gz");
			$trimmed_hs2 = $parser->tempFile("_trimmed_hs2.fastq.gz");
			$parser->exec(get_path("ngs-bits")."FastqTrim", "-in $trimmed1 -out $trimmed_hs1 -end 1",true);
			$parser->exec(get_path("ngs-bits")."FastqTrim", "-in $trimmed2 -out $trimmed_hs2 -start 1",true);
			$trimmed1 = $trimmed_hs1;
			$trimmed2 = $trimmed_hs2;
		}
		
		//merge index files (in case sample was distributed over several lanes)
		$index_files = glob("$out_folder/*_index_*.fastq.gz");
		if (count($index_files)>0) $parser->exec("cat",implode(" ",$index_files)." > $index_file",true);
	}
	else trigger_error("Unknown UMI-type ${sys['umi_type']}. Will be processed without molecular barcodes.",E_USER_WARNING);
	
	//3. annotate barcodes from index file to fastq header; fix read naming within fastq that is incorrect using barcode_to_header
	if(file_exists($index_file))
	{
		$barcode_correction = true;
		
		//move index to fastq header
		$trimmed1_bc = $parser->tempFile("_trimmed_bc1.fastq.gz");
		$parser->exec("python  ".repository_basedir()."/src/NGS/barcode_to_header.py", "-i $trimmed1 -bc1 $index_file -o $trimmed1_bc",true);
		$trimmed1 = $trimmed1_bc;
		$trimmed2_bc = $parser->tempFile("_trimmed_bc2.fastq.gz");
		$parser->exec("python  ".repository_basedir()."/src/NGS/barcode_to_header.py", "-i $trimmed2 -bc1 $index_file -o $trimmed2_bc",true);
		$trimmed2 = $trimmed2_bc;
	}
	else trigger_error("Index file '$index_file' was not found. The data will be processed without molecular barcodes!", E_USER_WARNING);
}

//clean up MIPs - remove single reads and MIP arms, skip all reads without perfect match mip arm
if($sys['type']=="Panel MIPs")
{
	$count1 = 0;
	$count2 = 0;
	
	//remove ligation and extension arms from read start (overlap clipping removes the arms at the end of a read later), filter for perfect match ligation/extension arms
	$mip_file = isset($sys['mip_file']) ? $sys['mip_file'] : "/mnt/share/data/mipfiles/".$sys["name_short"].".txt";
	if(!is_file($mip_file))	trigger_error("MIP file '".$sys['mip_file']."' is not available.",E_USER_ERROR);
	$mips = Matrix::fromTSV($mip_file,"\t",">");
	$tmp1 = $parser->tempFile("_R1_woarms.fastq.gz");
	$tmp2 = $parser->tempFile("_R2_woarms.fastq.gz");
	$idx_ext_seq = $mips->getColumnIndex("ext_probe_sequence");
	$idx_lig_seq = $mips->getColumnIndex("lig_probe_sequence");
	$idx_strand = $mips->getColumnIndex("probe_strand");
	$handle1 = gzopen($trimmed1, "r");
	if($handle1===FALSE) trigger_error("Could not open file $trimmed1 for reading.",E_USER_ERROR);
	$handle2 = gzopen($trimmed2, "r");
	if($handle2===FALSE) trigger_error("Could not open file $trimmed2 for reading.",E_USER_ERROR);
	$handle4 = gzopen($tmp1, "w");
	if($handle4===FALSE) trigger_error("Could not open file $tmp1 for writing.",E_USER_ERROR);
	$handle5 = gzopen($tmp2, "w");
	if($handle5===FALSE) trigger_error("Could not open file $tmp2 for writing.",E_USER_ERROR);
	while(!feof($handle1) && !feof($handle2))
	{
		$line1 = array(fgets($handle1),fgets($handle1),fgets($handle1),fgets($handle1));
		$line2 = array(fgets($handle2),fgets($handle2),fgets($handle2),fgets($handle2));
		
		list($id1,) = explode(" ", $line1[0]);
		list($id2,) = explode(" ", $line2[0]);
		if($id1!=$id2)	trigger_error("This should not happen.",E_USER_ERROR);
		
		$valid = false;
		for($i=0;$i<$mips->rows();++$i)
		{
			$mip = $mips->getRow($i);

			$pos1 = strpos($line1[1], rev_comp($mip[$idx_lig_seq]));
			$pos2 = strpos($line2[1], $mip[$idx_ext_seq]);
			
			if($pos1===0 && $pos2===0)
			{
				$line1[1] = substr($line1[1],strlen($mip[$idx_lig_seq]));
				$line1[3] = substr($line1[3],strlen($mip[$idx_lig_seq]));
				$line2[1] = substr($line2[1],strlen($mip[$idx_ext_seq]));
				$line2[3] = substr($line2[3],strlen($mip[$idx_ext_seq]));
				$valid = true;
			}
		}
		
		if($valid)
		{
			fwrite($handle4, implode("",$line1));
			fwrite($handle5, implode("",$line2));
			++$count1;
		}
		++$count2;
	}
	gzclose($handle1);
	gzclose($handle2);
	gzclose($handle4);
	gzclose($handle5);
	
	$trimmed1 = $tmp1;
	$trimmed2 = $tmp2;
	
	trigger_error(($count2-$count1)."/$count2 read pairs were removed by filtering for MIP arms.",E_USER_NOTICE);
}

// mapping
$bam_current = $parser->tempFile("_mapping_bwa.bam");
$args = array();
$args[] = "-in1 $trimmed1";
$args[] = "-in2 $trimmed2";
$args[] = "-out $bam_current";
$args[] = "-sample $out_name";
$args[] = "-build ".$sys['build'];
$args[] = "-threads $threads";
if ($sys['shotgun'] && !$barcode_correction && $sys['umi_type']!="ThruPLEX") $args[] = "-dedup";
$parser->execTool("NGS/mapping_bwa.php", implode(" ", $args));

//UMIs - ThruPlex; will be done by barcode_correction in the future
if($sys['umi_type']=="ThruPLEX")
{
	//keep bam before deduplication
	$parser->copyFile($bam_current, $basename."_before_dedup.bam");
	$parser->copyFile($bam_current.".bai", $basename."_before_dedup.bam.bai");
	
	$tmp_bam = $parser->tempFile("_connor_dedup.bam");
	$parser->exec(get_path("connor"), "$bam_current $tmp_bam --log_file ".$parser->tempFile("_connor_dedup.log"), true);
	$parser->indexBam($tmp_bam, $threads);
	$bam_current = $tmp_bam;
}

//UMIs - remove duplicates by molecular barcode
if($barcode_correction)
{
	//keep bam before deduplication
	$parser->copyFile($bam_current, $basename."_before_dedup.bam");
	$parser->copyFile($bam_current.".bai", $basename."_before_dedup.bam.bai");
	
	//additional mismatch filter for MIP experiments
	if($sys['umi_type']=="MIPs")
	{
		$tmp_bam1 = $parser->tempFile("_dedup1.bam");
		$parser->sortBam($bam_current,$tmp_bam1,$threads,TRUE);
		$tmp_bam2 = $parser->tempFile("_dedup2.bam");
		$parser->exec("python  ".repository_basedir()."/src/NGS/filter_bam.py", "--infile $tmp_bam1 --outfile $tmp_bam2",true);
		$tmp_bam3 = $parser->tempFile("_dedup3.bam");
		$parser->sortBam($tmp_bam2,$tmp_bam3,$threads);
		$parser->indexBam($tmp_bam3, $threads);
		$bam_current = $tmp_bam3;
	}

	//barcode correction
	$tmp_bam4 = $parser->tempFile("_dedup4.bam");
	$parser->exec("python  ".repository_basedir()."/src/NGS/barcode_correction.py", "--infile $bam_current --outfile $tmp_bam4 --step 3",true);
	$parser->indexBam($tmp_bam4, $threads);
	$parser->moveFile($tmp_bam4.".log", $out."_stats_dedup.tsv");
	
	$bam_current = $tmp_bam4;
}

//perform indel realignment
if (!$no_abra && ($sys['target_file']!="" || $sys['type']=="WGS"))
{
	$tmp_bam = $parser->tempFile("_indel_realign.bam");
	
	$args = array();
	$args[] = "-in $bam_current";
	$args[] = "-out $tmp_bam";
	$args[] = "-build ".$sys['build'];
	$args[] = "-threads ".$threads;
	if ($sys['type']=="WGS") //for WGS use exome target region
	{
		$args[] = "-roi ".get_path("data_folder")."/enrichment/ssHAEv6_2017_01_05.bed";
	}
	else if ($sys['target_file']!="")
	{
		$args[] = "-roi ".$sys['target_file'];
	}
	$parser->execTool("NGS/indel_realign_abra.php", implode(" ", $args));
	$parser->indexBam($tmp_bam, $threads);
	
	$bam_current = $tmp_bam;
}

//remove reads from HaloPlex data that are short
if ($sys['type']=="Panel Haloplex")
{
	$tmp_bam = $parser->tempFile("_clean.bam");
	$parser->exec(get_path("ngs-bits")."BamCleanHaloplex -min_match 30", "-in $bam_current -out $tmp_bam", true);
	$parser->indexBam($tmp_bam, $threads);
	
	$bam_current = $tmp_bam;
}

//clip overlapping reads
if($clip_overlap)
{
	$tmp_bam = $parser->tempFile("_clip_overlap_unsorted.bam");
	$parser->exec(get_path("ngs-bits")."BamClipOverlap", "-in $bam_current -out $tmp_bam", true);	
	$tmp_bam2 = $parser->tempFile("_clip_overlap_sorted.bam");
	$parser->sortBam($tmp_bam, $tmp_bam2, $threads);
	$parser->indexBam($tmp_bam2, $threads);

	$bam_current = $tmp_bam2;
}

//move BAM to final output location
$parser->moveFile($bam_current, $out);
$parser->moveFile($bam_current.".bai", $out.".bai");

//add baf file
$params = array();
if(!empty($sys['target_file'])) $params[] = "-target ".$sys['target_file'];
$parser->execTool("NGS/mapping_baf.php", "-in ${out} -out ${basename}_bafs.igv ".implode(" ", $params));

//run mapping QC
$stafile2 = $basename."_stats_map.qcML";
$params = array("-in $out", "-out $stafile2");
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
else
{
	$params[] = "-no_cont";
}
$parser->exec(get_path("ngs-bits")."MappingQC", implode(" ", $params), true);

?>
