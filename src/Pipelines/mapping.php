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
$parser->addInfileArray("in_index",  "Index reads FASTQ file(s).", true);
$parser->addInfile("system",  "Processing system INI file (automatically determined from NGSD if 'out_name' is a valid processed sample name).", true);
$parser->addInt("threads", "The maximum number of threads used.", true, 2);
$parser->addFlag("clip_overlap", "Soft-clip overlapping read pairs.", true);
$parser->addFlag("no_abra", "Skip realignment with ABRA.", true);
$parser->addFlag("correction_n", "Use Ns for barcode correction.", true);
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

// barcode handling
$barcode_correction = false;
if ($sys['umi_type'] === "HaloPlex HS" || $sys['umi_type'] === "SureSelect HS" )
{
	$index_files = glob("$out_folder/*_index_*.fastq.gz");

	if ($in_index === NULL || empty($in_index))
	{
		trigger_error("Processing system ".$sys['name_short']." has UMI type ".$sys['umi_type'].", but no index files are specified => UMI-based de-duplication skipped!", E_USER_WARNING);
		$parser->exec(get_path("ngs-bits")."SeqPurge", "-in1 ".implode(" ", $in_for)." -in2 ".implode(" ", $in_rev)." -out1 $trimmed1 -out2 $trimmed2 -a1 ".$sys["adapter1_p5"]." -a2 ".$sys["adapter2_p7"]." -qc $stafile1 -qcut 0 -ncut 0 -threads ".bound($threads, 1, 4), true);
	}
	else
	{
		// add barcodes to header
		$merged1_bc = $parser->tempFile("_bc1.fastq.gz");
		$merged2_bc = $parser->tempFile("_bc2.fastq.gz");
		$parser->exec(get_path("ngs-bits")."FastqAddBarcode", "-in1 ".implode(" ", $in_for)." -in2 ".implode(" ", $in_rev)." -in_barcode ".implode(" ", $in_index)." -out1 $merged1_bc -out2 $merged2_bc", true);
		// run SeqPurge
		$parser->exec(get_path("ngs-bits")."SeqPurge", "-in1 $merged1_bc -in2 $merged2_bc -out1 $trimmed1 -out2 $trimmed2 -a1 ".$sys["adapter1_p5"]." -a2 ".$sys["adapter2_p7"]." -qc $stafile1 -threads ".($threads>2 ? 2 : 1), true);

		// trim 1 base from end of R1 and 1 base from start of R2 (only HaloPlex HS)
		if ($sys['umi_type'] === "HaloPlex HS")
		{
			//cut extra C at beginning of read2 and end of read1
			$trimmed_hs1 = $parser->tempFile("_merged_hs1.fastq.gz");
			$trimmed_hs2 = $parser->tempFile("_merged_hs2.fastq.gz");

			$parser->exec(get_path("ngs-bits")."FastqTrim", "-in $trimmed1 -out $trimmed_hs1 -end 1",true);
			$parser->exec(get_path("ngs-bits")."FastqTrim", "-in $trimmed2 -out $trimmed_hs2 -start 1",true);
			

			$parser->deleteTempFile($trimmed1);
			$parser->deleteTempFile($trimmed2);
			$trimmed1 = $trimmed_hs1;
			$trimmed2 = $trimmed_hs2;
		}
		
		$barcode_correction = true;
	}
}
else if (in_array($sys['umi_type'], [ "MIPs", "ThruPLEX", "Safe-SeqS", "QIAseq" ]))
{
	//handle UMI protocols where the UMI is placed at the head of read 1 and/or read 2

	//remove sequencing adapter
	$parser->exec(get_path("ngs-bits")."SeqPurge", "-in1 ".implode(" ", $in_for)." -in2 ".implode(" ", $in_rev)." -out1 $trimmed1 -out2 $trimmed2 -a1 ".$sys["adapter1_p5"]." -a2 ".$sys["adapter2_p7"]." -qc $stafile1 -qcut 0 -ncut 0 -threads ".($threads>2 ? 2 : 1), true);

	//set protocol specific UMI lengths
	switch ($sys['umi_type'])
	{
		case "MIPS":		//8bp from R1
			$cut1 = 8;
			$cut2 = 0;
			break;
		case "ThruPLEX":	//8bp each from R1, R2
			$cut1 = 8;
			$cut2 = 8;
			break;
		case "Safe-SeqS":	//12bp from R1
			$cut1 = 12;
			$cut2 = 0;
			break;
		case "QIAseq":		//12bp from R2
			$cut1 = 0;
			$cut2 = 12;
			break;
	}

	$trimmed1_bc = $parser->tempFile("_trimmed1_bc.fastq.gz");
	$trimmed2_bc = $parser->tempFile("_trimmed2_bc.fastq.gz");
	$parser->exec(get_path("ngs-bits")."FastqExtractUMI", "-in1 $trimmed1 -in2 $trimmed2 -out1 $trimmed1_bc -out2 $trimmed2_bc -cut1 $cut1 -cut2 $cut2", true);
	
	$parser->deleteTempFile($trimmed1);
	$parser->deleteTempFile($trimmed2);
	$trimmed1 = $trimmed1_bc;
	$trimmed2 = $trimmed2_bc;

	$barcode_correction = true;
}
else
{
	if ($sys['umi_type']!=="n/a") trigger_error("Unknown UMI-type ".$sys['umi_type'].". No barcode correction.",E_USER_WARNING);
	$parser->exec(get_path("ngs-bits")."SeqPurge", "-in1 ".implode(" ", $in_for)." -in2 ".implode(" ", $in_rev)." -out1 $trimmed1 -out2 $trimmed2 -a1 ".$sys["adapter1_p5"]." -a2 ".$sys["adapter2_p7"]." -qc $stafile1 -threads ".($threads>2 ? 2 : 1), true);
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

// remove temporary FASTQs (they can be really large for WGS)
$parser->deleteTempFile($trimmed1);
$parser->deleteTempFile($trimmed2);

//UMIs - remove duplicates by molecular barcode
if($barcode_correction)
{
	//keep bam before deduplication
	$parser->copyFile($bam_current, $basename."_before_dedup.bam");
	$parser->copyFile($bam_current.".bai", $basename."_before_dedup.bam.bai");
	
	//additional mismatch filter for MIP experiments
	if($sys['umi_type']=="MIPs")
	{
		$tmp_bam_filtered = $parser->tempFile("_filtered.bam");
		$parser->exec(get_path("ngs-bits")."BamFilter", "-in $bam_current -out $tmp_bam_filtered", true);

		$tmp_bam_filtered_sorted = $parser->tempFile("_filtered_sorted.bam");
		$parser->sortBam($tmp_bam_filtered, $tmp_bam_filtered_sorted, $threads);
		$parser->indexBam($tmp_bam_filtered_sorted, $threads);
		
		$parser->deleteTempFile($bam_current);
		$bam_current = $tmp_bam_filtered_sorted;
	}

	//barcode correction
	$tmp_bam4 = $parser->tempFile("_dedup4.bam");
	$args = "";
	if($correction_n) $args .= "--n ";
	$parser->exec("python  ".repository_basedir()."/src/NGS/barcode_correction.py", "--infile $bam_current --outfile $tmp_bam4 --step 3 ".$args,true);
	$parser->indexBam($tmp_bam4, $threads);
	$parser->moveFile($tmp_bam4.".log", $out."_stats_dedup.tsv");
	
	$parser->deleteTempFile($bam_current);
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
		$args[] = "-roi ".get_path("data_folder")."/gene_lists/genes_exons.bed";
	}
	else if ($sys['target_file']!="")
	{
		$args[] = "-roi ".$sys['target_file'];
	}
	$parser->execTool("NGS/indel_realign_abra.php", implode(" ", $args));
	$parser->indexBam($tmp_bam, $threads);
	
	$parser->deleteTempFile($bam_current);
	$bam_current = $tmp_bam;
}

//remove reads from HaloPlex data that are short
if ($sys['type']=="Panel Haloplex")
{
	$tmp_bam = $parser->tempFile("_clean.bam");
	$parser->exec(get_path("ngs-bits")."BamCleanHaloplex -min_match 30", "-in $bam_current -out $tmp_bam", true);
	$parser->indexBam($tmp_bam, $threads);
	
	$parser->deleteTempFile($bam_current);
	$bam_current = $tmp_bam;
}

//clip overlapping reads
if($clip_overlap)
{
	$tmp_bam = $parser->tempFile("_clip_overlap_unsorted.bam");
	$parser->exec(get_path("ngs-bits")."BamClipOverlap", "-in $bam_current -out $tmp_bam -overlap_mismatch_basen", true);
	$tmp_bam2 = $parser->tempFile("_clip_overlap_sorted.bam");
	$parser->sortBam($tmp_bam, $tmp_bam2, $threads);
	$parser->indexBam($tmp_bam2, $threads);

	$parser->deleteTempFile($bam_current);
	$bam_current = $tmp_bam2;
}

//move BAM to final output location
$parser->moveFile($bam_current, $out);
$parser->moveFile($bam_current.".bai", $out.".bai");

//add baf file
$params = array();
if ($sys['type']!="WGS" && !empty($sys['target_file']) && in_array($sys['build'], [ "hg19", "GRCh37" ]))
{
	$params[] = "-target ".$sys['target_file'];
}
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
if ($sys['build']!="hg19" && $sys['build']!="GRCh37")
{
	$params[] = "-no_cont";
}
$parser->exec(get_path("ngs-bits")."MappingQC", implode(" ", $params), true);

?>
