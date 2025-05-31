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
$parser->addInfile("out_folder", "Output folder.", false);
$parser->addString("out_name", "Output file base name (e.g. 'GS120001_01').", false);
//optional
$parser->addInfileArray("in_index",  "Index reads FASTQ file(s).", true);
$parser->addInfile("system",  "Processing system INI file (automatically determined from NGSD if 'out_name' is a valid processed sample name).", true);
$parser->addOutfile("local_bam", "Filename the local BAM file is written to. It can be used to speed up variant calling etc. Not created if unset.", true);
$parser->addInt("threads", "The maximum number of threads used.", true, 2);
$parser->addFlag("clip_overlap", "Soft-clip overlapping read pairs.");
$parser->addFlag("no_abra", "Skip realignment with ABRA.");
$parser->addFlag("no_trim", "Skip adapter trimming with SeqPurge.");
$parser->addFlag("correction_n", "Use Ns for barcode correction.");
$parser->addFlag("filter_bam", "Filter alignments prior to barcode correction.");
$parser->addFlag("bam_output", "Output is BAM instead of CRAM.");
$parser->addFlag("somatic_custom_map", "Calculate mapping QC metrics for somatic custom subpanel");
$parser->addInt("min_mapq", "The minimum mapping quality for reads to be considered for barcode correction (cfDNA only).", true, 0);
extract($parser->parse($argv));


//extract processing system information from DB
$sys = load_system($system, $out_name);
$build = $sys['build'];

// determine output file base name
$basename = $out_folder."/".$out_name;

// check FASTQ quality encoding
$files = array_merge($in_for, $in_rev);
foreach($files as $file)
{
	list($stdout, $stderr) = $parser->execApptainer("ngs-bits", "FastqFormat", "-in $file", [$file]);
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
$stafile1_dir = dirname($stafile1);


$trimmed1 = $parser->tempFile("_trimmed1.fastq.gz");
$trimmed2 = $parser->tempFile("_trimmed2.fastq.gz");

// barcode handling
$barcode_correction = false;

// input files for ngs-bits container
$in_files = $in_for;
$in_files = array_merge($in_files, $in_rev);

if (in_array($sys['umi_type'], ["HaloPlex HS", "SureSelect HS", "IDT-UDI-UMI"]))
{
	if ($in_index === false || empty($in_index))
	{
		trigger_error("Processing system ".$sys['name_short']." has UMI type ".$sys['umi_type'].", but no index files are specified => UMI-based de-duplication skipped!", E_USER_WARNING);
		$parser->execApptainer("ngs-bits", "SeqPurge", "-in1 ".implode(" ", $in_for)." -in2 ".implode(" ", $in_rev)." -out1 $trimmed1 -out2 $trimmed2 -a1 ".$sys["adapter1_p5"]." -a2 ".$sys["adapter2_p7"]." -qc $stafile1 -qcut 0 -ncut 0 -threads ".bound($threads-2, 1, 12), $in_files, [$stafile1_dir]);
	}
	else
	{
		// add barcodes to header
		$merged1_bc = $parser->tempFile("_bc1.fastq.gz");
		$merged2_bc = $parser->tempFile("_bc2.fastq.gz");
		$in_files = array_merge($in_files, $in_index);
		$parser->execApptainer("ngs-bits", "FastqAddBarcode", "-in1 ".implode(" ", $in_for)." -in2 ".implode(" ", $in_rev)." -in_barcode ".implode(" ", $in_index)." -out1 $merged1_bc -out2 $merged2_bc", $in_files);
		// run SeqPurge
		$parser->execApptainer("ngs-bits", "SeqPurge", "-in1 $merged1_bc -in2 $merged2_bc -out1 $trimmed1 -out2 $trimmed2 -a1 ".$sys["adapter1_p5"]." -a2 ".$sys["adapter2_p7"]." -qc $stafile1 -threads ".bound($threads-2, 1, 12), [], [$stafile1_dir]);

		// trim 1 base from end of R1 and 1 base from start of R2 (only HaloPlex HS)
		if ($sys['umi_type'] === "HaloPlex HS")
		{
			//cut extra C at beginning of read2 and end of read1
			$trimmed_hs1 = $parser->tempFile("_merged_hs1.fastq.gz");
			$trimmed_hs2 = $parser->tempFile("_merged_hs2.fastq.gz");

			$parser->execApptainer("ngs-bits", "FastqTrim", "-in $trimmed1 -out $trimmed_hs1 -end 1");
			$parser->execApptainer("ngs-bits", "FastqTrim", "-in $trimmed2 -out $trimmed_hs2 -start 1");
			

			$parser->deleteTempFile($trimmed1);
			$parser->deleteTempFile($trimmed2);
			$trimmed1 = $trimmed_hs1;
			$trimmed2 = $trimmed_hs2;
		}
		
		$barcode_correction = true;
	}
}
else if (in_array($sys['umi_type'], [ "MIPs", "ThruPLEX", "Safe-SeqS", "QIAseq", "IDT-xGen-Prism", "Twist"]))
{

	//handle UMI protocols where the UMI is placed at the head of read 1 and/or read 2

	//remove sequencing adapter
	$parser->execApptainer("ngs-bits", "SeqPurge", "-in1 ".implode(" ", $in_for)." -in2 ".implode(" ", $in_rev)." -out1 $trimmed1 -out2 $trimmed2 -a1 ".$sys["adapter1_p5"]." -a2 ".$sys["adapter2_p7"]." -qc $stafile1 -qcut 0 -ncut 0 -threads ".bound($threads-2, 1, 12), $in_files, [$stafile1_dir]);

	//set protocol specific UMI lengths
	switch ($sys['umi_type'])
	{
		case "MIPs":		//8bp from R2
			$cut1 = 0;
			$cut2 = 8;
			break;
		case "ThruPLEX":	//8bp each from R1, R2
		case "IDT-xGen-Prism":
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
		case "Twist":		//5bp each from R1, R2
			$cut1 = 5;
			$cut2 = 5;
			break;
	}

	$trimmed1_bc = $parser->tempFile("_trimmed1_bc.fastq.gz");
	$trimmed2_bc = $parser->tempFile("_trimmed2_bc.fastq.gz");
	$parser->execApptainer("ngs-bits", "FastqExtractUMI", "-in1 $trimmed1 -in2 $trimmed2 -out1 $trimmed1_bc -out2 $trimmed2_bc -cut1 $cut1 -cut2 $cut2");

	$parser->deleteTempFile($trimmed1);
	$parser->deleteTempFile($trimmed2);
	$trimmed1 = $trimmed1_bc;
	$trimmed2 = $trimmed2_bc;

	//special handling of Twist UMIs: cut off first and last 2bp
	if ($sys['umi_type'] == "Twist")
	{
		$trimmed1_bc = $parser->tempFile("_trimmed1_bc.fastq.gz");
		$trimmed2_bc = $parser->tempFile("_trimmed2_bc.fastq.gz");
		$parser->execApptainer("ngs-bits", "FastqTrim", "-in $trimmed1 -out $trimmed1_bc -start 2");
		$parser->execApptainer("ngs-bits", "FastqTrim", "-in $trimmed2 -out $trimmed2_bc -start 2");

		$parser->deleteTempFile($trimmed1);
		$parser->deleteTempFile($trimmed2);
		$trimmed1 = $trimmed1_bc;
		$trimmed2 = $trimmed2_bc;
	}

	if ($sys['umi_type'] === "QIAseq")
	{
		$trimmed2_trim = $parser->tempFile("_trimmed2_trim.fastq.gz");
		$parser->execApptainer("ngs-bits", "FastqTrim", "-in $trimmed2 -out $trimmed2_trim -start 11");
		$parser->deleteTempFile($trimmed2);
		$trimmed2 = $trimmed2_trim;
	}

	$barcode_correction = true;
}
else //normal analysis without UMIs
{
	if ($sys['umi_type']!=="n/a") trigger_error("Unknown UMI-type ".$sys['umi_type'].". No barcode correction.",E_USER_WARNING);

	if (!$no_trim)
	{
		$parser->execApptainer("ngs-bits", "SeqPurge", "-in1 ".implode(" ", $in_for)." -in2 ".implode(" ", $in_rev)." -out1 $trimmed1 -out2 $trimmed2 -a1 ".$sys["adapter1_p5"]." -a2 ".$sys["adapter2_p7"]." -qc $stafile1 -threads ".bound($threads-2, 1, 12), $in_files, [$stafile1_dir]);
	}
}

//clean up MIPs - remove single reads and MIP arms, skip all reads without perfect match mip arm
if($sys['type']=="Panel MIPs")
{
	$count1 = 0;
	$count2 = 0;
	
	//remove ligation and extension arms from read start (overlap clipping removes the arms at the end of a read later), filter for perfect match ligation/extension arms
	$mip_file = isset($sys['mip_file']) ? $sys['mip_file'] : get_path("data_folder")."/mipfiles/".$sys["name_short"].".txt";
	if(!file_exists($mip_file))	trigger_error("MIP file '{$mip_file}' is not available.",E_USER_ERROR);
	$mips = Matrix::fromTSV($mip_file,"\t",">");
	$tmp1 = $parser->tempFile("_R1_woarms.fastq.gz");
	$tmp2 = $parser->tempFile("_R2_woarms.fastq.gz");
	$idx_ext_seq = $mips->getColumnIndex("ext_probe_sequence");
	$idx_lig_seq = $mips->getColumnIndex("lig_probe_sequence");
	$handle1 = gzopen2($trimmed1, "r");
	$handle2 = gzopen2($trimmed2, "r");
	$handle4 = gzopen2($tmp1, "w");
	$handle5 = gzopen2($tmp2, "w");
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
			if ($pos1===FALSE) continue;
			$pos2 = strpos($line2[1], $mip[$idx_ext_seq]);
			if ($pos2===FALSE) continue;
			
			//$parser->log(__LINE__.": mip_index=$i pos1=$pos1 pos2=$pos2");
			
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

//mapping with BAW-mem2
$bam_current = $parser->tempFile("_mapping_bwa.bam");
$args = array();
$args[] = "-in1 $trimmed1";
$args[] = "-in2 $trimmed2";
$args[] = "-out $bam_current";
$args[] = "-sample $out_name";
$args[] = "-build ".$build;
$args[] = "-threads $threads";
if ($sys['shotgun'] && !$barcode_correction && $sys['umi_type']!="ThruPLEX") $args[] = "-dedup";
$parser->execTool("Tools/mapping_bwa.php", implode(" ", $args));

// remove temporary FASTQs (they can be really large for WGS)
$parser->deleteTempFile($trimmed1);
$parser->deleteTempFile($trimmed2);

//perform indel realignment
if (!$no_abra && ($sys['target_file']!="" || $sys['type']=="WGS"))
{
	$tmp_bam = $parser->tempFile("_indel_realign.bam");

	$args = array();
	$args[] = "-in $bam_current";
	$args[] = "-out $tmp_bam";
	$args[] = "-build ".$build;
	$args[] = "-threads ".$threads;
	if ($sys['type']=="WGS") //for WGS use exome target region
	{
		$args[] = "-roi ".repository_basedir()."/data/gene_lists/gene_exons_pad20.bed";
	}
	else if ($sys['target_file']!="")
	{
		$args[] = "-roi ".$sys['target_file'];
	}
	$parser->execTool("Tools/indel_realign_abra.php", implode(" ", $args));
	$parser->indexBam($tmp_bam, $threads);

	// delete local bam file or DRAGEN bam file
	$parser->deleteTempFile($bam_current);
	$bam_current = $tmp_bam;
}

//UMIs - remove duplicates by molecular barcode
if($barcode_correction)
{
	//keep bam before deduplication
	$parser->copyFile($bam_current, $basename."_before_dedup.bam");
	$parser->copyFile($bam_current.".bai", $basename."_before_dedup.bam.bai");
	
	//additional mismatch filter for MIP experiments
	if($sys['umi_type']=="MIPs" || $filter_bam)
	{
		$tmp_bam_filtered = $parser->tempFile("_filtered.bam");
		$parser->execApptainer("ngs-bits", "BamFilter", "-in $bam_current -out $tmp_bam_filtered", [$bam_current]);

		$tmp_bam_filtered_sorted = $parser->tempFile("_filtered_sorted.bam");
		$parser->sortBam($tmp_bam_filtered, $tmp_bam_filtered_sorted, $threads, $build);
		$parser->indexBam($tmp_bam_filtered_sorted, $threads);
		
		$parser->deleteTempFile($bam_current);
		$bam_current = $tmp_bam_filtered_sorted;
	}

	// filter BAM file by minMQ
	if($min_mapq > 0)
	{
		$tmp_bam_filtered = $parser->tempFile("_filtered.bam");
		$parser->execApptainer("ngs-bits", "BamFilter", "-minMQ $min_mapq -in $bam_current -out $tmp_bam_filtered", [$bam_current]);

		$tmp_bam_filtered_sorted = $parser->tempFile("_filtered_sorted.bam");
		$parser->sortBam($tmp_bam_filtered, $tmp_bam_filtered_sorted, $threads, $build);
		$parser->indexBam($tmp_bam_filtered_sorted, $threads);
		
		$parser->deleteTempFile($bam_current);
		$bam_current = $tmp_bam_filtered_sorted;
	}

	//barcode correction
	$tmp_bam4 = $parser->tempFile("_dedup4.bam");
	$tmp_bam4_sorted = $parser->tempFile("_dedup4_sorted.bam");
	$args = array();
	if($correction_n) $args[] = "--n ";
	
	// allow UMI barcode errors for IDT/Twist
	if($sys['umi_type'] == "Twist") $args[] = "--barcode_error 2";
	elseif($sys['umi_type'] == "IDT-xGen-Prism") $args[] = "--barcode_error 3";
	
	// use the barcode correction of umiVar2
	$parser->execApptainer("umiVar", "barcode_correction.py", "--infile $bam_current --outfile $tmp_bam4 ".implode(" ", $args));
	$parser->sortBam($tmp_bam4, $tmp_bam4_sorted, $threads, $build);
	$parser->indexBam($tmp_bam4_sorted, $threads);
	
	$parser->deleteTempFile($bam_current);
	$bam_current = $tmp_bam4_sorted;
}

//remove reads from HaloPlex data that are short
if ($sys['type']=="Panel Haloplex")
{
	$tmp_bam = $parser->tempFile("_clean.bam");
	$parser->execApptainer("ngs-bits", "BamCleanHaloplex", "-min_match 30 -in $bam_current -out $tmp_bam", [$bam_current]);
	$parser->indexBam($tmp_bam, $threads);
	
	$parser->deleteTempFile($bam_current);
	$bam_current = $tmp_bam;
}

//clip overlapping reads
if($clip_overlap)
{
	$tmp_bam = $parser->tempFile("_clip_overlap_unsorted.bam");
	$parser->execApptainer("ngs-bits", "BamClipOverlap", "-in $bam_current -out $tmp_bam -overlap_mismatch_basen", [$bam_current]);
	$tmp_bam2 = $parser->tempFile("_clip_overlap_sorted.bam");
	$parser->sortBam($tmp_bam, $tmp_bam2, $threads, $build);
	$parser->indexBam($tmp_bam2, $threads);


	// delete local bam file or DRAGEN bam file
	$parser->deleteTempFile($bam_current);
	$bam_current = $tmp_bam2;
}

//run mapping QC
$stafile2 = $basename."_stats_map.qcML";
$in_files = [
	$bam_current,
	genome_fasta($build)
];
$out_files = [dirname($stafile2)];

$params = array("-in $bam_current", "-out $stafile2", "-ref ".genome_fasta($build), "-build ".ngsbits_build($build));
if ($sys['target_file']=="" || $sys['type']=="WGS" || $sys['type']=="WGS (shallow)")
{
	$params[] = "-wgs";
}
else
{
	$params[] = "-roi ".realpath($sys['target_file']);
	$in_files[] = $sys['target_file'];
}
if ($build!="GRCh38")
{
	$params[] = "-no_cont";
}
$somatic_custom_panel = get_path("data_folder") . "/enrichment/somatic_VirtualPanel_v5.bed";
if ($somatic_custom_map && file_exists($somatic_custom_panel))
{
	$params[] = "-somatic_custom_bed $somatic_custom_panel";
	$in_files[] = $somatic_custom_panel;
}
if (!file_exists($stafile1))
{
	$params[] = "-read_qc $stafile1";
	$out_files[] = $stafile1_dir;
}
$parser->execApptainer("ngs-bits", "MappingQC", implode(" ", $params), $in_files, $out_files);

//create CRAM/BAM in output folder
if ($bam_output)
{
	$out_bam = $basename.".bam";
	$parser->copyFile($bam_current, $out_bam);
	$parser->copyFile($bam_current.".bai", $out_bam.".bai");
}
else
{
	$out_cram = $basename.".cram";
	$parser->execTool("Tools/bam_to_cram.php", "-bam {$bam_current} -cram {$out_cram} -build {$build} -threads {$threads}");
	
	//in case we re-map an old analysis with BAM output, we need to delete the BAM file
	$bam = $basename.".bam";
	if (file_exists($bam)) unlink($bam);
	$bai = $basename.".bam.bai";
	if (file_exists($bai)) unlink($bai);
}

// rename tmp BAM to allow using it for variant calling etc
if (!empty($local_bam))
{
	$parser->moveFile($bam_current, $local_bam);
	$parser->moveFile($bam_current.".bai", $local_bam.".bai");
}


?>
