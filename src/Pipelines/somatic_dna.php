<?php

/**
 * @page somatic_dna
 */

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("somatic_dna", "Differential analysis of tumor and normal DNA samples.");

//mandatory
$parser->addInfile("t_bam", "Tumor sample BAM file.", false);
$parser->addString("out_folder", "Output folder.", false);

//optional
$parser->addInfile("n_bam", "Normal sample BAM file.", true);
$parser->addInfile("t_rna_bam", "Tumor RNA sample BAM file.", true);

$parser->addString("prefix", "Output file prefix.", true, "somatic");

$steps_all = array("vc", "cn", "an", "ci", "msi","db");
$parser->addString("steps", "Comma-separated list of steps to perform:\n" .
	"vc=variant calling, an=annotation, ci=CGI annotation,\n" .
	"cn=copy-number analysis, msi=microsatellite analysis, db=database import",
	true, "vc,cn,an,msi,db");

$parser->addString("cancer_type", "Tumor type, see CancerGenomeInterpreter.org for nomenclature (resolved from GENLAB if not set).", true);

$parser->addInfile("system",  "Processing system file used for tumor DNA sample (resolved from NGSD via tumor BAM by default).", true);
$parser->addInfile("n_system",  "Processing system file used for normal DNA sample (resolved from NGSD via normal BAM by default).", true);

$parser->addFlag("skip_correlation", "Skip sample correlation check.");

$parser->addFlag("include_germline", "Include germline variant annotation with CGI.");

//default cut-offs
$parser->addFloat("min_af", "Allele frequency detection limit.", true, 0.05);
$parser->addFloat("min_correlation", "Minimum correlation for tumor/control pair.", true, 0.8);
$parser->addFloat("min_depth_t", "Tumor sample coverage cut-off for low coverage statistic.", true, 100);
$parser->addFloat("min_depth_n", "Control sample coverage cut-off for low coverage statistic.", true, 100);
$parser->addInt("min_cov_files", "Minimum number of required coverage files for CNV calling.", true, 20);

$parser->addInt("threads", "The maximum number of threads to use.", true, 4);
extract($parser->parse($argv));


###################################### AUXILARY FUNCTIONS ######################################

//Creates BAF file from "$gsvar" and "bam" file and writes content to "out_file"
function create_baf_file($gsvar,$bam,$out_file, $ref_genome)
{
	if(!file_exists($gsvar) || !file_exists($bam)) return;
	
	//Abort if out_file exists to prevent interference with other jobs
	if(file_exists($out_file)) return;
	
	$tmp_out = temp_file(".tsv");
	exec2(get_path("ngs-bits")."/VariantAnnotateFrequency -in {$gsvar} -bam {$bam} -depth -out {$tmp_out} -ref {$ref_genome}", true);
	
	$in_handle  = fopen($tmp_out,"r");
	$out_handle = fopen($out_file,"w");
	
	if($out_handle === false)
	{
		trigger_error("Could not open BAF file for writing: {$out_file}.",E_USER_WARNING);
	}
	
	while(!feof($in_handle))
	{
		$line = trim(fgets($in_handle));
		if(starts_with($line,"#")) continue;
		if(empty($line)) continue;
		
		$parts = explode("\t",$line);
		list($chr,$start,$end,$ref,$alt) = $parts;
		if(strlen($ref) != 1 || strlen($alt) != 1 ) continue; //Skip INDELs

		$freq = $parts[count($parts)-2]; //freq annotation ist 2nd last col
		$depth= $parts[count($parts)-1]; //depth annotation is last col
		fputs($out_handle,"{$chr}\t{$start}\t{$end}\t{$chr}_{$start}\t{$freq}\t{$depth}\n");
	}
	fclose($in_handle);
	fclose($out_handle);
}

//Checks $baf_folder for missing B-AF files and creates them if neccessary
function complement_baf_folder($t_n_id_file,$baf_folder,&$db_conn, $ref_genome)
{
	$ids = file($t_n_id_file);
	foreach($ids as $line)
	{
		if(starts_with($line,"#")) continue;
		list($tid,$nid) = explode(",",trim($line));
		if(!file_exists("{$baf_folder}/{$nid}.tsv"))
		{
			$ninfo = get_processed_sample_info($db_conn,$nid);
			$n_gsvar = $ninfo["ps_folder"] ."/{$nid}.GSvar";
			$n_bam = $ninfo["ps_bam"];
			create_baf_file($n_gsvar,$n_bam,"{$baf_folder}/{$nid}.tsv", $ref_genome);
		}
		if(!file_exists("{$baf_folder}/{$tid}.tsv"))
		{
			$ninfo = get_processed_sample_info($db_conn,$nid);
			$n_gsvar = $ninfo["ps_folder"] ."/{$nid}.GSvar";
			$tinfo = get_processed_sample_info($db_conn,$tid);
			$t_bam = $tinfo["ps_bam"];
			create_baf_file($n_gsvar,$t_bam,"{$baf_folder}/{$tid}.tsv", $ref_genome);
		}
	}
}

//parse and copy single CGI results files
function parse_cgi_single_result($from,$to)
{
	if(!file_exists($from)) return;
	addCommentCharInHeader($from);
	copy($from,$to);
}

function extract_and_move_cgi_results($zip_file,$dest_folder,$prefix = "CGI")
{
	$tmp_cgi_results_folder = temp_folder();
	exec2("unzip -n $zip_file -d $tmp_cgi_results_folder");
	
	parse_cgi_single_result("{$tmp_cgi_results_folder}/mutation_analysis.tsv","{$dest_folder}/{$prefix}_cgi_mutation_analysis.tsv");
	parse_cgi_single_result("{$tmp_cgi_results_folder}/cna_analysis.tsv","{$dest_folder}/{$prefix}_cgi_cnv_analysis.tsv");
	parse_cgi_single_result("{$tmp_cgi_results_folder}/fusion_analysis.txt","{$dest_folder}/{$prefix}_cgi_fusion_analysis.tsv");
	parse_cgi_single_result("{$tmp_cgi_results_folder}/drug_prescription.tsv","{$dest_folder}/{$prefix}_cgi_drug_prescription.tsv");
	parse_cgi_single_result("{$tmp_cgi_results_folder}/drug_prescription_bioactivities.tsv","{$dest_folder}/{$prefix}_cgi_drug_prescription_bioactivities.tsv");
	parse_cgi_single_result("{$tmp_cgi_results_folder}/malformed_cnas.txt","{$dest_folder}/{$prefix}_cgi_malformed_cnas.txt");
	parse_cgi_single_result("{$tmp_cgi_results_folder}/not_mapped_entries.txt","{$dest_folder}/{$prefix}_cgi_not_mapped_entries.txt");
}

//Check whether cancer acronym appears in acronym list
function is_valid_cgi_acronym($name)
{
	$cancer_acronyms_list = file(repository_basedir()."/data/misc/cancer_types.tsv");
	$valid_cgi_acronym = false;
	foreach($cancer_acronyms_list as $line)
	{
		list($acronym,,,) = explode("\t",$line);
		if($acronym == $name)
		{
			return true;
		}
	}
	return false;
}

//check steps
$steps = explode(",", $steps);
foreach($steps as $step)
{
	if (!in_array($step, $steps_all))
	{
		trigger_error("Unknown processing step '$step'!", E_USER_ERROR);
	}
}

###################################### SCRIPT START ######################################

//output prefix
$full_prefix = "{$out_folder}/{$prefix}";

//IDs, system and target region
$t_id = basename($t_bam, ".bam");
$sys = load_system($system, $t_id);
$roi = $sys["target_file"];
$ref_genome = genome_fasta($sys['build']);

//normal sample data (if not single sample analysis)
$single_sample = !isset($n_bam);
if (!$single_sample)
{
	$n_id = basename($n_bam, ".bam");
	$n_sys = load_system($n_system, $n_id);
	$ref_folder_n = get_path("data_folder")."/coverage/".$n_sys['name_short'];
}

//sample similiarty check
$bams = array_filter([$t_bam, $n_bam, $t_rna_bam]);
if (count($bams) > 1)
{
    if ($skip_correlation)
    {
        trigger_error("Genotype correlation check has been disabled!", E_USER_WARNING);
    }
    else
    {
    	$args_similarity = [
    		"-in", implode(" ", $bams),
			"-mode", "bam",
			"-max_snps", 4000
		];
		if (!empty($roi))
		{
			$args_similarity[] = "-roi";
			$args_similarity[] = $roi;
		}
        $output = $parser->exec(get_path("ngs-bits")."SampleSimilarity", implode(" ", $args_similarity), true);

		//extract columen 3 from output
		$table = array_map(
			function($str) { return explode("\t", $str); },
			array_slice($output[0], 1)
		);
		$correlation = array_column($table, 3);
        if (min($correlation) < $min_correlation)
        {
            trigger_error("Genotype correlation lower than {$min_correlation}!\n" . implode("\n", $output[0]), E_USER_ERROR);
        }
    }
}

//low coverage statistics

$low_cov = "{$full_prefix}_stat_lowcov.bed";					// low coverage BED file
if ($sys['type'] !== "WGS" && !empty($roi))
{
	if(!file_exists($low_cov)) $parser->exec(get_path("ngs-bits")."BedLowCoverage", "-in $roi -bam $t_bam -out $low_cov -cutoff $min_depth_t", true);
	//combined tumor and normal low coverage files
	//normal coverage is calculated only for tumor target region
	if(!$single_sample)
	{
		$low_cov_n = $parser->tempFile("_nlowcov.bed");
		if(!file_exists($low_cov_n)) $parser->exec(get_path("ngs-bits")."BedLowCoverage", "-in $roi -bam $n_bam -out $low_cov_n -cutoff $min_depth_n", true);
		$parser->execPipeline([
			[get_path("ngs-bits")."BedAdd", "-in $low_cov $low_cov_n"],
			[get_path("ngs-bits")."BedMerge", "-out $low_cov"]
		], "merge low coverage BED files");
	}
	// annotate with gene names
	if (db_is_enabled("NGSD"))
	{
		$parser->exec(get_path("ngs-bits") . "BedAnnotateGenes", "-in $low_cov -extend 25 -out $low_cov", true);
	}
}

//variant calling
$manta_indels = $full_prefix . "_var_smallIndels.vcf.gz";		// small indels from manta
$manta_sv     = $full_prefix . "_var_structural.vcf.gz";		// structural variants (vcf)
$manta_sv_tsv = $full_prefix . "_var_structural.tsv";			// structural variants (tsv)
$variants     = $full_prefix . "_var.vcf.gz";					// variants
$ballele      = $full_prefix . "_bafs.igv";						// B-allele frequencies
if (in_array("vc", $steps))
{
	// structural variant calling
	if (!$sys['shotgun'])
	{
		trigger_error("Structural variant calling deactivated for amplicon samples.", E_USER_NOTICE);
	}
	else if ($sys['umi_type'] === "ThruPLEX")
	{
		trigger_error("Structural variant calling deactivated for ThruPLEX samples.", E_USER_NOTICE);
	}
	else
	{
		$args_manta = [
			"-t_bam {$t_bam}",
			"-out {$manta_sv}",
			"-build ".$sys['build'],
			"-smallIndels {$manta_indels}",
			"-threads {$threads}",
			"-fix_bam"
		];
		if (!$single_sample)
		{
			$args_manta[] = "-bam $n_bam";
		}
		if ($sys['type'] === "WES")
		{
			$args_manta[] = "-exome";
		}
		if (!empty($roi))
		{
			$args_manta[] = "-target {$roi}";
		}
		$parser->execTool("NGS/vc_manta.php", implode(" ", $args_manta));
		$parser->execTool("Tools/converter_manta2tsv.php", "-in $manta_sv -out $manta_sv_tsv -tumor_id $t_id");
	}

	// variant calling
	if ($single_sample)	// variant calling using freebayes
	{
		// vc_freebayes uses ngs-bits that can not handle multi-sample vcfs
		$tmp1 = $parser->tempFile();
		$bams = $single_sample ? "$t_bam" : "$t_bam $n_bam";
		$parser->execTool("NGS/vc_freebayes.php", "-bam $bams -out $tmp1 -build ".$sys['build']." -no_ploidy -min_af $min_af -target ".$roi." -threads ".$threads);

		// find and rewrite header (tumor, normal sample columns)
		$tmp2 = $parser->tempFile();
		$s = Matrix::fromTSV($tmp1);
		$tmp_headers = $s->getHeaders();
		$tumor_col = NULL;
		$normal_col = NULL;
		for ($i=0; $i < count($tmp_headers); ++$i)
		{
			if (contains($tmp_headers[$i], $t_id))
			{
				$tumor_col = $i;
				$tmp_headers[$tumor_col] = $t_id;
			}
			else if (!$single_sample && contains($tmp_headers[$i],$n_id))
			{
				$normal_col = $i;
				$tmp_headers[$normal_col] = $n_id;
			}
		}
		$s->setHeaders($tmp_headers);

		// zip and index output file
		$s->toTSV($tmp2);
		$parser->exec("bgzip", "-c $tmp2 > $variants", true);
		$parser->exec("tabix", "-f -p vcf $variants", true);
	}
	else
	{
		$args_strelka = [
			"-t_bam {$t_bam}",
			"-n_bam {$n_bam}",
			"-out {$variants}",
			"-build ".$sys['build'],
			"-threads {$threads}"
		];
		if (!empty($roi))
		{
			$args_strelka[] = "-target {$roi}";
		}
		if (is_file($manta_indels))
		{
			$args_strelka[] = "-smallIndels {$manta_indels}";
		}
		$parser->execTool("NGS/vc_strelka2.php", implode(" ", $args_strelka));
	}

	//add somatic BAF file
	$baf_args = [
		"-in {$t_bam}",
		"-out {$ballele}",
		"-build ".$sys['build']
	];
	if (!$single_sample)
	{
		$baf_args[] = "-n_in $n_bam";
	}
	if (!empty($roi))
	{
		$baf_args[] = "-target {$roi}";
	}
	if(!$single_sample)
	{
		$variants_germline_vcf = dirname($n_bam)."/{$n_id}_var.vcf.gz";
		if (file_exists($variants_germline_vcf))
		{
			$baf_args[] = "-sites {$variants_germline_vcf}";
		}
	}
	$parser->execTool("NGS/mapping_baf.php", implode(" ", $baf_args));
}

//CNV calling
$som_cnv = $full_prefix . "_cnvs.tsv"; //CNVHunter output file
$som_clincnv = $full_prefix . "_clincnv.tsv"; //ClinCNV output file
if(in_array("cn",$steps))
{
	// copy number variant calling
	$tmp_folder = $parser->tempFolder();
	
	//Generate file with off-target region
	$off_target_bed = get_path("data_folder")."/coverage/off_target_beds/".$sys['name_short'].".bed";
	if(!file_exists($off_target_bed))
	{
		//generate bed file that contains edges of whole reference genome
		$ref_content = array();
		$ref_borders = file("{$ref_genome}.fai"); 
		foreach($ref_borders as $line)
		{
			$line = trim($line);
			list($chr,$chr_length) = explode("\t",$line);
			if(chr_check($chr,22,false) === false) continue;
			$chr_length--; //0-based coordinates
			$ref_content[] = "{$chr}\t0\t{$chr_length}\n";
		}
		$ref_bed = temp_file(".bed");
		file_put_contents($ref_bed,$ref_content);
		$ngs_bits = get_path("ngs-bits");
		$tmp_bed = temp_file(".bed");
		exec2("{$ngs_bits}BedExtend -in ".$sys['target_file']." -n 1000 -fai {$ref_genome}.fai | {$ngs_bits}BedMerge -out {$tmp_bed}");
		exec2("{$ngs_bits}BedSubtract -in ".$ref_bed." -in2 {$tmp_bed} | {$ngs_bits}BedChunk -n 100000 | {$ngs_bits}BedShrink -n 25000 | {$ngs_bits}BedExtend -n 25000 -fai {$ref_genome}.fai | {$ngs_bits}BedAnnotateGC -ref {$ref_genome} | {$ngs_bits}BedAnnotateGenes -out {$off_target_bed}"); //@TODO: Review parameters after a few samples have been analyzed
	}
	
	/***************************************************
	 * GENERATE AND COPY COVERAGE FILES TO DATA FOLDER *
	 ***************************************************/
	// coverage for tumor sample
	$t_cov = "{$tmp_folder}/{$t_id}.cov";
	$parser->exec(get_path("ngs-bits")."BedCoverage", "-min_mapq 0 -decimals 4 -bam $t_bam -in $roi -out $t_cov",true); //@TODO Leave min_mapq at 0 as long as we have CNVHunter, when replaced, set to 2 or 3 and recalculate all .cov files
	
	// coverage for tumor sample (off-target)
	$t_cov_off_target = "{$tmp_folder}/{$t_id}_off_target.cov";
	$parser->exec(get_path("ngs-bits")."BedCoverage","-min_mapq 10 -decimals 4 -in $off_target_bed -bam $t_bam -out $t_cov_off_target",true);
	
	//folders with tumor reference coverage files (of same processing system)
	$ref_folder_t = get_path("data_folder")."/coverage/".$sys['name_short']."-tumor";
	$ref_folder_t_off_target = $ref_folder_t . "_off_target";
	
	//copy tumor sample coverage file to reference folder (has to be done before ClinCNV call to avoid analyzing the same sample twice)
	if (db_is_enabled("NGSD") && is_valid_ref_tumor_sample_for_cnv_analysis($t_id))
	{
		//create reference folder if it does not exist
		create_directory($ref_folder_t);
		
		//copy file
		$ref_file_t = "{$ref_folder_t}/{$t_id}.cov";
		
		if(!file_exists($ref_file_t)) //Do not overwrite existing reference files in cov folder
		{
			$parser->copyFile($t_cov, $ref_file_t); 
			$t_cov = $ref_file_t;
		}
		
		//create reference folder for tumor off target coverage files
		create_directory($ref_folder_t_off_target);
		
		$ref_file_t_off_target = "{$ref_folder_t_off_target}/{$t_id}.cov";
		if(!file_exists($ref_file_t_off_target))
		{
			$parser->copyFile($t_cov_off_target,$ref_file_t_off_target);
			$t_cov_off_target = $ref_file_t_off_target;
		}
	}

	if($single_sample) //use CNVHunter in case of single samples
	{
		trigger_error("Tumor without normal sample is not supported by ClinCNV. Skipping 'cn' step.", E_USER_WARNING);
	}
	else //ClinCNV for differential sample
	{
		// coverage for normal sample
		$n_cov = "{$tmp_folder}/{$n_id}.cov";
		$parser->exec(get_path("ngs-bits")."BedCoverage", "-min_mapq 0 -decimals 4 -bam $n_bam -in ".$n_sys['target_file']." -out $n_cov", true);
		
		// coverage for normal sample (off-target)
		$n_cov_off_target = "{$tmp_folder}/{$n_id}_off_target.cov";
		$parser->exec(get_path("ngs-bits")."BedCoverage","-min_mapq 10 -decimals 4 -in $off_target_bed -bam $n_bam -out $n_cov_off_target",true);
		
		// copy normal sample coverage file to reference folder (only if valid).
		if (db_is_enabled("NGSD") && is_valid_ref_sample_for_cnv_analysis($n_id))
		{
			//create reference folder if it does not exist
			$ref_folder_n = get_path("data_folder")."/coverage/".$n_sys['name_short'];
			create_directory($ref_folder_n);
			
			//copy file
			$ref_file_n = $ref_folder_n."/".$n_id.".cov";
			if (!file_exists($ref_file_n)) //do not overwrite existing coverage files in ref folder
			{
				$parser->copyFile($n_cov, $ref_file_n);
				$n_cov = $ref_file_n;
			}
			
			$ref_folder_n_off_target = $ref_folder_n . "_off_target";
			create_directory($ref_folder_n_off_target);
			
			$ref_file_n_off_target = "{$ref_folder_n_off_target}/{$n_id}.cov";
			
			if(!file_exists($ref_file_n_off_target))
			{
				$parser->copyFile($n_cov_off_target,$ref_file_n_off_target);
				$t_cov_off_target = $ref_file_t_off_target;
			}
		}
		
		//append tumor-normal IDs to list with tumor normal IDs (stored in same folder as tumor coverage files)
		$t_n_list_file = $ref_folder_t . "/" . "list_tid-nid.csv";
		
		if (!file_exists($t_n_list_file))
		{
			$header = "##THIS FILE CONTAINS TUMOR AND NORMAL IDS OF PROCESSING SYSTEM ".$n_sys['name_short']."\n";
			$header .= "#tumor_id,normal_id\n";
			file_put_contents($t_n_list_file,$header);
		}
		
		//use temporary list file if n or t cov files are not valid
		if(!db_is_enabled("NGSD") || !is_valid_ref_sample_for_cnv_analysis($n_id) || !is_valid_ref_tumor_sample_for_cnv_analysis($t_id))
		{
			$tmp_file_name = temp_file(".csv");
			$parser->copyFile($t_n_list_file,$tmp_file_name);
			$t_n_list_file = $tmp_file_name;
		}
		
		//Append tumor-normal pair to csv list in tumor coverage folder if not included yet
		$already_in_list = false;
		foreach(file($t_n_list_file) as $line)
		{
			if(starts_with($line,'#')) continue;
			if(empty(trim($line))) continue;
			list($t,$n) = explode(",",trim($line));
			if($t == $t_id && $n == $n_id)
			{
				$already_in_list = true;
				break;
			}
		}
		if(!$already_in_list)
		{
			file_put_contents($t_n_list_file,"{$t_id},{$n_id}\n", FILE_APPEND | LOCK_EX);
		}

		$baf_folder = get_path("data_folder")."/coverage/". $sys['name_short']."_bafs";
		create_directory($baf_folder);
		if(db_is_enabled("NGSD"))
		{
			$db_conn = DB::getInstance("NGSD");
			//Create BAF file for each sample with the same processing system if not existing
			complement_baf_folder($t_n_list_file,$baf_folder,$db_conn, $ref_genome);
		}
		else
		{
			$normal_gsvar = dirname($n_bam) ."/{$n_id}.GSvar";
			create_baf_file($normal_gsvar,$n_bam,"{$baf_folder}/{$n_id}.tsv", $ref_genome);
			create_baf_file($normal_gsvar,$t_bam,"{$baf_folder}/{$t_id}.tsv", $ref_genome);
		}

		//Skip CNV Calling if there are less than 5 tumor-normal coverage pairs
		if(count(file($t_n_list_file)) > 7)
		{
			/*******************
			 * EXECUTE CLINCNV *
			 *******************/
			$args_clincnv = [
			"-t_id", $t_id,
			"-n_id", $n_id,
			"-t_cov", $t_cov,
			"-n_cov", $n_cov,
			"-out", $som_clincnv,
			"-cov_pairs",$t_n_list_file,
			"-system", $system,
			"-bed", $sys['target_file'],
			"-t_cov_off", $t_cov_off_target,
			"-n_cov_off", $n_cov_off_target,
			"-bed_off", $off_target_bed,
			"-baf_folder", $baf_folder,
			"-threads {$threads}"
			];
			$parser->execTool("NGS/vc_clincnv_somatic.php",implode(" ",$args_clincnv));
		}
		else
		{
			print("Not enough reference tumor-normal coverage files for processing system {$system} found. Skipping CNV calling.\n");
		}
	}
}

//annotation
$variants_annotated = $full_prefix . "_var_annotated.vcf.gz";	// annotated variants
$variants_gsvar     = $full_prefix . ".GSvar";					// GSvar variants
$somaticqc          = $full_prefix . "_stats_som.qcML";			// SomaticQC qcML
if (in_array("an", $steps))
{
	// annotate vcf (in temp folder)
	$tmp_folder1 = $parser->tempFolder();
	$tmp_vcf = "{$tmp_folder1}/{$prefix}_var_annotated.vcf.gz";
	$parser->execTool("Pipelines/annotate.php", "-out_name $prefix -out_folder $tmp_folder1 -system $system -vcf $variants -somatic -updown");

	// run somatic QC
	if (!$single_sample)
	{
		$links = array_filter([
			dirname($t_bam)."/{$t_id}_stats_fastq.qcML",
			dirname($t_bam)."/{$t_id}_stats_map.qcML",
			dirname($n_bam)."/{$n_id}_stats_fastq.qcML",
			dirname($n_bam)."/{$n_id}_stats_map.qcML"
		], "file_exists");

		$args_somaticqc = [
			"-tumor_bam", $t_bam,
			"-normal_bam", $n_bam,
			"-somatic_vcf", $tmp_vcf,
			"-target_bed", $roi,
			"-ref_fasta", $ref_genome,
			"-out", $somaticqc
		];
		if (!empty($links))
		{
			$args_somaticqc[] = "-links";
			$args_somaticqc[] = implode(" ", $links);
		}

		$parser->exec(get_path("ngs-bits")."SomaticQC", implode(" ", $args_somaticqc), true);
	}

	//add sample info to VCF header
	$s = Matrix::fromTSV($tmp_vcf);
	$comments = $s->getComments();
	$comments[] = gsvar_sample_header($t_id, array("IsTumor" => "yes"), "#", "");
	if (!$single_sample)
	{
		$comments[] = gsvar_sample_header($n_id, array("IsTumor" => "no"), "#", "");
	}
	$s->setComments(sort_vcf_comments($comments));
	$s->toTSV($tmp_vcf);

	// zip and index vcf file
	$parser->exec("bgzip", "-c $tmp_vcf > $variants_annotated", true);
	$parser->exec("tabix", "-f -p vcf $variants_annotated", true);

	// convert vcf to GSvar
	$args = array("-in $tmp_vcf", "-out $variants_gsvar", "-t_col $t_id");
	if (!$single_sample) $args[] = "-n_col $n_id";
	$parser->execTool("NGS/vcf2gsvar_somatic.php", implode(" ", $args));

	// annotate NGSD data
	if (db_is_enabled("NGSD"))
	{
		$parser->exec(get_path("ngs-bits")."VariantAnnotateNGSD", "-in $variants_gsvar -out $variants_gsvar -psname $t_id -mode somatic", true);
		$parser->exec(get_path("ngs-bits")."VariantAnnotateNGSD", "-in $variants_gsvar -out $variants_gsvar -psname $t_id -mode germline", true);
	}

	//annotate vcf/GSvar with frequency/depth from tumor RNA sample
	if (isset($t_rna_bam))
	{
		$tmp_vcf_rna = $parser->tempFile("_var_rna.vcf");
		$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency", "-in $tmp_vcf -bam $t_rna_bam -out $tmp_vcf_rna -name rna_tum -depth -ref $ref_genome", true);
		$parser->exec("bgzip", "-c $tmp_vcf_rna > $variants_annotated", true);
		$parser->exec("tabix", "-f -p vcf $variants_annotated", true);

		$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency", "-in $variants_gsvar -bam $t_rna_bam -out $variants_gsvar -name rna_tum -depth -ref $ref_genome", true);
	}
}

//QCI/CGI annotation
//@TODO: implementation for translocation files
$variants_qci = $full_prefix . "_var_qci.vcf.gz";				//CGI annotated vcf file
if (in_array("ci", $steps))
{
	// add QCI output
	$parser->execTool("Tools/converter_vcf2qci.php", "-in $variants_annotated -t_id $t_id -n_id $n_id -out $variants_qci -pass");
	
	/*********************************
	 * GET CGI CANCER_TYPE FROM NGSD *
	 *********************************/
	if(db_is_enabled("NGSD"))
	{
		$db_conn = DB::getInstance("NGSD");
		$t_name = explode("_",$t_id)[0];
		$sample_id = $db_conn->getValue("SELECT id FROM sample WHERE name = '{$t_name}'");
		$disease_info = $db_conn->executeQuery("SELECT * FROM sample_disease_info WHERE sample_id = '{$sample_id}'","");
		$cancer_type_from_ngsd = "";
		foreach($disease_info as $res)
		{
			if(array_key_exists("type",$res) && $res["type"] == "CGI cancer type")
			{
				$cancer_type_from_ngsd = $res["disease_info"];
				break;
			}
		}
		
		//Abort if no cancer type is set
		if($cancer_type_from_ngsd == "" && !isset($cancer_type))
		{
			trigger_error("There is no CGI cancer type set in NGSD and no parameter \"-cancer_type\" given. Aborting CGI analysis.",E_USER_ERROR);
		}
		//Abort if cancer type from parameter and NGSD differ
		if($cancer_type_from_ngsd != "" && isset($cancer_type) && $cancer_type_from_ngsd != $cancer_type)
		{
			trigger_error("Cancer type \"{$cancer_type_from_ngsd}\" in NGSD and -cancer_type \"{$cancer_type}\" differ.", E_USER_ERROR);
		}
		
		//Insert CGI cancer_type into NGSD if not set
		if($cancer_type_from_ngsd == "" && isset($cancer_type))
		{
			$time = get_timestamp(false);
			$user = exec2("whoami")[0][0];
			$user_id = $db_conn->getValue("SELECT id FROM user WHERE user_id = '{$user}' AND active='1'",-1);
			
			//If user is not in NGSD, we use the user "noone"
			if($user_id == -1) $user_id = 1;
			
			if(is_valid_cgi_acronym($cancer_type))
			{
				$db_conn->executeStmt("INSERT INTO sample_disease_info (`sample_id`,`disease_info`,`type`,`user_id`,`date`) VALUES ('{$sample_id}','{$cancer_type}','CGI cancer type','{$user_id}','{$time}')");
			}
			else
			{
				trigger_error("Could not insert CGI cancer type \"{$cancer_type}\" to NGSD table \"sample_disease_info\" because \"{$cancer_type}\" was not recognized as a valid CGI acronym.",E_USER_WARNING);
			}
		}
		
		if($cancer_type_from_ngsd != "") $cancer_type = $cancer_type_from_ngsd;
	}
	if(!is_valid_cgi_acronym($cancer_type))
	{
		trigger_error("CGI cancer_type {$cancer_type} is invalid. Aborting CGI analysis.",E_USER_ERROR);
	}
	
	/*************************
	 * EXECUTE CGI_SEND_DATA *
	 *************************/
	$tmp_cgi_zipped_outfile = temp_file(".zip");
	$parameters = [
		"-out",$tmp_cgi_zipped_outfile,
		"-cancertype", $cancer_type
	];
	if (file_exists($variants_annotated)) $parameters[] = "-mutations $variants_annotated";
	if (file_exists($som_clincnv)) $parameters[] = "-cnas $som_clincnv";
	if (file_exists($manta_sv_tsv)) $parameters[] = "-translocations $manta_sv_tsv";
	//if we know genes in target region we set parameter for this file
	$genes_in_target_region =  dirname($sys["target_file"]) . "/" .basename($sys["target_file"],".bed")."_genes.txt";
	if(file_exists($genes_in_target_region)) $parameters[] = "-t_region $genes_in_target_region";
	//execution will not stop pipeline if it fails but display several errors
	$result_send_data = $parser->execTool("NGS/cgi_send_data.php", implode(" ",$parameters),false);
	$error_code = $result_send_data[2];
	if($error_code != 0)
	{
		trigger_error("step \"ci\" did not exit cleanly. Please check sample CGI files manually. cgi_send_data returned code ".$error_code,E_USER_WARNING);
	}
	
	/******************************
	 * UNZIP AND COPY CGI RESULTS *
	 ******************************/
	extract_and_move_cgi_results($tmp_cgi_zipped_outfile,$out_folder,$prefix);
	
	/*************************************************
	 * ANNOTATE GSVAR AND CLINCNV FILE WITH CGI DATA *
	 *************************************************/
	//only try annotate SNVS to GSVar file if $som_vann (variant file) was uploaded to CGI
	$cgi_snv_result_file = $full_prefix . "_cgi_mutation_analysis.tsv";
	if(file_exists($cgi_snv_result_file) && file_exists($variants_gsvar))
	{
		$parser->execTool("NGS/cgi_snvs_to_gsvar.php","-gsvar_in $variants_gsvar -cgi_snv_in $cgi_snv_result_file -out $variants_gsvar");
	}
	//annotate CGI cnv genes to cnv input file (which was originally created by CNVHunter)
	$cgi_cnv_result_file = $full_prefix . "_cgi_cnv_analysis.tsv";
	if(file_exists($som_clincnv) && file_exists($cgi_cnv_result_file))
	{
		$parser->execTool("NGS/cgi_annotate_cnvs.php","-cnv_in $som_clincnv -cnv_in_cgi $cgi_cnv_result_file -out $som_clincnv");
	}
	
	/*************************
	 * GERMLINE CGI ANALYSIS *
	 *************************/
	$variants_germline = dirname($n_bam)."/{$n_id}_var_annotated.vcf.gz";
	if(file_exists($variants_germline) && $include_germline)
	{
		$tmp_file = temp_file(".vcf");
		exec2("gzip -d -k -c  $variants_germline > $tmp_file");
		$tmp_file_content = file($tmp_file, FILE_IGNORE_NEW_LINES);
		$tmp_new_file = array("#CHROM\tPOS\tID\tREF\tALT\n");
		
		//Remove variants for CGI upload with more than 3% AF in 1000G/GnomAD
		foreach($tmp_file_content as $line)
		{
			if(starts_with($line,"#")) continue;
			list($chr,$pos,$id,$ref,$alt,$qual,$filter,$info) = explode("\t",$line);
		
			list($qual,$filters) = explode(";",$info);
			$filters = explode(",",$filters)[0]; //take 1st transcript (is equal in all transcripts?)
			$af_1000g = explode("|",$filters)[14];
			$af_gnomad = explode("|",$filters)[15];
			if(!empty($af_1000g) && floatval($af_1000g) > 0.03) continue;
			if(!empty($af_gnomad) && floatval($af_gnomad) > 0.03) continue;
			
			$tmp_new_file[] = "{$chr}\t{$pos}\t.\t{$ref}\t{$alt}\n";
		}
		$filtered_germline_variants = temp_file(".vcf");
		file_put_contents($filtered_germline_variants,$tmp_new_file);
		$tmp_cgi_zipped_outfile = temp_file();
		
		$params = [
			"-out",$tmp_cgi_zipped_outfile,
			"-mutations",$filtered_germline_variants,
			"-cancertype","CANCER"
		];
		
		//Include germline CNVs if available
		$germline_cnv_file = dirname($n_bam)."/{$n_id}_clincnv.tsv";
		if(!file_exists($germline_cnv_file)) $germline_cnv_file = dirname($n_bam)."/{$n_id}_cnvs.tsv";
		if(file_exists($germline_cnv_file))
		{
			$params[] = "-cnas $germline_cnv_file";
			if(file_exists($genes_in_target_region)) $params[] = "-t_region $genes_in_target_region";
		}
		
		$result_send_data = $parser->execTool("NGS/cgi_send_data.php",implode(" ",$params),false);
		
		//Error handling
		$error_code = $result_send_data[2];
		if($error_code != 0)
		{
			trigger_error("step \"ci\" for germline variants did not exit cleanly. Please check sample CGI files manually. cgi_send_data returned code ".$error_code,E_USER_WARNING);
		}

		//Move to germline Sample data
		$out_folder_germline = dirname($n_bam);
		extract_and_move_cgi_results($tmp_cgi_zipped_outfile,$out_folder_germline,$n_id);

		//Annotate
		$germline_gsvar_file = dirname($n_bam)."/{$n_id}.GSvar";
		if(file_exists("{$out_folder_germline}/{$n_id}_cgi_mutation_analysis.tsv") && file_exists($germline_gsvar_file))
		{
			$parameters = " -gsvar_in $germline_gsvar_file -cgi_snv_in {$out_folder_germline}/{$n_id}_cgi_mutation_analysis.tsv -out $germline_gsvar_file";
			$parser->execTool("NGS/cgi_snvs_to_gsvar.php",$parameters,false);
		}
		if(file_exists($germline_cnv_file) && file_exists("{$out_folder_germline}/{$n_id}_cgi_cnv_analysis.tsv"))
		{
			$parameters = "-cnv_in $germline_cnv_file -out $germline_cnv_file -cnv_in_cgi {$out_folder_germline}/{$n_id}_cgi_cnv_analysis.tsv";
			$parser->execTool("NGS/cgi_annotate_cnvs.php",$parameters,false);
		}
		
	}
	elseif($include_germline)
	{
		print("Could not create CGI report for germline variants because germline variant file $variants_germline does not exist.\n");
	}
}

//MSI calling
$msi_o_file = $full_prefix . "_msi.tsv";						//MSI
if (in_array("msi", $steps) && !$single_sample)
{
	//check whether file with loci exists in output folder
	//if not: intersect with loci file of reference
	$reference_loci_file = get_path("data_folder") . "/dbs/MANTIS/".$n_sys['build']."_msi_loci.bed";
	if(!file_exists($reference_loci_file))
	{
		print("Could not find loci reference file $reference_loci_file. Trying to generate it.\n");
		$parser->exec(get_path("mantis")."/tools/RepeatFinder","-i $ref_genome -o $reference_loci_file",false);
	}

	//file that contains MSI in target region -> is intersection of loci reference with target region
	$target_bed_file = $n_sys['target_file'];
	$target_loci_file = $full_prefix . "_msi_loci_target.bed";
	
	//target loci file is intersection of reference loci with target region
	if(!file_exists($target_loci_file))
	{
		$parameters = "-in ".$reference_loci_file." -in2 ".$target_bed_file ." -mode in -out ".$target_loci_file;
		$parser->exec(get_path("ngs-bits")."BedIntersect",$parameters,false);
	}

	$parameters = "-n_bam $n_bam -t_bam $t_bam -threads $threads -bed_file $target_loci_file -out " .$msi_o_file. " -build ".$n_sys['build'];
	
	if($sys['type'] == "WES") $parameters .= " -is_exome";
	
	$parser->execTool("NGS/detect_msi.php",$parameters);
}
elseif ($single_sample && in_array("msi",$steps))
{
	trigger_error("Calling microsatellite instabilities is only possible for tumor normal pairs",E_USER_NOTICE);
}

//DB import
if (in_array("db", $steps) && db_is_enabled("NGSD"))
{
	$db_conn = DB::getInstance("NGSD");
	
	$t_info = get_processed_sample_info($db_conn, $t_id, false);
	$n_info = get_processed_sample_info($db_conn, $n_id, false);
	
	if (is_null($t_info) || (!$single_sample && is_null($n_info)))
	{
		trigger_error("No database import since no valid processing ID (T: {$t_id}".($single_sample ? "" : " /N: {$n_id}").")", E_USER_WARNING);
	}
	else
	{
		// check sex
		$parser->execTool("NGS/db_check_gender.php", "-in $t_bam -pid $t_id");

		// import qcML files
		$log_db = dirname($t_bam)."/{$t_id}_log4_db.log";
		$qcmls = implode(" ", array_filter([
			dirname($t_bam)."/{$t_id}_stats_fastq.qcML",
			dirname($t_bam)."/{$t_id}_stats_map.qcML",
			$somaticqc
		], "file_exists"));
		$parser->execTool("NGS/db_import_qc.php", "-id $t_id -files $qcmls -force -min_depth 0 --log $log_db");

		// check tumor/normal flag
		if (!$t_info['is_tumor'])
		{
			trigger_error("Tumor sample $t_id is not flagged as tumor in NGSD!", E_USER_WARNING);
		}

		// analogous steps for normal sample, plus additional
		if (!$single_sample)
		{
			// check sex
			$parser->execTool("NGS/db_check_gender.php", "-in $n_bam -pid $n_id");

			// import qcML files
			$log_db = dirname($n_bam)."/{$n_id}_log4_db.log";
			$qcmls = implode(" ", array_filter([
				dirname($n_bam)."/{$n_id}_stats_fastq.qcML",
				dirname($n_bam)."/{$n_id}_stats_map.qcML",
				dirname($n_bam)."/{$n_id}_stats_vc.qcML"
			], "file_exists"));
			$parser->execTool("NGS/db_import_qc.php", "-id $n_id -files $qcmls -force -min_depth 0 --log $log_db");

			// check tumor/normal flag
			if ($n_info['is_tumor'])
			{
				trigger_error("Normal sample $n_id is flagged as tumor in NGSD!", E_USER_WARNING);
			}

			// update normal sample entry for tumor
			if (updateNormalSample($db_conn, $t_id, $n_id))
			{
				trigger_error("Updated normal sample ($n_id) for tumor ($t_id) in NGSD.", E_USER_NOTICE);
			}
			
			// import sample relation
			$s_id_t = $t_info['s_id'];
			$s_id_n = $n_info['s_id'];
			$db_conn->executeStmt("INSERT IGNORE INTO `sample_relations`(`sample1_id`, `relation`, `sample2_id`) VALUES ({$s_id_t},'tumor-normal',{$s_id_n})");

			// import variants (not for WGS)
			if (file_exists($variants_gsvar) && $sys['type'] !== "WGS")
			{
				$parser->execTool("NGS/db_import_variants.php", "-id {$t_id}-{$n_id} -var {$variants_gsvar} -build ".$sys['build']." -force -mode somatic");
			}
		}
	}
}
?>
