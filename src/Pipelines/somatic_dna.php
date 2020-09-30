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

$steps_all = array("vc", "vi", "cn", "an", "ci", "msi", "an_rna", "db");
$parser->addString("steps", "Comma-separated list of steps to perform:\n" .
	"vc=variant calling, an=annotation, ci=CGI annotation,\n" .
	"cn=copy-number analysis, msi=microsatellite analysis,\n".
	"an_rna=annotate data from somatic RNA files,\n".
	"vi=virus detection, db=database import",
	true, "vc,cn,an,msi,vi,db");

$parser->addString("cancer_type", "Tumor type, see CancerGenomeInterpreter.org for nomenclature (resolved from GENLAB if not set).", true);
$parser->addInfile("system",  "Processing system file used for tumor DNA sample (resolved from NGSD via tumor BAM by default).", true);
$parser->addInfile("n_system",  "Processing system file used for normal DNA sample (resolved from NGSD via normal BAM by default).", true);

$parser->addFlag("skip_correlation", "Skip sample correlation check.");
$parser->addFlag("skip_low_cov", "Skip low coverage statistics.");
$parser->addFlag("include_germline", "Include germline variant annotation with CGI.");

//default cut-offs
$parser->addFloat("min_af", "Allele frequency detection limit (for tumor-only calling only).", true, 0.05);
$parser->addFloat("min_correlation", "Minimum correlation for tumor/normal pair.", true, 0.8);
$parser->addFloat("min_depth_t", "Tumor sample coverage cut-off for low coverage statistics.", true, 100);
$parser->addFloat("min_depth_n", "Normal sample coverage cut-off for low coverage statistics.", true, 100);
$parser->addInt("min_cov_files", "Minimum number of required coverage files for CNV calling.", true, 20);
$parser->addString("cnv_baseline_pos","baseline region for ClinCNV, format e.g. chr1:12-12532",true);
$parser->addString("rna_ref_tissue", "Reference data for RNA annotation", true);
$parser->addInt("threads", "The maximum number of threads to use.", true, 4);
extract($parser->parse($argv));


###################################### AUXILARY FUNCTIONS ######################################

//Creates BAF file from "$gsvar" and "bam" file and writes content to "out_file"
function create_baf_file($gsvar,$bam,$out_file, $ref_genome)
{
	global $parser;
	
	if(!file_exists($gsvar) || !file_exists($bam)) return;
	
	//Abort if out_file exists to prevent interference with other jobs
	if(file_exists($out_file)) return;
	
	$tmp_out = $parser->tempFile(".tsv");
	exec2(get_path("ngs-bits")."/VariantAnnotateFrequency -in {$gsvar} -bam {$bam} -depth -out {$tmp_out} -ref {$ref_genome}", true);
	
	$in_handle  = fopen2($tmp_out,"r");
	$out_handle = fopen2($out_file,"w");
	
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
	global $parser;
	
	$tmp_cgi_results_folder = $parser->tempFolder();
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

//set up local NGS data copy (to reduce network traffic and speed up analysis)
$parser->execTool("Tools/data_setup.php", "-build ".$sys['build']);

//normal sample data (if not single sample analysis)
$single_sample = !isset($n_bam);
if (!$single_sample)
{
	$n_id = basename($n_bam, ".bam");
	$n_sys = load_system($n_system, $n_id);
	$ref_folder_n = get_path("data_folder")."/coverage/".$n_sys['name_short'];
	
	//Check whether both samples have same processing system
	if($roi != $n_sys["target_file"])
	{
		trigger_error("Tumor sample $t_id and normal sample $n_id have different target regions. Aborting...",E_USER_ERROR);
	}
}

//Disable steps "vc" and "cnv" if somatic report config exists in NGSD
if (!$single_sample && db_is_enabled("NGSD"))
{
	$db = DB::getInstance("NGSD", false);
	list($config_id, $config_vars_exist, $config_cnvs_exist) = somatic_report_config($db, $t_id, $n_id);
	if (in_array("vc", $steps) && $config_vars_exist)
	{
		trigger_error("Somatic report configuration with SNVs exists in NGSD! Delete somatic report configuration for reanalysis of step 'vc'.", E_USER_ERROR);
	}
	if (in_array("cn", $steps) && $config_cnvs_exist)
	{
		trigger_error("Somatic report configuration with CNVs exists in NGSD! Delete somatic report configuration for reanalysis of step 'cn'.", E_USER_ERROR);
	}
}

//sample similiarty check
$bams = array_filter([$t_bam, $n_bam]);
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


//Check samples are flagged correctly in NGSD
if( db_is_enabled("NGSD") && count($bams) > 1 )
{
	$db = DB::getInstance("NGSD");
	
	$tinfo = get_processed_sample_info($db, $t_id, false);

	if(!is_null($tinfo) && $tinfo["is_tumor"] != 1)
	{
		trigger_error("Please check tumor processed sample {$t_id} in NGSD. The sample is not flagged as tumor tissue.", E_USER_ERROR);
	}
	
	$ninfo = get_processed_sample_info($db, $n_id, false);
	if(!is_null($ninfo) && $ninfo["is_tumor"] != 0)
	{
		trigger_error("Please check normal processed sample {$n_id} in NGSD. The sample is flagged as tumor tissue.", E_USER_ERROR);
	}
}

//low coverage statistics

$low_cov = "{$full_prefix}_stat_lowcov.bed";					// low coverage BED file
if ($sys['type'] !== "WGS" && !empty($roi) && !$skip_low_cov)
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
$manta_indels  = $full_prefix . "_manta_var_smallIndels.vcf.gz";		// small indels from manta
$manta_sv      = $full_prefix . "_manta_var_structural.vcf.gz";		// structural variants (vcf)
$manta_sv_bedpe= $full_prefix . "_manta_var_structural.bedpe"; 		// structural variants (bedpe)
$variants      = $full_prefix . "_var.vcf.gz";					// variants
$ballele       = $full_prefix . "_bafs.igv";					// B-allele frequencies

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
			"-threads {$threads}"
		];
		if (!$single_sample)
		{
			$args_manta[] = "-bam $n_bam";
		}
		if ($sys['type'] !== "WGS") //use exome flag for non targeted / exome samples (i.e. non WGS samples)
		{
			$args_manta[] = "-exome";
		}
		if (!empty($roi))
		{
			$args_manta[] = "-target {$roi}";
		}
		$parser->execTool("NGS/vc_manta.php", implode(" ", $args_manta));
		
		exec2(get_path("ngs-bits") . "VcfToBedpe -in $manta_sv -out $manta_sv_bedpe");
		if(!$single_sample)
		{
			$parser->execTool("Tools/bedpe2somatic.php", "-in $manta_sv_bedpe -out $manta_sv_bedpe -tid $t_id -nid $n_id");
		}
		
		if( db_is_enabled("NGSD") )
		{
			$parser->exec(get_path("ngs-bits") . "BedpeGeneAnnotation", "-in $manta_sv_bedpe -out $manta_sv_bedpe -add_simple_gene_names", true );
		}
	}

	// variant calling
	if ($single_sample)	// variant calling using varscan2
	{
		$parser->execTool("NGS/vc_varscan2.php", "-bam $t_bam -out $variants -build " .$sys['build']. " -target ". $roi. " -name $t_id -min_af $min_af");
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
		if ($sys['type'] === "WGS")
		{
			$args_strelka[] = "-wgs";
		}
		if (is_file($manta_indels))
		{
			$args_strelka[] = "-smallIndels {$manta_indels}";
		}
		$parser->execTool("NGS/vc_strelka2.php", implode(" ", $args_strelka));
	}

	//add somatic BAF file
	if (!$single_sample)
	{
		$variants_germline_vcf = dirname($n_bam)."/{$n_id}_var.vcf.gz";
		if (file_exists($variants_germline_vcf))
		{
			$baf_args = [
				"-bam_t {$t_bam}",
				"-bam_n {$n_bam}",
				"-vcf {$variants_germline_vcf}",
				"-out {$ballele}",
				"-build ".$sys['build']
			];
			
			$parser->execTool("NGS/baf_somatic.php", implode(" ", $baf_args));
		}
	}
}

//Viral sequences alignment
$tumor_prefix = dirname($t_bam) . "/" . basename($t_bam, ".bam");

$viral         = $tumor_prefix . "_viral.tsv";					// viral sequences results
$viral_bam     = $tumor_prefix . "_viral.bam";					// viral sequences alignment
$viral_bam_raw = $tumor_prefix . "_viral_before_dedup.bam";		// viral sequences alignment (no deduplication)
$viral_bed = get_path("data_folder") . "/enrichment/somatic_viral.bed"; //viral enrichment
$viral_igv     = $tumor_prefix . "_viral.xml";					// IGV session
if (in_array("vi", $steps))
{
	$genome_rel = relative_path(dirname($viral_igv), get_path("data_folder") . "/genomes/somatic_viral.fa");
	if(!file_exists($genome_rel) || !file_exists($viral_bed))
	{
		trigger_error("Could not find reference genome {$genome_rel}. Skipping step \"vi\".", E_USER_WARNING);
	}
	else
	{
		//detection of viral sequences
		$t_bam_dedup = dirname($t_bam) . "/" . basename($t_bam, ".bam") . "_before_dedup.bam";
		$t_bam_map_qc = dirname($t_bam) . "/" . basename($t_bam, ".bam") . "_stats_map.qcML";
		$dedup_used = file_exists($t_bam_dedup);
		$vc_viral_args = [
			"-in", $dedup_used ? $t_bam_dedup : $t_bam,
			"-viral_bam", $viral_bam,
			"-viral_bam_raw", $viral_bam_raw,
			"-viral_cov", $viral,
			"-viral_chrs", "chrNC_007605",
			"-build_viral", "somatic_viral",
			"-in_qcml", $t_bam_map_qc
		];
		if ($dedup_used)
		{
			$vc_viral_args[] = "-barcode_correction";
		}
		$parser->execTool("NGS/vc_viral_load.php", implode(" ", $vc_viral_args));

		$igv_tracks = implode(" ", array_filter([
			$viral_bam,
			$viral_bam_raw,
			$viral_bed],
			"file_exists"));
		$parser->execTool("NGS/igv_session.php", "-genome {$genome_rel} -out {$viral_igv} -in {$igv_tracks} -relative");
	}
}

//CNV calling
$som_clincnv = $full_prefix . "_clincnv.tsv"; //ClinCNV output file
if(in_array("cn",$steps))
{
	// copy number variant calling
	$tmp_folder = $parser->tempFolder();
	
	//Generate file with off-target region
	$off_target_bed = get_path("data_folder")."/coverage/off_target_beds/".$sys['name_short'].".bed";
	if(!file_exists($off_target_bed))	create_off_target_bed_file($off_target_bed,$sys['target_file'], $ref_genome);
	
	/***************************************************
	 * GENERATE AND COPY COVERAGE FILES TO DATA FOLDER *
	 ***************************************************/
	// coverage for tumor sample
	$t_cov = "{$tmp_folder}/{$t_id}.cov";
	$parser->exec(get_path("ngs-bits")."BedCoverage", "-min_mapq 0 -decimals 4 -bam $t_bam -in $roi -out $t_cov",true);
	
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

	/*if($single_sample) //use CNVHunter in case of single samples
	{
		trigger_error("Currently only tumor normal pairs are supported for ClinCNV calls. Using CNVHunter on tumor sample instead.", E_USER_WARNING);
		$parser->execTool("NGS/vc_cnvhunter.php","-cov {$t_cov} -out {$som_clincnv} -system {$system} -min_corr 0 -seg {$t_id} -n {$min_cov_files}");

		// annotate CNV file
		$repository_basedir = repository_basedir();
		$data_folder = get_path("data_folder");
		$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$som_clincnv} -in2 {$repository_basedir}/data/misc/cn_pathogenic.bed -no_duplicates -url_decode -out {$som_clincnv}", true);
		$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$som_clincnv} -in2 {$data_folder}/dbs/ClinGen/dosage_sensitive_disease_genes.bed -no_duplicates -url_decode -out {$som_clincnv}", true);
		$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$som_clincnv} -in2 {$data_folder}/dbs/ClinVar/clinvar_cnvs_2020-05.bed -name clinvar_cnvs -no_duplicates -url_decode -out {$som_clincnv}", true);

		$hgmd_file = "{$data_folder}/dbs/HGMD/HGMD_CNVS_2020_1.bed"; //optional because of license
		if (file_exists($hgmd_file))
		{
			$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$som_clincnv} -in2 {$hgmd_file} -name hgmd_cnvs -no_duplicates -url_decode -out {$som_clincnv}", true);
		}
		$omim_file = "{$data_folder}/dbs/OMIM/omim.bed"; //optional because of license
		if (file_exists($omim_file))
		{
			$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$som_clincnv} -in2 {$omim_file} -no_duplicates -url_decode -out {$som_clincnv}", true);
		}
		$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$som_clincnv} -in2 {$repository_basedir}/data/gene_lists/genes.bed -no_duplicates -url_decode -out {$som_clincnv}", true);

		//annotate additional gene info
		$parser->exec(get_path("ngs-bits")."CnvGeneAnnotation", "-in {$som_clincnv} -out {$som_clincnv}", true);
		// skip annotation if no connection to the NGSD is possible
		if (db_is_enabled("NGSD"))
		{
			//annotate overlap with pathogenic CNVs
			$parser->exec(get_path("ngs-bits")."NGSDAnnotateCNV", "-in {$som_clincnv} -out {$som_clincnv}", true);
		}
	}*/
	if($single_sample) //use clincnv_germline with -mosaicism flag in case of single samples
	{		
		//create BED file with GC and gene annotations - if missing
		$bed = $ref_folder_t."/roi_annotated.bed";
		if (!file_exists($bed))
		{
			$ngsbits = get_path("ngs-bits");
			$pipeline = [
					["{$ngsbits}BedAnnotateGC", "-in ".$sys['target_file']." -ref ".genome_fasta($sys['build'])],
					["{$ngsbits}BedAnnotateGenes", "-out {$bed}"],
				];
			$parser->execPipeline($pipeline, "creating annotated BED file for ClinCNV");
		}
		//perform CNV analysis
		$cnv_out = $tmp_folder."/output.tsv";
		$cnv_out2 = $tmp_folder."/output.seg";
		$args = array(
			"-cov {$t_cov}",
			"-cov_folder {$ref_folder_t}",
			"-bed {$bed}",
			"-out {$som_clincnv}",
			"--log ".$parser->getLogFile(),
			"-mosaicism",
			"-cov_max 200",
			"-max_cnvs 200"
		);
		$parser->execTool("NGS/vc_clincnv_germline.php", implode(" ", $args), true);

		// annotate CNV file
		$repository_basedir = repository_basedir();
		$data_folder = get_path("data_folder");
		$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$som_clincnv} -in2 {$repository_basedir}/data/misc/cn_pathogenic.bed -no_duplicates -url_decode -out {$som_clincnv}", true);
		$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$som_clincnv} -in2 {$data_folder}/dbs/ClinGen/dosage_sensitive_disease_genes.bed -no_duplicates -url_decode -out {$som_clincnv}", true);
		$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$som_clincnv} -in2 {$data_folder}/dbs/ClinVar/clinvar_cnvs_2020-05.bed -name clinvar_cnvs -no_duplicates -url_decode -out {$som_clincnv}", true);

		$hgmd_file = "{$data_folder}/dbs/HGMD/HGMD_CNVS_2020_1.bed"; //optional because of license
		if (file_exists($hgmd_file))
		{
			$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$som_clincnv} -in2 {$hgmd_file} -name hgmd_cnvs -no_duplicates -url_decode -out {$som_clincnv}", true);
		}
		$omim_file = "{$data_folder}/dbs/OMIM/omim.bed"; //optional because of license
		if (file_exists($omim_file))
		{
			$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$som_clincnv} -in2 {$omim_file} -no_duplicates -url_decode -out {$som_clincnv}", true);
		}
		$parser->exec(get_path("ngs-bits")."BedAnnotateFromBed", "-in {$som_clincnv} -in2 {$repository_basedir}/data/gene_lists/genes.bed -no_duplicates -url_decode -out {$som_clincnv}", true);

		//annotate additional gene info
		$parser->exec(get_path("ngs-bits")."CnvGeneAnnotation", "-in {$som_clincnv} -out {$som_clincnv}", true);
		// skip annotation if no connection to the NGSD is possible
		if (db_is_enabled("NGSD"))
		{
			//annotate overlap with pathogenic CNVs
			$parser->exec(get_path("ngs-bits")."NGSDAnnotateCNV", "-in {$som_clincnv} -out {$som_clincnv}", true);
		}
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
				$n_cov_off_target = $ref_file_n_off_target;
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
			$tmp_file_name = $parser->tempFile(".csv");
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

		//Skip CNV Calling if there are less than 7 tumor-normal coverage pairs
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
			
			if(isset($cnv_baseline_pos))
			{
				$args_clincnv[] = "-guide_baseline $cnv_baseline_pos";
			}
			
			$parser->execTool("NGS/vc_clincnv_somatic.php",implode(" ",$args_clincnv));
		}
		else
		{
			trigger_error("Not enough reference tumor-normal coverage files for processing system {$system} found. Skipping CNV calling.\n", E_USER_NOTICE);
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
	$parser->execTool("Pipelines/annotate.php", "-out_name $prefix -out_folder $tmp_folder1 -system $system -vcf $variants -somatic -updown -threads $threads");

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
			"-target_exons", repository_basedir()."/data/gene_lists/gene_exons.bed", //file containing all human exons to determine exonic variants in TMB calculation
			"-blacklist", repository_basedir() ."/data/gene_lists/somatic_tmb_blacklist.bed", //Blacklisted genes that are not included in TMB calculation (e.g. HLA-A and HLA-B)
			"-tsg_bed", repository_basedir() ."/data/gene_lists/somatic_tmb_tsg.bed", //TSG genes whose mutations are treated specially in TMB calculation
			"-ref", $ref_genome,
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
			
			//If user is not in NGSD, we use the user "unknown"
			if($user_id == -1) $user_id = $db_conn->getValue("SELECT id FROM user WHERE user_id='unknown'");
			
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
	$tmp_cgi_zipped_outfile = $parser->tempFile(".zip");
	$parameters = [
		"-out",$tmp_cgi_zipped_outfile,
		"-cancertype", $cancer_type
	];
	if (file_exists($variants_annotated)) $parameters[] = "-mutations $variants_gsvar";
	if (file_exists($som_clincnv)) $parameters[] = "-cnas $som_clincnv";
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
		$parser->execTool("NGS/an_somatic_gsvar.php","-gsvar_in $variants_gsvar -cgi_snv_in $cgi_snv_result_file -out $variants_gsvar -include_ncg");
	}
	//annotate CGI cnv genes to cnv input file (which was originally created by CNVHunter)
	$cgi_cnv_result_file = $full_prefix . "_cgi_cnv_analysis.tsv";
	if(file_exists($som_clincnv) && file_exists($cgi_cnv_result_file))
	{
		$parser->execTool("NGS/an_somatic_cnvs.php","-cnv_in $som_clincnv -cnv_in_cgi $cgi_cnv_result_file -out $som_clincnv -include_ncg -include_cytoband");
	}
	
	/*************************
	 * GERMLINE CGI ANALYSIS *
	 *************************/
	$variants_germline = dirname($n_bam)."/{$n_id}.GSvar";
	if(file_exists($variants_germline) && $include_germline)
	{
		$germl_gsvar_content = file($variants_germline, FILE_IGNORE_NEW_LINES);
		$tmp_new_file = array();
		
		//Create temporary GSvar file filtered for variants with < 3% AF in gnomAD for CGI annotation
		$i_gnomad = -1;
		$skip_germline_snvs = false;
		foreach($germl_gsvar_content as $line)
		{
			if(starts_with($line, "#chr") )
			{
				$parts = explode("\t",$line);
				for($i=0; $i<count($parts); ++$i)
				{
					if($parts[$i] == "gnomAD")
					{
						$i_gnomad = $i;
						break;
					}
				}
				$tmp_new_file[] = $line;
				continue;
			}
			if(starts_with($line,"#"))
			{
				$tmp_new_file[]= $line;
				continue;
			}
			
			if($i_gnomad < -1)
			{
				trigger_error("Could not find column \"gnomAD\" for filtering germline variants for CGI annotation. Skipping germline SNV CGI annotation.", E_USER_WARNING);
				$skip_germline_snvs = true;
				break;
			}
			
			$parts = explode("\t",$line);
			
			$af_gnomad = $parts[$i_gnomad];
			if(!empty($af_gnomad) && floatval($af_gnomad) > 0.03) continue;
			
			$tmp_new_file[] = $line;
		}
		$filtered_germline_variants = $parser->tempFile(".GSvar");
		file_put_contents($filtered_germline_variants, implode("\n",$tmp_new_file) );
		
		$tmp_cgi_zipped_outfile = $parser->tempFile();
		
		$params = [
			"-out",$tmp_cgi_zipped_outfile,
			"-cancertype","CANCER",
			"-no_snv_limit"
		];
		
		if(!$skip_germline_snvs) $params[] = "-mutations $filtered_germline_variants";

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
			$parser->execTool("NGS/an_somatic_gsvar.php",$parameters,false);
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

//RNA annotation
if (in_array("an_rna", $steps))
{
	list($t_name) = explode("_", $t_id);
	
	//Determine tumor rna sample from NGSD via sample_relations
	$ps_rna_bams = array();
	
	if( db_is_enabled("NGSD") && !isset($t_rna_bam) )
	{
		$db = DB::getInstance("NGSD");
		$res = $db->executeQuery("SELECT id, name FROM sample WHERE name=:name", array("name" => $t_name));
		$sample_id = $res[0]["id"];
		
		//get related RNA samples from NGSD
		if (count($res) >= 1)
		{
			$res = $db->executeQuery("SELECT * FROM sample_relations WHERE relation='same sample' AND (sample1_id=:sid OR sample2_id=:sid)", array("sid" => $sample_id));
			foreach($res as $row)
			{
				$sample_id_annotation = $row['sample1_id'] != $sample_id ? $row['sample1_id'] : $row['sample2_id'];
				
				$res = $db->executeQuery("SELECT ps.sample_id, ps.process_id, ps.processing_system_id, ps.quality, sys.id, sys.type, CONCAT(s.name, '_', LPAD(ps.process_id, 2, '0')) as psample FROM processed_sample as ps, processing_system as sys, sample as s WHERE ps.sample_id=:sid AND sys.type='RNA' AND ps.processing_system_id=sys.id AND ps.sample_id=s.id AND (NOT ps.quality='bad')", array("sid" => $sample_id_annotation));
				
				$rna_ids = array_column($res, 'psample');
				foreach($rna_ids as $rna_id)
				{
					$ps_rna_bams[$rna_id] = get_processed_sample_info($db, $rna_id)["ps_bam"];
				}
			}
		}
	}
	elseif(isset($t_rna_bam))
	{
		$ps_rna_bams[basename($t_rna_bam)] = $t_rna_bam; 
	}
	
	if(count($ps_rna_bams) < 1)
	{
		trigger_error("For annotation step \"an_rna\" tumor RNA bam file must be specified (via paramter -t_rna_bam or determined via sample_relations).", E_USER_ERROR);
	}
	
	//Determine reference tissue type 1.) from parameter -rna_ref_tissue or 2.) from NGSD (if available)
	if(!isset($rna_ref_tissue) && !db_is_enabled("NGSD"))
	{
		trigger_error("For annotation step \"an_rna\" a tissue type for RNA reference data has to be specified.", E_USER_ERROR);
	}
	if(!isset($rna_ref_tissue) && db_is_enabled("NGSD"))
	{
		$db = DB::getInstance("NGSD");
		$res = $db->getValues("SELECT di.disease_info FROM sample as s, sample_disease_info as di WHERE s.name = '{$t_name}' AND s.id = di.sample_id AND di.type = 'RNA reference tissue'");
		if(count($res) == 1)
		{
			list($rna_ref_tissue) = $res;
		}
		else
		{
			trigger_error("Found multiple or no RNA reference tissue in NGSD. Aborting...", E_USER_ERROR);
		}
	}
	
	$rna_counts = array(); //file contains transcript counts
	foreach($ps_rna_bams as $rna_id => $rna_bam)
	{
		$rna_counts_tmp = glob(dirname($rna_bam)."/*_counts.tsv");
		if(count($rna_counts_tmp) != 1)
		{
			trigger_error("Could not find or found multiple RNA count files in sample folder dirname($rna_bam)", E_USER_ERROR);
		}
		
		$rna_counts[$rna_id] = $rna_counts_tmp[0];
	}

	//Annotate data from all detected RNA files
	foreach($ps_rna_bams as $rna_id => $rna_bam)
	{
		$rna_count = $rna_counts[$rna_id];
		
		//SNVs
		$args = [
		"-gsvar_in $variants_gsvar",
		"-out $variants_gsvar",
		"-rna_id $rna_id",
		"-rna_counts $rna_count",
		"-rna_bam $rna_bam"];
		$parser->execTool("NGS/an_somatic_gsvar.php", implode(" ", $args));
		
		//CNVs
		if(file_exists($som_clincnv))
		{
			$parser->execTool("NGS/an_somatic_cnvs.php", " -cnv_in $som_clincnv -out $som_clincnv -rna_counts $rna_count -rna_id $rna_id -rna_ref_tissue " .str_replace(" ", 0, $rna_ref_tissue));
		}
	}
	
	//Reference tissue SNVs
	$args = [
		"-gsvar_in $variants_gsvar",
		"-out $variants_gsvar",
		"-rna_ref_tissue " .str_replace(" ", 0, $rna_ref_tissue)//Replace spaces by 0 because it is diffcult to pass spaces via command line.
	];
	$parser->execTool("NGS/an_somatic_gsvar.php", implode(" ", $args));
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
			// check sex using control sample
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

			
			// import SNVS (not for WGS)
			if (file_exists($variants_gsvar) && $sys['type'] !== "WGS")
			{
				
				$parser->exec(get_path("ngs-bits") . "/NGSDAddVariantsSomatic", " -t_ps $t_id -n_ps $n_id -var $variants_gsvar -var_force");
			}
			
			if(file_exists($som_clincnv) && $sys['type'] !== "WGS")
			{
				$parser->exec(get_path("ngs-bits") . "/NGSDAddVariantsSomatic", " -t_ps $t_id -n_ps $n_id -cnv $som_clincnv -cnv_force");
			}
		}
	}
}
?>
