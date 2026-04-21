<?php

/**
	@page multisample_longread
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("multisample_longread", "Multisample analysis pipeline for ONT/PacBio long-read data.");
$parser->addInfileArray("bams", "Input BAM/CRAM files.", false);
$parser->addStringArray("status", "List of affected status of the input samples (BAMs) - can be 'affected' or 'control'.", false);
$parser->addInfile("out_folder", "Output folder name.", false);
//optional
$parser->addString("prefix", "Output file prefix.", true, "multi");
$parser->addInfile("system",  "Processing system INI file used for all samples (automatically determined from NGSD if the basename of 'c' is a valid processed sample name).", true);
$steps_all = array("vc", "cn", "sv", "an", "db");
$parser->addString("steps", "Comma-separated list of steps to perform:\nvc=variant calling, cn=copy-number analysis, sv=structural variant calling, an=annotation, db=database import.", true, implode(",", $steps_all));
$parser->addInt("threads", "The maximum number of threads used.", true, 2);
$parser->addFlag("no_sync", "Skip syncing annotation databases and genomes to the local tmp folder (Needed only when starting many short-running jobs in parallel).");
$parser->addInfile("ped", "(PLINK-compatible) Pedigree file to phase related samples.", true);
$parser->addFlag("no_splice", "Skip SpliceAI scoring of variants that are not precalculated.");
extract($parser->parse($argv));

//create output folder if missing
$out_folder .= "/";
if (!file_exists($out_folder))
{
	if (!mkdir($out_folder))
	{
		trigger_error("Could not create output folder '{$out_folder}'!", E_USER_ERROR);
	}
	if (!chmod($out_folder, 0777))
	{
		trigger_error("Could not change privileges of output folder '{$out_folder}'!", E_USER_ERROR);
	}
}

//create log file in output folder if none is provided
if ($parser->getLogFile()=="") $parser->setLogFile($out_folder."/multi_".date("YmdHis").".log");

//log server, user, etc.
$parser->logServerEnvronment();

//check steps
$steps = explode(",", $steps);
foreach($steps as $step)
{
	if (!in_array($step, $steps_all)) trigger_error("Unknown processing step '$step'!", E_USER_ERROR);
}

//file names
$gsvar = "{$out_folder}/{$prefix}.GSvar";
$vcf_file = "{$out_folder}/{$prefix}_var.vcf.gz";
$cnv_file = "{$out_folder}/{$prefix}_cnvs_clincnv.tsv";
$sv_vcf_file = "{$out_folder}/{$prefix}_var_structural_variants.vcf.gz";
$bedpe_out = "{$out_folder}/{$prefix}_var_structural_variants.bedpe";

//check input counts
if (count($bams)!=count($status))
{
	trigger_error("Different number of arguments for 'bams' and 'status' parameters! BAMs: ".count($bams)." Status: ".count($status)."!", E_USER_ERROR);
}

//check input status
$counts = array_count_values($status);
if (isset($counts["affected"]) && $counts["affected"] === 1)
{
	$affected_bam = $bams[array_search("affected", $status)];
}
else $affected_bam="";

foreach($status as $stat)
{
	$valid = array("affected", "control");
	if (!in_array($stat, $valid))
	{
		trigger_error("Invalid status '$stat' given for '$bam'. Valid are '".implode("', '", $valid)."'", E_USER_ERROR);
	}
}
$status = array_combine($bams, $status);

//make BAMs absolute paths
for ($i=0; $i<count($bams); ++$i)
{
	$bams[$i] = realpath($bams[$i]);
}

//check input sample names
$names = [];
foreach($bams as $bam)
{
	$name = basename2($bam);
	if (in_array($name, array_values($names))) trigger_error("Sample file name '$name' occurs twice in input file list. Each sample must be uniquely identifiable by name!", E_USER_ERROR);
	$names[$bam] = $name;
}

//determine system/platform (from first sample)
$sys = load_system($system, $names[$bams[0]]);
$genome = genome_fasta($sys['build']);
$platform = $sys['platform'];
trigger_error("Long-read platform: {$platform}", E_USER_NOTICE);
if ($platform!="PacBio" && $platform!="ONT") trigger_error("Unsupported platform '{$platform}'! Only 'ONT' or 'PacBio' are supported!", E_USER_ERROR);
if ($sys['target_file']=="") trigger_error("Unsupported processing system '".$sys2["name_short"]."' (has no target region)!", E_USER_ERROR);
	
//check genome build and processing systems of all samples
foreach($bams as $bam)
{
	check_genome_build($bam, $sys['build']);
	$sys2 = load_system($system, $names[$bam]);
	if($sys2['type']!="lrGS") trigger_error("Unsupported system type '".$sys2['type']."' of {$bam} (only 'lrGS' is supported)!", E_USER_ERROR);
	if($sys2['platform']!=$platform) trigger_error("Unsupported mixed platforms '{$platform}' and '".$sys2['platform']."'!", E_USER_ERROR);
}

//enable check if chromosomes are complete only if the ROI covers the whole genome (pipeline test are not)
$check_chrs = bed_size(realpath($sys['target_file'])) > 3e9;

//set up local NGS data copy (to reduce network traffic and speed up analysis)
if (!$no_sync) $parser->execTool("Tools/data_setup.php", "-build ".$sys['build']);

//variant calling (and phasing)
if (in_array("vc", $steps))
{
	//determine gVCFs	
	$gvcfs = [];
	foreach($bams as $bam)
	{
		$gvcf = dirname($bam)."/".basename2($bam)."_var.gvcf.gz";
		if (!file_exists($gvcf)) trigger_error("gVCF file '{$gvcf}' not found! Please make sure the single-sample variant calling has run successfully.", E_USER_ERROR);
		$gvcfs[] = $gvcf;
	}
	
	// merge GVCFs
	$analysis_type = ($prefix=="trio" ? "GERMLINE_TRIO" : "GERMLINE_MULTISAMPLE");
	if($platform=="ONT")
	{
		$args = [];
		$args[] = "-gvcfs ".implode(" ", $gvcfs);
		$args[] = "-status ".implode(" ", $status);
		$args[] = "-out {$vcf_file}";
		$args[] = "-threads {$threads}";
		$args[] = "-analysis_type {$analysis_type}";
		$args[] = "-mode clair3";
		$parser->execTool("Tools/merge_gvcf.php", implode(" ", $args));
	}
	else //PacBio
	{
		//merge gVCFs with GLnexus
		$tmp_vcf_merged = $parser->tempFile(".vcf.gz");
		$args = [];
		$args[] = "--dir ".$parser->tempFolder()."/GLnexus.DB/";
		$args[] = "--config DeepVariantWGS";
		$args[] = "--threads {$threads}";
		$args[] = implode(" ", $gvcfs);
		$pipeline = [];
		$pipeline[] = ["", $parser->execApptainer("glnexus", "glnexus_cli", implode(" ", $args), $gvcfs, [], true)];
		$pipeline[] = ["", $parser->execApptainer("glnexus", "bcftools view", "", [], [], true)];
		$pipeline[] = ["", $parser->execApptainer("htslib", "bgzip", "-@ -c > {$tmp_vcf_merged}", [], [], true)];
		$parser->execPipeline($pipeline, "GLnexus gVCF merging");
		
		//post-processing pipeline
		$tmp_vcf_post = $parser->tempFile(".vcf");
		$pipeline = [];
		$pipeline[] = array("zcat", $tmp_vcf_merged);
		//split complex variants to primitives - this step has to be performed before VcfBreakMulti - otherwise mulitallelic variants that contain both 'hom' and 'het' genotypes fail - see NA12878 amplicon test chr2:215632236-215632276
		$pipeline[] = ["", $parser->execApptainer("vcflib", "vcfallelicprimitives", "-kg", [], [], true)];
		//split multi-allelic variants
		$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfBreakMulti", "-no_errors", [], [], true)];
		//remove invalid variants
		$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfFilter", "-remove_invalid -ref {$genome}", [$genome], [], true)];
		//normalize all variants and align INDELs to the left
		$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfLeftNormalize", "-stream -ref {$genome}", [$genome], [], true)];
		//sort variants by genomic position
		$pipeline[] = ["", $parser->execApptainer("ngs-bits", "VcfStreamSort", "> {$tmp_vcf_post}", [], [], true)];
		//execute post-processing pipeline
		$parser->execPipeline($pipeline, "DeepVariant multi-sample post processing");
		
		//add pipeline, sample infos, etc. to VCF header
		list($comments, $samples) = vcf_load_header($tmp_vcf_post);
		$comments[] = "##reference={$genome}\n";
		$comments[] = "##ANALYSISTYPE={$analysis_type}\n";
		$comments[] = "##PIPELINE=".repository_revision(true)."\n";
		for ($i=0; $i < count($gvcfs); $i++) 
		{ 
			$comments[] = gsvar_sample_header($samples[$i], array("DiseaseStatus"=>array_values($status)[$i]), "##", "");
		}
		vcf_replace_comments($tmp_vcf_post, $comments);
		
		//bgzip and index
		$parser->execApptainer("htslib", "bgzip", "-@ {$threads} -c {$tmp_vcf_post} > {$vcf_file}", [], [dirname($vcf_file)]);
		$parser->execApptainer("htslib", "tabix", "-f -p vcf {$vcf_file}", [], [dirname($vcf_file)]);
	}

	//phasing (WhatsHap)
	$vcf_file_phased = $parser->tempFile("_phased.vcf");
	$args = [];
	$args[] = "phase";
	if ($ped != "")	$args[] = "--ped {$ped}";
	$args[] = "--reference={$genome}";
	$args[] = "-o {$vcf_file_phased}";
	$args[] = "{$vcf_file}";
	$args[] = implode(" ", $bams);
	$out_files = [];
	$out_files[] = dirname($vcf_file);
	foreach($bams as $bam)
	{
		$out_files[] = dirname($bam);
	}
	$parser->execApptainer("whatshap", "whatshap", implode(" ", $args), [$ped, $genome], $out_files);

	//create compressed file and index and replace original VCF
	$parser->execApptainer("htslib", "bgzip", "-c {$vcf_file_phased} > {$vcf_file}", [], [dirname($vcf_file)]);
	$parser->execApptainer("htslib", "tabix", "-f -p vcf $vcf_file", [], [dirname($vcf_file)]);
}

//copy-number
if (in_array("cn", $steps))
{
	//prepare multi-sample parameters
	$args_multisample = [
		"-bams ".implode(" ", $bams),
		"-status ".implode(" ", $status),
		"-out_folder $out_folder",
		"-system $system",
		"-prefix {$prefix}",
		"-threads $threads",
		"-no_sync" //already done if needed
		];
	$parser->execTool("Pipelines/multisample.php", implode(" ", $args_multisample)." -steps cn", true);
}

//sv calling
if (in_array("sv", $steps))
{
	//run Sniffles
	$parser->execTool("Tools/vc_sniffles.php", "-bam  ".implode(" ", $bams)." -include_mosaic -out {$sv_vcf_file} -threads {$threads} -build ".$sys['build']);
}

//annotation
if (in_array("an", $steps))
{
	$status_map = [];
	foreach ($status as $bam => $disease_status) 
	{
		$status_map[basename2($bam)] = $disease_status;
	}

	//annotate small variant VCF file
	if (file_exists($vcf_file))
	{
		//basic annotation
		$args = [];
		$args[] = "-out_name {$prefix}";
		$args[] = "-out_folder {$out_folder}";
		$args[] = "-system {$system}";
		$args[] = "-threads {$threads}";
		$args[] = "-multi";
		if ($no_splice) $args[] = "-no_splice";
		$parser->execTool("Tools/annotate.php", "-out_name $prefix -out_folder $out_folder -system $system -threads $threads -multi");
		
		//check for truncated output
		if ($check_chrs) $parser->execTool("Tools/check_for_missing_chromosomes.php", "-in {$out_folder}/{$prefix}_var_annotated.vcf.gz -max_missing_perc 5");
		
		//update sample entry 
		update_gsvar_sample_header($gsvar, $status_map);
	}
	else
	{
		trigger_error("Small variant file {$vcf_file} does not exist, skipping small variant annotation!", E_USER_WARNING);
	}

	// annotate CNV file
	if (file_exists($cnv_file))
	{
		$repository_basedir = repository_basedir();
		$data_folder = get_path("data_folder");
		$parser->execApptainer("ngs-bits", "BedAnnotateFromBed", "-in {$cnv_file} -in2 {$repository_basedir}/data/misc/af_genomes_imgag.bed -overlap -out {$cnv_file}", [$cnv_file, "{$repository_basedir}/data/misc/af_genomes_imgag.bed"]);
		$parser->execApptainer("ngs-bits", "BedAnnotateFromBed", "-in {$cnv_file} -in2 {$repository_basedir}/data/misc/cn_pathogenic.bed -no_duplicates -url_decode -out {$cnv_file}", [$cnv_file, "{$repository_basedir}/data/misc/cn_pathogenic.bed"]);
		$parser->execApptainer("ngs-bits", "BedAnnotateFromBed", "-in {$cnv_file} -in2 {$data_folder}/dbs/ClinGen/dosage_sensitive_disease_genes_GRCh38.bed -no_duplicates -url_decode -out {$cnv_file}", [$cnv_file, "{$data_folder}/dbs/ClinGen/dosage_sensitive_disease_genes_GRCh38.bed"]);
		$parser->execApptainer("ngs-bits", "BedAnnotateFromBed", "-in {$cnv_file} -in2 {$data_folder}/dbs/ClinVar/clinvar_cnvs_2026-04.bed -name clinvar_cnvs -no_duplicates -url_decode -out {$cnv_file}", [$cnv_file, "{$data_folder}/dbs/ClinVar/clinvar_cnvs_2026-04.bed"]);


		$hgmd_file = "{$data_folder}/dbs/HGMD/HGMD_CNVS_2026_1.bed"; //optional because of license
		if (file_exists($hgmd_file))
		{
			$parser->execApptainer("ngs-bits", "BedAnnotateFromBed", "-in {$cnv_file} -in2 {$hgmd_file} -name hgmd_cnvs -no_duplicates -url_decode -out {$cnv_file}", [$cnv_file, $hgmd_file]);
		}
		$omim_file = "{$data_folder}/dbs/OMIM/omim.bed"; //optional because of license
		if (file_exists($omim_file))
		{
			$parser->execApptainer("ngs-bits", "BedAnnotateFromBed", "-in {$cnv_file} -in2 {$omim_file} -no_duplicates -url_decode -out {$cnv_file}", [$cnv_file, $omim_file]);
		}

		//annotate gene info
		if(db_is_enabled("NGSD"))
		{
			$parser->execApptainer("ngs-bits", "CnvGeneAnnotation", "-in {$cnv_file} -add_simple_gene_names -out {$cnv_file}", [$cnv_file]);
		}
		
		//annotate overlap with pathogenic CNVs
		if (db_is_enabled("NGSD"))
		{
			$parser->execApptainer("ngs-bits", "NGSDAnnotateCNV", "-in {$cnv_file} -out {$cnv_file}", [$cnv_file]);
		}
	}
	else
	{
		trigger_error("CNV file {$cnv_file} does not exist, skipping CNV annotation!", E_USER_WARNING);
	}

	//annotate SV VCF file
	if (file_exists($sv_vcf_file))
	{
		check_genome_build($sv_vcf_file, $sys['build']);

		//create BEDPE files
		$parser->execApptainer("ngs-bits", "VcfToBedpe", "-in {$sv_vcf_file} -out {$bedpe_out}", [$sv_vcf_file], [dirname($bedpe_out)]);

		// correct filetype
		$bedpe_table = Matrix::fromTSV($bedpe_out);
		$bedpe_table->removeComment("#fileformat=BEDPE");
		if ($prefix == "trio") $bedpe_table->prependComment("#fileformat=BEDPE_GERMLINE_TRIO");
		else $bedpe_table->prependComment("#fileformat=BEDPE_GERMLINE_MULTI");
		$bedpe_table->toTSV($bedpe_out);

		//add gene info annotation
		if (db_is_enabled("NGSD"))
		{
			$parser->execApptainer("ngs-bits", "BedpeGeneAnnotation", "-in {$bedpe_out} -out {$bedpe_out} -add_simple_gene_names", [$bedpe_out]);
		}

		//add NGSD counts from flat file
		$ngsd_annotation_folder = get_path("data_folder")."/dbs/NGSD/";
		$ngsd_sv_files = array("sv_deletion.bedpe.gz", "sv_duplication.bedpe.gz", "sv_insertion.bedpe.gz", "sv_inversion.bedpe.gz", "sv_translocation.bedpe.gz");
		$db_file_dates = [];

		// check file existance
		$all_files_available = file_exists($ngsd_annotation_folder."sv_breakpoint_density.igv");
		foreach ($ngsd_sv_files as $filename) 
		{
			if(!(file_exists($ngsd_annotation_folder.$filename) && file_exists($ngsd_annotation_folder.$filename.".tbi")))
			{
				$all_files_available = false;
				break;
			}
		}
		if ($all_files_available)
		{
			// store flat file modification date to detect changes during annotation 
			foreach ($ngsd_sv_files as $filename)
			{
				$db_file_dates[$filename] = filemtime($ngsd_annotation_folder.$filename);
				if ($db_file_dates[$filename] == false)
				{
					trigger_error("Cannot get modification date of '".$ngsd_annotation_folder.$filename."'!",E_USER_ERROR);
				}
			}
			
			//perform annotation
			$args = array(
				"-in $bedpe_out",
				"-out $bedpe_out",
				"-processing_system ".$sys["name_short"],
				"-ann_folder $ngsd_annotation_folder",
			);

			if (db_is_enabled("NGSD") && $affected_bam!="")
			{
				$db = DB::getInstance("NGSD", false);
				$ps_id = get_processed_sample_id($db, $names[$affected_bam], false);
	
				if ($ps_id != -1)
				{
					$args[] = "-ps_name ".$names[$affected_bam];
				}
				else 
				{
					trigger_error("No processed sample ID found for sample ".$names[$affected_bam].", skipping count annotation by disease group!", E_USER_WARNING);
				}
			}
			
			$parser->execApptainer("ngs-bits", "BedpeAnnotateCounts", implode(" ", $args), [$bedpe_out, $ngsd_annotation_folder]);

			$sys_specific_density_file = $ngsd_annotation_folder."sv_breakpoint_density_".$sys["name_short"].".igv";
			if (file_exists($sys_specific_density_file))
			{
				$parser->execApptainer("ngs-bits", "BedpeAnnotateBreakpointDensity", "-in {$bedpe_out} -out {$bedpe_out} -density {$ngsd_annotation_folder}sv_breakpoint_density.igv -density_sys {$sys_specific_density_file}", [$bedpe_out, $ngsd_annotation_folder, $sys_specific_density_file]);
			}
			else
			{
				$parser->execApptainer("ngs-bits", "BedpeAnnotateBreakpointDensity", "-in {$bedpe_out} -out {$bedpe_out} -density {$ngsd_annotation_folder}sv_breakpoint_density.igv", [$bedpe_out, $ngsd_annotation_folder]);
			}

			// check if files changed during annotation
			foreach ($ngsd_sv_files as $filename)
			{
				if ($db_file_dates[$filename] != filemtime($ngsd_annotation_folder.$filename))
				{
					trigger_error("Annotation file '".$ngsd_annotation_folder.$filename."' has changed during annotation!",E_USER_ERROR);
				}
			}

			//annotate class 4 and 5 pathogenic SVs			
			if (db_is_enabled("NGSD"))
			{
				$parser->execApptainer("ngs-bits", "NGSDAnnotateSV", "-in {$bedpe_out} -out {$bedpe_out}", [$bedpe_out]);
			}
		}
		else
		{
			trigger_error("Cannot annotate NGSD counts! At least one required file in '{$ngsd_annotation_folder}' is missing!", E_USER_WARNING);
		}
		
		//add optional OMIM annotation
		$omim_file = get_path("data_folder")."/dbs/OMIM/omim.bed"; 
		if(file_exists($omim_file))//OMIM annotation (optional because of license)
		{
			$parser->execApptainer("ngs-bits", "BedpeAnnotateFromBed", "-in $bedpe_out -out $bedpe_out -bed $omim_file -url_decode -replace_underscore -col_name OMIM", [$bedpe_out, $omim_file]);
		}

		//add CNV overlap annotation
		if (file_exists($cnv_file))
		{
			$parser->execApptainer("ngs-bits", "BedpeAnnotateCnvOverlap", "-in $bedpe_out -out $bedpe_out -cnv $cnv_file", [$bedpe_out, $cnv_file]);
		}

		//write genotype in own column
		$parser->execApptainer("ngs-bits", "BedpeExtractGenotype", "-in $bedpe_out -out $bedpe_out -include_unphased", [$bedpe_out]);

		//extract columns
		$parser->execApptainer("ngs-bits", "BedpeExtractInfoField", "-in $bedpe_out -out $bedpe_out -info_fields SVLEN,SUPPORT,COVERAGE", [$bedpe_out]);

		//update sample entry 
		update_gsvar_sample_header($bedpe_out, $status_map);
		
		//check for truncated output
		if ($check_chrs) $parser->execTool("Tools/check_for_missing_chromosomes.php", "-in {$bedpe_out}");
	}

}

//everything worked > import/update sample data in NGSD
if (in_array("db", $steps) && db_is_enabled("NGSD"))
{
	$db_conn = DB::getInstance("NGSD", false);
	
	//add secondary analysis (if missing)
	$parser->execTool("Tools/db_import_secondary_analysis.php", "-type 'multi sample' -gsvar {$gsvar}");
}

?>
