<?php

/** 
	@page vc_dragen_somatic
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_dragen_somatic", "Performs small/structural variant calling for tumor or tumor-normal samples using the Illumina DRAGEN server.");
$parser->addInfile("t_bam", "Input file tumor sample BAM.", false);
$parser->addOutfile("out", "Output .VCF file for small variants.", false);

//optional
$parser->addInfile("n_bam", "Input file normal sample BAM.", true);
$parser->addOutfile("out_sv", "Outfile for dragen somatic SV calls. Activates dragen SV calling when given.", true);
$parser->addString("build", "The genome build to use. The genome must be indexed for BWA!", true, "GRCh38");
$parser->addString("tumor", "Sample name of the tumor sample. If unset the basename of the 'tumor_bam' file is used.", true, "");
$parser->addString("normal", "Sample name of the normal Sample. If unset the basename of the 'normal_bam' file is used.", true, "");
$parser->addFlag("is_targeted", "If the sequencing is a targeted sequencing.");
$parser->addFlag("debug", "Add debug output to the log file.");

extract($parser->parse($argv));

// ********************************* init *********************************//

if ($debug)
{
	list($stdout) = exec2("hostname");
	$parser->log("server:", $stdout);
	
	list($stdout) = exec2("whoami");
	$parser->log("user:", $stdout);
	
	list($stdout) = exec2("ulimit -a");
	$parser->log("limits:", $stdout);
	
	list($stdout) = exec2("env");
	$parser->log("env:", $stdout);

	list($stdout) = exec2("dragen_lic");
	$parser->log("dragen licence status before:", $stdout);
}


//if no sample name is given use output name
if ($tumor=="") $tumor = basename($t_bam, ".bam");
if ($n_bam != "" && $normal=="") $normal = basename($n_bam, ".bam");

//check if valid reference genome is provided
$dragen_genome_path = get_path("dragen_genomes")."/".$build."/dragen/";
if (!file_exists($dragen_genome_path)) 
{
	trigger_error("Invalid genome build '".$build."' given. Path '".$dragen_genome_path."' not found on Dragen!", E_USER_ERROR);
}

//check if input files are readable and output file is writeable
if (!is_readable($t_bam)) trigger_error("Input file '$t_bam' is not readable!", E_USER_ERROR);
if ($n_bam != "" && !is_readable($n_bam)) trigger_error("Input file '$n_bam' is not readable!", E_USER_ERROR);
if (!is_writable2($out)) trigger_error("Output file '$out' is not writable!", E_USER_ERROR);
if ($out_sv != "" && !is_writable2($out_sv)) trigger_error("Output file '$out_sv' is not writable!", E_USER_ERROR);

// create empty folder for analysis
$working_dir = get_path("dragen_data")."/megSAP_working_dir/";
if (file_exists($working_dir))	
{
	$parser->exec("rm", "-rf $working_dir");
}
if (!mkdir($working_dir, 0700))
{
	trigger_error("Could not create working directory '".$working_dir."'!", E_USER_ERROR);
}

// ********************************* call dragen *********************************//

//parameters
$dragen_parameter = [];
$dragen_parameter[] = "-r ".$dragen_genome_path;
$dragen_parameter[] = "--tumor-bam-input ".$t_bam;
if ($n_bam != "") $dragen_parameter[] = "--bam-input ".$n_bam;
$dragen_parameter[] = "--output-directory $working_dir";
$dragen_parameter[] = "--output-file-prefix output";
$dragen_parameter[] = "--enable-map-align false"; # cannot map multiple (tumor, normal) inputs at once
$dragen_parameter[] = "--pair-by-name true";

//small variant calling
$dragen_parameter[] = "--enable-variant-caller true";
$dragen_parameter[] = "--vc-min-read-qual 1";
$dragen_parameter[] = "--vc-min-tumor-read-qual 3"; #default 3 for t-n, 20 for t-only 
$dragen_parameter[] = "--vc-min-base-qual 15";
$dragen_parameter[] = "--vc-callability-tumor-thresh 15"; # default 15 - minimum coverage in bam to try variant calling at that position
$dragen_parameter[] = "--vc-callability-normal-thresh 5"; # default 5
$dragen_parameter[] = "--vc-enable-unequal-ntd-errors false"; # disables model to correct FFPE errors.. TODO get to work with model
$dragen_parameter[] = "--vc-combine-phased-variants-distance 1"; # Merge variants if they are directly adjecent on the same strand (2 SNVs -> 1 MNP)

//SV calling:
if ($out_sv != "")
{
	$dragen_parameter[] = "--enable-sv true";
	$dragen_parameter[] = "--sv-use-overlap-pair-evidence true";
	if ($is_targeted)
	{
		$dragen_parameter[] = "--sv-exome true";
	}
}

$parser->log("DRAGEN parameters:", $dragen_parameter);

//run
$parser->exec("dragen_reset", "");
$parser->exec("LANG=en_US.UTF-8 dragen", implode(" ", $dragen_parameter)); //LANG is necessary to avoid the error "locale::facet::_S_create_c_locale name not valid" if the locale from the ssh source shell is not available on the Dragen server 
if ($debug)
{
	list($stdout) = exec2("ls $working_dir");
	$parser->log("working_dir content:", $stdout);
}

// ********************************* add filters and copy data back *********************************//

$vcf = $parser->tempFile("_unpacked.vcf");
$parser->exec("bgzip", "-d -c "$working_dir."output.vcf.gz > $vcf", true);

//left-align
$vcf_aligned = $parser->tempFile("_aligned.vcf");
$parser->exec(get_path("ngs-bits")."VcfLeftNormalize", "-stream -in $vcf -out $vcf_aligned -ref ".genome_fasta("GRCh38"), true);

//sort
$vcf_sorted = $parser->tempFile("_sorted.vcf");
$parser->exec(get_path("ngs-bits")."VcfSort","-in $vcf_aligned -out $vcf_sorted", true);


//################################################################################################
//Filter variants
//################################################################################################

$variants = Matrix::fromTSV($vcf_sorted);
$variants_filtered = new Matrix();

//get indicies
$colnames = $variants->getHeaders();
$colidx_tumor = array_search($tumor_name, $colnames);
$colidx_normal = array_search($normal_name, $colnames);

//quality cutoffs (taken from Strelka)
$min_td = 20;
$min_taf = 0.05;
$min_tsupp = 3;
$min_nd = 20;
$max_naf_rel = 1/6;

//set comments and column names
$filter_format = '#FILTER=<ID=%s,Description="%s">';
$comments = [
	sprintf($filter_format, "all-unknown", "Allele unknown"),
	sprintf($filter_format, "special-chromosome", "Special chromosome"),
	sprintf($filter_format, "depth-tum", "Sequencing depth in tumor is too low (< {$min_td})"),
	sprintf($filter_format, "freq-tum", "Allele frequency in tumor < {$min_taf}"),
	sprintf($filter_format, "depth-nor", "Sequencing depth in normal is too low (< {$min_nd})"),
	sprintf($filter_format, "freq-nor", "Allele frequency in normal > ".number_format($max_naf_rel, 2)." * allele frequency in tumor"),
	sprintf($filter_format, "lt-3-reads", "Less than {$min_tsupp} supporting tumor reads")
	];

$variants_filtered->setComments(array_merge($variants->getComments(), $comments));
$variants_filtered->setHeaders($colnames);

$count_passing = 0;
for($i = 0; $i < $variants->rows(); ++$i)
{
	$row = $variants->getRow($i);

	$ref = $row[3];
	$alt = $row[4];
	$format = $row[8];
	$tumor = $row[$colidx_tumor];
	$normal = $row[$colidx_normal];

	$filters = [];

	$filter = array_diff(explode(";", $row[6]), ["."]);

	if (!preg_match("/^[acgtACGT]*$/", $alt))
	{
		$filter[] = "all-unknown";
	}
	if (chr_check($row[0], 22, false) === FALSE)
	{
		$filter[] = "special-chromosome";
	}
	$calls = [];

	list($td, $tf) = vcf_dragen_var($format, $tumor, $alt);
	list($nd, $nf) = vcf_dragen_var($format, $normal, $alt);
	$calls[] = [ $alt, $td, $tf, $nd, $nf, $filter ];

	foreach ($calls as $call)
	{
		$variant = $row;
		list($alt, $td, $tf, $nd, $nf, $filter) = $call;
		$variant[4] = $alt;

		if ($td * $tf < $min_tsupp) $filter[] = "lt-3-reads";
		if ($td < $min_td) $filter[] = "depth-tum";
		if ($nd < $min_nd) $filter[] = "depth-nor";
		if ($tf < $min_taf) $filter[] = "freq-tum";
		if ($nf > $max_naf_rel * $tf) $filter[] = "freq-nor";

		if (empty($filter))
		{
			$filter[] = "PASS";
			$count_passing++;
		}
	}
	
	$variant[6] = implode(";", $filter);
	$variants_filtered->addRow($variant);
}
$vcf_filtered = $parser->tempFile("_filtered.vcf");
$variants_filtered->toTSV($vcf_filtered);

//remove invalid variant
$vcf_invalid = $parser->tempFile("_filtered_invalid.vcf");
$parser->exec(get_path("ngs-bits")."VcfFilter", "-remove_invalid -in $vcf_filtered -out $vcf_invalid", true);
$final = $vcf_invalid;

//flag off-target variants
if (!empty($target))
{
	$vcf_offtarget = $parser->tempFile("_filtered.vcf");
	$parser->exec(get_path("ngs-bits")."VariantFilterRegions", "-in $final -mark off-target -reg $target -out $vcf_offtarget", true);
	$final = $vcf_offtarget;
}

//remove artefacts specific for processing system (blacklist)
//artefacts are caused e.g. by hairpin sequences when using enzymatic digestion
//how artefacts are determined is documented in /mnt/storage3/users/ahsturm1/Sandbox/2023_01_31_twist_indel_artefacts/
if (!empty($target))
{
	$artefact_vcf = repository_basedir()."/data/misc/enzymatic_digestion_artefacts/".basename($target, ".bed").".vcf";
	if (file_exists($artefact_vcf))
	{
		$vcf_with_artefacts = dirname($out)."/".basename($out, ".vcf.gz")."_with_enzymatic_artefacts.vcf.gz";
		$parser->exec("bgzip", "-c $final > $vcf_with_artefacts", true);
		$vcf_no_artefacts = $parser->tempFile("_filtered_no_artefacts.vcf");
		$parser->exec(get_path("ngs-bits")."VcfSubstract", "-in $final -in2 $artefact_vcf -out $vcf_no_artefacts");
		$final = $vcf_no_artefacts;
	}
}


//zip and index output file
$parser->exec("bgzip", "-c $final > $out", true);
$parser->exec("tabix", "-p vcf $out", true);

if ($out_sv != "")
{
	$parser->log("Copying SVs output folder");
	$parser->copyFile($working_dir."output.sv.vcf.gz", $out_sv);
	$parser->copyFile($working_dir."output.sv.vcf.gz.tbi", $out_sv.".tbi");
}

// ********************************* cleanup *********************************//

//delete working directory
$parser->exec("rm", "-rf $working_dir");

if ($debug)
{
	list($stdout) = exec2("dragen_lic");
	$parser->log("dragen licence status after:", $stdout);
}

//print to STDOUT executed successfully (because there is no exit code from SGE after a job has finished)
print "DRAGEN successfully finished!";

?>
