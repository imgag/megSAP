<?php 
/** 
	@page vc_cnvkit
*/
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_cnvkit", "Copy-number calling with CNVkit.");

//bam (references)
//include region
//exclude region
//binsize
//minbinsize

//antitargets?
//ref-genome

//optional
$parser->addInfile("target_region",  "Enrichment targets BED file.", true);
$parser->addInfile("exclude_region",  "BED file containing unmappable/variable/poorly sequenced regions which should be excluded.", true);
$parser->addInt("binsize", "", true, 267);
$parser->addInt("min_binsize", "", true, 267);
$parser->addInt("min_gapsize", "", true, -1);
//center/method/ploidy for cnvkit call?
//params for plot
$parser->addString("build", "The genome build to use.", true, "GRCh38");
extract($parser->parse($argv));

//create targets

//access

//create anti-targets

//calculate target coverage

//calculate anti-target coverage

//create reference

//fix coverage

//calculate discrete copy-numbers

//create plots

//determine gender

//create pdf






$parser->addInfile("cov", "Coverage file for sample (tab-separated file with columns chr, start, end, coverage).", false);
$parser->addInfile("cov_folder", "Coverage files folder.", false);
$parser->addInfile("bed", "BED file with annotations e.g. GC-content and gene names.", false);
$parser->addOutFile("out", "Output file in TSV format.", false);
$parser->addInt("threads", "The maximum number of threads used.", true, 1);
//optional
$parser->addInt("cov_min", "Minimum number of referece coverage files required for CNV analysis.", true, 10);
$parser->addInt("cov_max", "Maximum number of referece coverage files used for CNV analysis. This parameter is needed to keep run-time and RAM requirement manageable.", true, 150);
$parser->addInt("cov_compare_max", "Maximum number of coverage files to compare during similarity calculation. Only possible with NGSD support enabled.", true, 600);
$parser->addInt("max_cnvs", "Number of expected CNVs (~200 for WES and ~2000 for WGS).", true, 2000);
$parser->addInt("max_tries", "Maximum number of tries for calling ClinCNV (R parallelization sometimes breaks with no reason", true, 10);
$parser->addInt("regions", "Number of subsequent regions that must show a signal for a call.", true, 2);
$parser->addFlag("skip_super_recall", "Skip super-recall (down to one region and log-likelihood 3).");
$parser->addFlag("mosaic","Detect additionally large mosaic regions (for WES, WGS and shallow WGS");
//tumor only flags (use for somatic with no normal sample)
$parser->addFlag("tumor_only", "Analyze tumor sample without a paired normal sample.");
$parser->addInfile("bed_off","Off-target bed file.",true); //s_dna
$parser->addInfile("cov_off","Off-target coverage file for normal sample",true);
$parser->addString("cov_folder_off", "Folder with off-target normal coverage files (if different from [data_folder]/coverage/[system_short_name]_off_target/.", true, "auto");
$parser->addString("baf_folder","Folder containing files with B-Allele frequencies.",true);
$parser->addString("cov_off_min","Minimum number of off-target coverage files required for CNV analysis.",true, 10);





?>
