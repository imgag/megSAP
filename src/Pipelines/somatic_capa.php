<?php
/**
	@page capa_diagnostic
	@todo use best practice filter set for somatic_dna, do extra synonymous filtering in this script
	@todo report - add pharmacogenomics analysis of germline
	@todo report - change depth information based on target (exome - 20x, panel - 100x)
	@todo report - change filter settings based on target (exome vs. panel - 100x)
	@todo report - add tumor content information and HPO terms from GenLAB
	@todo report - add CN status to CNV caption
	@todo move cgi annotation to annotation step
	@todo use single letter code for amino acids
	@todo cap TumorContentEstimate at 100 %
	@todo build drug table based on additional cgi columns
	@todo variant filter region adds off-target several times, filter duplicated
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("somatic_capa", "Differential analysis of tumor/reference.");
$parser->addString("p_folder",  "Project folder.", false);
$parser->addString("t_id",  "Tumor DNA-sample processed sample ID.", false);
$parser->addString("o_folder", "Folder where output will be generated.", false);
//optional
$parser->addInfile("t_sys",  "Tumor sample processing system INI file (determined from 't_id' by default).", true);
$parser->addString("n_id",  "Reference DNA-sample processing ID.", true, "");
$parser->addInfile("n_sys",  "Reference sample processing system INI file (determined from 'n_id' by default).", true);
$parser->addInt("td",  "Min-depth for tumor low-coverage / reports.", true, 100);
$parser->addInt("nd",  "Min-depth for normal low-coverage / reports.", true, 100);
$steps_all = array("ma", "vc", "an", "ci", "db","re");
$parser->addString("steps", "Comma-separated list of processing steps to perform. Available are: ".implode(",", $steps_all), true, implode(",", $steps_all));
$parser->addFlag("abra", "Turn on ABRA realignment.");
$parser->addFlag("amplicon", "Turn on amplicon mode.");
$parser->addFlag("nsc", "Skip sample correlation check (only in pair mode).");
$parser->addFlag("all_variants", "Do not use strelka filter.", true);
$parser->addInfile("filter",  "Bed-file for output filtering. Fourth column should be gene name.", true);
$parser->addFlag("strelka1","",true);
extract($parser->parse($argv));

//init
$single_sample = empty($n_id) || $n_id=="na";
$tmp_steps = array();
foreach(explode(",", $steps) as $step)
{
	if (!in_array($step, $steps_all)) trigger_error("Unknown processing step '$step'!", E_USER_ERROR);
	$tmp_steps[] = $step;
}
$steps = $tmp_steps;

// run somatic pipeline
$s_txt = "";
$s_tsv = "";
$s_vcf = "";
$s_seg = "";
$t_bam = $p_folder."/Sample_".$t_id."/".$t_id.".bam";
$n_bam = $p_folder."/Sample_".$n_id."/".$n_id.".bam";
$n_gsvar = $p_folder."/Sample_".$n_id."/".$n_id.".GSvar";
$n_cnvs = $p_folder."/Sample_".$n_id."/".$n_id."_cnvs.tsv";
if($single_sample)
{
	$parser->log("Single sample mode.");
	$s_tsv = $o_folder."/".$t_id.".GSvar";
	$s_vcf = $o_folder."/".$t_id."_var.vcf.gz";
	$s_seg = $o_folder."/".$t_id."_cnvs.seg";
	$s_cnvs = $o_folder."/".$t_id."_cnvs.tsv";
	$s_txt = $o_folder."/".$t_id."_report.txt";
}
else
{
	$parser->log("Paired Sample mode.");
	$s_tsv = $o_folder."/".$t_id."-".$n_id.".GSvar";
	$s_vcf = $o_folder."/".$t_id."-".$n_id."_var.vcf.gz";
	$s_seg = $o_folder."/".$t_id."-".$n_id."_cnvs.seg";
	$s_cnvs = $o_folder."/".$t_id."-".$n_id."_cnvs.tsv";
	$s_txt = $o_folder."/".$t_id."-".$n_id."_report.txt";
}

$system_t = load_system($t_sys, $t_id);
$system_n = array();
if(!$single_sample)	$system_n = load_system($n_sys, $n_id);

$panel = true;
if($system_t['type'] == "WES")	$panel = false;	// exome sequencing

// run somatic data analysis
$available_steps = array("ma", "vc", "an", "ci", "db");	//steps that can be used for somatic_dna script
if(count($tmp_steps=array_intersect($available_steps,$steps))>0)
{
	// run somatic_dna pipeline
	$extras = array("-steps ".implode(",",$tmp_steps));
	$extras[] = "-filter_set not-coding-splicing,synonymous";
	if($abra) $extras[] = "-abra";
	if($amplicon) $extras[] = "-amplicon";
	if($nsc) $extras[] = "-nsc";
	if($all_variants) $extras[] = "-keep_all_variants_strelka";
	if (isset($t_sys)) $extras[] = "-t_sys $t_sys";
	if($strelka1)	$extras[] = "-strelka1";
	$extras[] = $single_sample ? "-n_id na" : "-n_id $n_id";
	if (!$single_sample && isset($n_sys)) $extras[] = "-n_sys $n_sys";
	$parser->execTool("Pipelines/somatic_dna.php", "-p_folder $p_folder -t_id $t_id -o_folder $o_folder ".implode(" ", $extras));

	// add target region statistics (needed by GSvar)
	$target = $system_t['target_file'];
	$target_name = $system_t['name_short'];
	list($stdout) = $parser->exec(get_path("ngs-bits")."BedInfo", "-in $target", false);
	$bases_overall = explode(":", $stdout[1]);
	$bases_overall = trim($bases_overall[1]);
	$regions_overall = explode(":", $stdout[0]);
	$regions_overall = trim($regions_overall[1]);	
	//get low_cov_statistics, needed for GSvar
	$target_merged = $parser->tempFile("_merged.bed");
	$parser->exec(get_path("ngs-bits")."BedMerge", "-in $target -out $target_merged", true);
	//calculate low-coverage regions
	$low_cov = $o_folder."/".$t_id.($single_sample ? "" : "-".$n_id)."_stat_lowcov.bed";
	$parser->exec(get_path("ngs-bits")."BedLowCoverage", "-in $target_merged -bam $t_bam -out $low_cov -cutoff ".$td, true);
	if(!$single_sample)
	{
		$low_cov_n = $parser->tempFile("_nlowcov.bed");
		$parser->exec(get_path("ngs-bits")."BedLowCoverage", "-in $target_merged -bam $n_bam -out $low_cov_n -cutoff ".$nd, true);
		$parser->exec(get_path("ngs-bits")."BedAdd", "-in $low_cov -in2 $low_cov_n -out $low_cov", true);
		$parser->exec(get_path("ngs-bits")."BedMerge", "-in $low_cov -out $low_cov", true);
	}
	//annotate gene names (with extended to annotate splicing regions as well)
	$parser->exec(get_path("ngs-bits")."BedAnnotateGenes", "-in $low_cov -extend 25 -out $low_cov", true);
}

// run germline data analysis
$available_steps = array("vc","an","db");	//steps that can be used for germline script, mapping already done previously
if(!$single_sample && count($tmp_steps=array_intersect($available_steps,$steps))>0)
{
	if(in_array("vc",$tmp_steps))	$tmp_steps[] = "cn";
	$parser->execTool("/Pipelines/analyze.php", "-folder $p_folder/Sample_$n_id/ -name $n_id -steps ".implode(",",$tmp_steps)." -system $n_sys --log $p_folder/Sample_$n_id//analyze_".date("Ymdhis").".log");

	// change output format of germline variant calling to format comparable to strelka e.g. with allele frequencies needed for tumor normal comparison
	// 1. annotate tumor allele frequencies, normal allele frequencies
	$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency", "-in $n_gsvar -bam $t_bam -out $n_gsvar -depth -name tumor",true);
	$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency", "-in $n_gsvar -bam $n_bam -out $n_gsvar -depth -name normal",true);
}

// TODO move this to annotation
// annotate drug association
// prepare drug association column
$path_cgi = "/mnt/share/data/dbs/CGI/2017_09_01_cgi_biomarkers_per_variant.tsv";
$cgi_database = Matrix::fromTSV($path_cgi);
$idx_assoc_genes = $cgi_database->getColumnIndex("Gene");
$idx_evide_level = $cgi_database->getColumnIndex("Evidence level");
$idx_trans = $cgi_database->getColumnIndex("transcript");

// add drug association to somatic variants
$annotate = array($s_tsv, $n_gsvar);
foreach($annotate as $a)
{
	$variants = Matrix::fromTSV($a);
	$idx_genes = $variants->getColumnIndex("gene");
	
	if(in_array("CGI_drug_assoc",$variants->getHeaders()))	
	{
		$idx_cgi1 = $variants->getColumnIndex("CGI_drug_assoc");
		$variants->removeCol($idx_cgi1);
	}
	if(in_array("CGI_evid_level",$variants->getHeaders()))	
	{
		$idx_cgi1 = $variants->getColumnIndex("CGI_evid_level");
		$variants->removeCol($idx_cgi1);
	}
	if(in_array("CGI_transcript",$variants->getHeaders()))	
	{
		$idx_cgi1 = $variants->getColumnIndex("CGI_transcript");
		$variants->removeCol($idx_cgi1);
	}

	$new_col1 = array();
	$new_col2 = array();
	$new_col3 = array();
	for($i=0;$i<$variants->rows();++$i)
	{
		// extract genes
		$genes = explode(",",$variants->getRow($i)[$idx_genes]);
		
		// compare to cgi
		$tmp1 = array();	
		$tmp2 = array();	
		$tmp3 = array();	
		foreach($genes as $gene)
		{
			for($j=0;$j<$cgi_database->rows();++$j)
			{
				$row = $cgi_database->getRow($j);
				if($gene == $row[$idx_assoc_genes])
				{
					$tmp1[] = $row[$idx_assoc_genes];
					$tmp2[] = $row[$idx_evide_level];
					$tmp3[] = $row[$idx_trans];
				}
			}
		}
		$new_col1[] = implode(", ",array_unique($tmp1));
		$new_col2[] = implode(", ",array_unique($tmp2));
		$new_col3[] = implode(", ",array_unique($tmp3));
	}
	$variants->addCol($new_col1, "CGI_drug_assoc", "Cancer Genome Interpreter - therapeutic biomarkers");
	$variants->addCol($new_col2, "CGI_evid_level", "Cancer Genome Interpreter - evidence level");
	$variants->addCol($new_col3, "CGI_transcript", "Cancer Genome Interpreter - transcript to use");
	$variants->toTSV($a);
}

// add drug association to somatic CNVs; filter pre-clinical data
$cnv_analysis = true;
if(!is_file($s_cnvs))	$cnv_analysis = false;
if($cnv_analysis)
{
	$annotate = array($s_cnvs, $n_cnvs);
	foreach($annotate as $a)
	{
		$variants = Matrix::fromTSV($a);
		$idx_genes = $variants->getColumnIndex("genes");
		if(in_array("CGI_drug_assoc",$variants->getHeaders()))	
		{
			$idx_cgi1 = $variants->getColumnIndex("CGI_drug_assoc");
			$variants->removeCol($idx_cgi1);
		}
		if(in_array("CGI_evid_level",$variants->getHeaders()))	
		{
			$idx_cgi1 = $variants->getColumnIndex("CGI_evid_level");
			$variants->removeCol($idx_cgi1);
		}
			
		$new_col1 = array();
		$new_col2 = array();
		for($i=0;$i<$variants->rows();++$i)
		{		
			// extract genes and compare to CGI
			$tmp1 = array();	
			$tmp2 = array();	
			$genes = explode(",",$variants->getRow($i)[$idx_genes]);
			foreach($genes as $gene)
			{
				for($j=0;$j<$cgi_database->rows();++$j)
				{
					$row = $cgi_database->getRow($j);
					// TODO reconsider evidence level filter
					if($gene == $row[$idx_assoc_genes] && $row[$idx_evide_level]!="Pre-clinical")
					{
						$tmp1[] = $row[$idx_assoc_genes];
						$tmp2[] = $row[$idx_evide_level];
					}
				}
			}
			$new_col1[] = implode(", ",array_unique($tmp1));
			$new_col2[] = implode(", ",array_unique($tmp2));
		}
		$variants->addCol($new_col1, "CGI_drug_assoc", "Cancer Genome Interpreter - therapeutic biomarkers");
		$variants->addCol($new_col2, "CGI_evid_level", "Cancer Genome Interpreter - evidence level");
		$variants->toTSV($a);
	}

	// filter variants if specific bed is available
	$s_cnvsfil = $s_cnvs;
	$n_cnvsfil = $n_cnvs;
	if(!empty($filter))
	{
		$s_cnvsfil = dirname($s_cnvs)."/".basename($s_cnvs,".tsv")."_filtered.tsv";
		$n_cnvsfil = dirname($n_cnvs)."/".basename($n_cnvs,".tsv")."_filtered.tsv";
		
		// SNVs - soft filter
		$parser->exec(get_path("ngs-bits")."VariantFilterRegions","-in $s_tsv -mark -reg $filter -out $s_tsv", true);
		$parser->exec(get_path("ngs-bits")."VariantFilterRegions","-in $n_gsvar -mark -reg $filter -out $n_gsvar", true);
		// CNVs - hard filter
		$matrix = Matrix::fromTSV($filter);
		foreach(array($s_cnvs, $n_cnvs) as $file)
		{
			$cnvs = Matrix::fromTSV($file);
			$idx_genes = $cnvs->getColumnIndex("genes");

			$tmp = new Matrix();
			$tmp->setHeaders($cnvs->getHeaders());
			$tmp->setComments($cnvs->getComments());
			for($i=0;$i<$cnvs->rows();++$i)
			{
				$cnv = $cnvs->getRow($i);
				
				$found = false;
				$tmp_genes = array();
				for($j=0;$j<$matrix->rows();++$j)
				{
					$chr = $matrix->getRow($j)[0];
					$start = $matrix->getRow($j)[1]+1;
					$end = $matrix->getRow($j)[2];
					$gene = $matrix->getRow($j)[3];
					
					// filter range
					if($chr!=$cnv[0])	continue;
					if(range_overlap($cnv[1],$cnv[2],$start,$end))	$found = true;
					
					// filter gene names
					if(in_array($gene,explode(",",$cnv[$idx_genes])))	$tmp_genes[] = $gene;
				}
				if($found)
				{
					$cnv[$idx_genes] = implode(",",array_unique($tmp_genes));
					$tmp->addRow($cnv);
				}
			}
			if($file==$s_cnvs)	$tmp->toTSV($s_cnvsfil);
			if($file==$n_cnvs)	$tmp->toTSV($n_cnvsfil);
		}
	}
}

if (in_array("re", $steps))
{
	if(empty($system_t['target_file']))	
	{
		trigger_error("Tumor target file empty; no report generation.", E_USER_WARNING);
	}
	else
	{
		$tumor_name = basename($t_bam, ".bam");
		list($tname) = explode("_", $tumor_name);
		$tex_name = get_external_sample_name($tname, false);
		$nex_name = "n/a";
		$t_qcml = $p_folder."/Sample_".$t_id."/".$t_id."_stats_map.qcML";
		$n_qcml = $p_folder."/Sample_".$n_id."/".$n_id."_stats_map.qcML";
		$s_qcml = $o_folder."/".$t_id."-".$n_id."_stats_som.qcML";
		if(!$single_sample)
		{
			$normal_name = basename($n_bam, ".bam");
			list($nname) = explode("_", $normal_name);		
			$nex_name = get_external_sample_name($nname, false);
		}

		$filter = array("BRCA1", "BRCA2", "TP53", "STK11", "PTEN", "MSH2", "MSH6", "MLH1", "PMS2", "APC", "MUTYH", "SMAD4", "VHL", "MEN1", "RET", "RB1", "TSC1", "TSC2", "NF2", "WT1","SDHB","SDHD","SDHC","SDHAF2","BMPR1A");
		sort($filter);
		$hpo_terms = array();
		$tumor_content = "n/a %";
		if(isset(get_path("db_host")['GL8']))
		{
			list($tn,) = explode("_",$t_id);
			$db = DB::getInstance("GL8");
			$res = $db->executeQuery("SELECT * FROM v_ngs_sap WHERE labornummer LIKE '$tn' ");
			if(!empty($res))
			{
				if(!empty($res[0]['TUMORANTEIL']))	$tumor_content = $res[0]['TUMORANTEIL'];
				if(!empty($res[0]['HPOTERM1']))	$hpo_terms[] = $res[0]['HPOTERM1'];
				if(!empty($res[0]['HPOTERM2']))	$hpo_terms[] = $res[0]['HPOTERM2'];
				if(!empty($res[0]['HPOTERM3']))	$hpo_terms[] = $res[0]['HPOTERM3'];
				if(!empty($res[0]['HPOTERM4']))	$hpo_terms[] = $res[0]['HPOTERM4'];
			}
		}
		$obo = Obo::fromObo("/mnt/share/data/dbs/HPO/hp.obo");
		foreach($hpo_terms as $k => $hpt)
		{
			if(empty($hpt))	continue;
			$hpo_terms[$k] = $obo->getTermNameByID($hpt);
		}
		$hpo = implode(", ", array_filter($hpo_terms));
		if(empty($hpo))	$hpo = "keine Angabe";
		
		// 1. prepare variant data for somatic report
		$snv_somatic_report = report_SNV_table($s_tsv, $single_sample, "tumor_af", "tumor_dp", "normal_af", "normal_dp", 0);
		if(!$single_sample)
		{
			$snv_germline_report = report_SNV_table($n_gsvar, $single_sample, "tumor_freq", "tumor_depth", "normal_freq", "normal_depth", 4, $filter);
		}
		
		if(isset($s_cnvsfil) && is_file($s_cnvsfil))
		{
			$min_zscore = 5;
			$min_regions = 10;
			$keep = array("MYC","MDM2","MDM4","CDKN2A","CDKN2A-AS1","CDK4","CDK6","PTEN","CCND1","RB1","CCND3","BRAF","KRAS","NRAS");
			sort($keep);
			$cnv_somatic_report = report_CNV_table($s_cnvsfil, $single_sample, $min_regions, $min_zscore, $keep);
			if(!$single_sample)
			{
				$min_zscore = 4;
				$min_regions = 3;
				$cnv_germline_report = report_CNV_table($n_cnvsfil, $single_sample, $min_regions, $min_zscore, array(), $filter);
			}
			
			$genes_amp = array();
			$genes_amp = array();
			$genes_loss = array();
			for($i=0;$i<$cnv_somatic_report->rows();++$i)
			{
				$row = $cnv_somatic_report->getRow($i);
				if($row[3]>2)	$genes_amp = array_merge($genes_amp, explode(",",$row[4]));
				if($row[3]<2)	$genes_loss = array_merge($genes_loss, explode(",",$row[4]));
			}
			//
			$comment = $cnv_somatic_report->getComments();
			if(!empty($comment) && !empty($comment[0]))
			{
				$cnvs = explode("; ", $comment[0]);
				foreach($cnvs as $cnv)
				{
					list($g,$cn) = explode(" ",$cnv, 2);
					$cn = trim($cn,'()');
					$cns = explode(", ",$cn);
					$cn = trim($cns[0],'CN ');
					if(!is_numeric($cn[0]))	trigger_error("'$cn' in '$cnv' is not numeric.", E_USER_ERROR);
					if($cn[0]>2)	$genes_amp[] = $g;
					if($cn[0]<2)	$genes_loss[] = $g;
				}
			}
		}

		// 2. prepare metadata and qc data for somatic report		
		$metadata = array();
		$metadata[] = array("Datum:", date("d.m.Y"));
		$metadata[] = array("Tumor:", $tumor_name);
		if(!$single_sample)	$metadata[] = array("Normal:", $normal_name);
		$metadata[] = array("Tumoranteil histol./molekular:","\highlight1 ca. $tumor_content / ".(!is_numeric(get_qc_from_qcml($s_qcml, "QC:2000054", "tumor content estimate"))?get_qc_from_qcml($s_qcml, "QC:2000054", "tumor content estimate"):" ca. ".number_format(get_qc_from_qcml($s_qcml, "QC:2000054", "tumor content estimate"), 0))." %\highlight0");
		$metadata[] = array("Diagnose:","\highlight1 $hpo\highlight0");
		
		$qc = array();
		$qc[] = array("Revision der Analysepipeline:", repository_revision(true));
		if($panel)	$qc[] = array("Coverage Tumor 100x:",get_qc_from_qcml($t_qcml, "QC:2000030", "target region 100x percentage")." %");
		else	$qc[] = array("Coverage Tumor 30x:",get_qc_from_qcml($t_qcml, "QC:2000028", "target region 30x percentage")." %");
		$qc[] = array("Durchschnittl. Tiefe Tumor:",get_qc_from_qcml($t_qcml, "QC:2000025", "target region read depth")."x");
		if(!$single_sample)	
		{
			if($panel)	$qc[] = array("Coverage Normal 100x:",get_qc_from_qcml($n_qcml, "QC:2000030", "target region 100x percentage")." %");
			else	$qc[] = array("Coverage Normal 30x:",get_qc_from_qcml($n_qcml, "QC:2000028", "target region 30x percentage")." %");
			$qc[] = array("Durchschnittl. Tiefe Normal:", get_qc_from_qcml($n_qcml, "QC:2000025", "target region read depth")."x");
		}
		$qc[] = array("Prozessierungssystem Tumor:", $system_t['name_manufacturer']);
		if(!$single_sample)	$qc[] = array("Prozessierungssystem Normal:", $system_n['name_manufacturer']);

		// 3. generate somatic rtf
		$report = array();
		$report[] = "{\rtf1\ansi\deff0 {\fonttbl {\f0 Times New Roman;}}{\colortbl;\red188\green230\blue138;}";
		$report[] = "{\header\pard\qr\fs14 Wiss. Bericht $tumor_name-$normal_name vom ".date("d.m.Y").", Seite \chpgn  von {\field{\*\fldinst  NUMPAGES }}\par}";
		$report[] = "\titlepg";
		$report[] = "{\headerf {\pard\sa2200 \par}}";
		$report[] = "{\footerf {\pard\sa1200 \par}}";
		// document head
		$report[] = "{\pard\sa270\fs28\b Wissenschaftlicher Bericht zu Tumorsequenzierung\par}";
		// metadata
		$report[] = "{\pard\sa45\fs24\b Allgemeine Informationen:\par}";
		$report[] = "{\pard\fs20".rtf_table_rows($metadata,array(3000,6500))."\par}";
		// variants
		$report[] = "{\pard\fs24\sa45\b Varianten\par}";
		$report[] = rtf_table_row(array("Die ausführlich gelisteten Varianten sind Biomarker, die in der CGI-Datenbank (www.cancergenomeinterpreter.org) gelistet werden. Die verbleibenden Gene mit Varianten werden im Block 'weitere Varianten' bzw. 'weitere Gene' gelistet. Details zu diesen Varianten sind auf Nachfrage erhältlich. Eine funktionelle Einschätzung der Varianten finden Sie im beigefügten QCI-Bericht. Bitte beachten Sie, dass bei geringem Tumoranteil die Ergebnisse aufgrund niedriger Allelfrequenzen ungenauer sein können. Veränderungen, die in der jeweils zugehörigen Blutprobe gefunden wurden, können ein Hinweis auf eine erbliche Erkrankung darstellen. In diesem Fall empfehlen wir eine Vorstellung in einer humangenetischen Beratungsstelle. Weitere Informationen zur Untersuchungstechnik, -umfang und Auswertung finden sich in der Rubrik Qualitätsparameter."), array(9500),array("\qj"));
		$report[] = "{\pard\fs20\sb180\sa45\b SNVs und kleine INDELs:\par}";
		$report[] = rtf_table_row(array("Mutationslast: ".get_qc_from_qcml($s_qcml, "QC:2000053", "somatic variant rate")), array(9500),array("\qj"));
		$report[] = rtf_table_row(array("Position","R","V","F/T Tumor","F/T Normal","cDNA"), array(2550,3200,3850,5050,6250,9500),array("\qc","\qc", "\qc","\qc","\qc","\qc"),true,true);
		$report[] = rtf_table_row(array("Somatische Veränderungen mit therapeutischer Relevanz"), array(9500),array("\qc"),true,true);
		if($snv_somatic_report->rows()>0)
		{
			for ($i=0; $i<$snv_somatic_report->rows(); ++$i)
			{
				$report[] = rtf_table_row($snv_somatic_report->getRow($i),array(2550,3200,3850,5050,6250,9500),array("","\qc", "\qc","\qc","\qc","\qj"),true);
			}
		}
		else	$report[] = rtf_table_row(array("keine"),array(9500),array("\qj"),true);
		$add_vars = implode(", ", $snv_somatic_report->getComments());
		if(!empty($add_vars))	$report[] = rtf_table_row(array("weitere Varianten:", $add_vars),array(2550,9500),array("\qr","\qj"),true);
		$report[] = rtf_table_row(array("Keimbahnvarianten*"), array(9500),array("\qc"),true,true);
		if($snv_germline_report->rows()>0)
		{
			for ($i=0; $i<$snv_germline_report->rows(); ++$i)
			{
				$report[] = rtf_table_row($snv_germline_report->getRow($i),array(2550,3200,3850,5050,6250,9500),array("","\qc", "\qc","\qc","\qc","\qj"),true);
			}
		}
		else	$report[] = rtf_table_row(array("keine"),array(9500),array("\qj"),true);
		$add_vars = implode(", ", $snv_germline_report->getComments());
		if(!empty($add_vars))	$report[] = rtf_table_row(array("weitere Varianten:", $add_vars),array(2550,9500),array("\qr","\qj"),true);
		$report[] = "{\pard\sa45\sb20";
		$report[] = "\fs14 Position - chromosomale Position (".$system_t['build']."); R - Allel Referenz; V - Allel Variante; F/T Tumor - Frequenz und Tiefe im Tumor; F/T Normal - Frequenz und Tiefe im Normal; cDNA - cDNA Position und Auswirkung Peptid.\fs20";
		$report[] = "\par}";
		// CNVs
		if($cnv_analysis)
		{
			$report[] = "{\pard\fs20\sb180\sa45\b CNVs:\par}";
			$report[] = rtf_table_row(array("Die folgenden Tabellen zeigen das wissenschaftliche Ergebnis der CNV-Analysen, die mit dem CNVHunter Tool durchgeführt wurden. Zur Validierung relevanter Veränderungen empfehlen wir eine zweite, unabhängige Methode. Die geschätzte Copy Number ist abhängig von Tumorgehalt und Verteilung der CNV-tragenden Zellen im Tumorgewebe."), array(9500), array("\qj"));
			$report[] = rtf_table_row(array("Position","Größe","Typ","CN","Gene"), array(2550,3800,4800,5300,9500),array("\qc","\qc","\qc","\qc","\qc"),true,true);
			$report[] = rtf_table_row(array("Somatische Veränderungen mit therapeutischer Relevanz"), array(9500),array("\qc"),true,true);
			if($cnv_somatic_report->rows()>0)
			{
				for ($i=0; $i<$cnv_somatic_report->rows(); ++$i)
				{
					$report[] = rtf_table_row($cnv_somatic_report->getRow($i), array(2550,3800,4800,5300,9500),array("","\qc","\qc","\qc","\qj"),true);
				}
			}
			else	$report[] = rtf_table_row(array("keine"),array(9500),array("\qj"),true);
			$add_vars = implode(", ", $cnv_somatic_report->getComments());
			if(!empty($add_vars))	$report[] = rtf_table_row(array("weitere Gene:", $add_vars),array(2550,9500),array("\qr","\qj"),true);
			$report[] = rtf_table_row(array("Keimbahn*"), array(9500),array("\qc"),true,true);
			if($cnv_germline_report->rows()>0)
			{
				for ($i=0; $i<$cnv_germline_report->rows(); ++$i)
				{
					$report[] = rtf_table_row($cnv_germline_report->getRow($i), array(2550,3800,4800,5300,9500),array("","\qc","\qc","\qc","\qj"),true);
				}
			}
			else	$report[] = rtf_table_row(array("keine"),array(9500),array("\qj"),true);
			$add_vars = implode(", ",$cnv_germline_report->getComments());
			if(!empty($add_vars))	$report[] = rtf_table_row(array("weitere Varianten:", "\highlight1".$add_vars."\highlight0"),array(2550,9500),array("\qr","\qj"),true);
			$report[] = "{\pard\sa180\sb20\fs14 Position - chromosomale Position (".$system_t['build']."); Größe - CNV-Größe in Basenpaaren; Typ - Verlust (LOSS) oder Amplifikation (AMP); CN - geschätzte mediane Copy Number in dieser chrom. Region;Gene - Gene in dieser Region.\par}";
			// qci
			$report[] = rtf_table_row(array("\highlight1Amplifzierte Gene:\highlight0","\highlight1 ".(count($genes_amp)>0 ? implode(", ",array_unique($genes_amp)) : "keine")."\highlight0"), array(1800,9000),array("","\qj"));
			$report[] = rtf_table_row(array("\highlight1Deletierte Gene:\highlight0","\highlight1 ".(count($genes_loss)>0 ? implode(", ",array_unique($genes_loss)) : "keine")."\highlight0"), array(1800,9000),array("","\qj"));
		}
		// qc
		$report[] = "{\pard\fs24\sa45\sb180\b Qualitätsparameter:\par}";
		$report[] = rtf_table_rows($qc,array(3500,9500));
		$report[] = rtf_table_row(array(""),array(9500));
		$report[] = rtf_table_row(array("* Folgende Gene wurden in die Auswertung der Blutprobe einbezogen: ".implode(", ",$filter)."."),array(9500),array("\qj"));
		$report[] = "}";
		
		//store output
		$report_file = $o_folder."/".$t_id."_report.rtf";
		if(!$single_sample)	$report_file = $o_folder."/".$t_id."-".$n_id."_report.rtf";
		file_put_contents($report_file, str_replace(array("\r","\t","\v","\e","\f"), array("\\r","\\t","\\v","\\e","\\f"), implode("", $report)));

		//  prepare IGV-session for all files
		if(is_dir($p_folder) && is_dir($o_folder))
		{
			$rel_path = relative_path($o_folder, $p_folder);
			$igv_session = array();
			$igv_session[] = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n";
			$igv_session[] = "<Session genome=\"1kg_v37\" hasGeneTrack=\"true\" hasSequenceTrack=\"true\" locus=\"all\" path=\".\" version=\"8\">\n";
			$igv_session[] = "\t<Resources>\n";
			if(is_file($s_seg))   $igv_session[] = "\t\t<Resource path=\"".$rel_path."/".$o_folder."/".basename($s_seg)."\"/>\n";
			if(is_file($s_tsv))   $igv_session[] = "\t\t<Resource path=\"".$rel_path."/".$o_folder."/".basename($s_vcf)."\"/>\n";
			if(is_file($t_bam))	$igv_session[] = "\t\t<Resource path=\"".$rel_path."/Sample_".basename($t_bam,".bam")."/".basename($t_bam)."\"/>\n";
			if(is_file($n_bam))	$igv_session[] = "\t\t<Resource path=\"".$rel_path."/Sample_".basename($n_bam,".bam")."/".basename($n_bam)."\"/>\n";
			$igv_session[] = "\t</Resources>\n";
			$igv_session[] = "\t<HiddenAttributes>\n";
			$igv_session[] = "\t\t<Attribute name=\"NAME\"/>\n";
			$igv_session[] = "\t\t<Attribute name=\"DATA FILE\"/>\n";
			$igv_session[] = "\t\t<Attribute name=\"DATA TYPE\"/>\n";
			$igv_session[] = "\t</HiddenAttributes>\n";
			$igv_session[] = "</Session>";
			$fn = $t_id;
			if(!$single_sample)	$fn = $t_id."-".$n_id;
			file_put_contents($o_folder."/".$fn."_igv.xml",$igv_session);
		}
		else
		{
			trigger_error("IGV-Session-File was not created. Folder $p_folder or $o_folder does not exist.",E_USER_WARNING);
		}
	}
}

function report_SNV_table($gsvar, $single_sample, $col_tum_freq, $col_tum_depth, $col_nor_freq, $col_nor_depth, $min_class = 0, $filter = array())
{
	$snv = Matrix::fromTSV($gsvar);
	
	$germline = false;
	if($min_class > 0)	$germline = true;

	$classification = false;
	if(in_array("classification",$snv->getHeaders()))	$classification = true;;

	$snv_report = new Matrix();
	if ($snv->rows()!=0)
	{
		$idx_chr = $snv->getColumnIndex("chr");
		$idx_start = $snv->getColumnIndex("start");
		$idx_end = $snv->getColumnIndex("end");
		$idx_ref = $snv->getColumnIndex("ref");
		$idx_obs = $snv->getColumnIndex("obs");
		$idx_ge = $snv->getColumnIndex("gene");
		$idx_cgi1 = $snv->getColumnIndex("CGI_drug_assoc");
		$idx_cgi2 = $snv->getColumnIndex("CGI_transcript");
		$idx_class = null;
		if($germline && $classification)	$idx_class = $snv->getColumnIndex("classification");		
		$idx_fi = $snv->getColumnIndex("filter");
		$idx_tvf = $snv->getColumnIndex($col_tum_freq);
		$idx_td = $snv->getColumnIndex($col_tum_depth);
		$idx_co = $snv->getColumnIndex("coding_and_splicing");
		if(!$single_sample)
		{
			$nvf_idx = $snv->getColumnIndex($col_nor_freq);
			$nd_idx = $snv->getColumnIndex($col_nor_depth);
		}
		$idx_gene = $snv->getColumnIndex("gene");

		$headers = array('Position', 'W', 'M', 'F/T Tumor', 'cDNA');
		if(!$single_sample)	$headers = array('Position', 'W', 'M', 'F/T Tumor', 'F/T Normal', 'cDNA');
		$snv_report->setHeaders($headers);
		
		//extract relevant MISO terms for filtering
		$obo = Obo::fromOBO(repository_basedir()."/data/dbs/Ontologies/so-xp_3_0_0.obo");
		$miso_terms_coding = array();
		$ids = array("SO:0001580","SO:0001568");
		foreach($ids as $id)
		{
			$miso_terms_coding[] = $obo->getTermNameByID($id);
			$tmp = $obo->getTermChildrenNames($id,true);
			$miso_terms_coding = array_merge($miso_terms_coding,$tmp);
		}
		$miso_terms_coding = array_unique($miso_terms_coding);
		
		$caption_genes = array();
		$caption_cod = array();
		for($i=0; $i<$snv->rows(); ++$i)
		{
			$row = $snv->getRow($i);
			
			// skip filtered variants
			if(!$germline && !empty($row[$idx_fi]))	continue;
			if($germline && $classification && $row[$idx_class]<$min_class)	continue;

			if(!empty($filter))
			{
				$skip = true;
				foreach($filter as $f)
				{
					$tmp = explode(",",$row[$idx_ge]);
					if(in_array($f,$tmp))	$skip = false;
				}
				if($skip)	continue;
			}
			
			// only report drug assoc genes in detail
			$skip = true;
			$cgi_drug_assoc = $row[$idx_cgi1];
			if(empty($cgi_drug_assoc))
			{
				// use first coding information to add p or c position
				list($cod_spli,) = explode(",",$row[$idx_co]);
				list($g,,,,,$c,$p) = explode(":",$cod_spli);		
				if(empty($g))	continue;

				$caption_genes[] = $g;
				continue;
			}

			// convert format: Gen c.DNA p.AA and use CGI transcript
			$tmp_row = array();
			$coding = $row[$idx_co];
			list($cgi_transcript,) = explode(", ", $row[$idx_cgi2]);
			if($coding=="" || empty($coding))
			{
				continue;
			}
			else
			{
				foreach(explode(",", $row[$idx_co]) as $c)
				{
					if(empty($c))	continue;
					
					// filter non coding / splicing codings
					$skip = true;
					$eff = explode(":", $c);
					foreach(explode("&",$eff[2]) as $t)
					{
						if(in_array($t,$miso_terms_coding) || in_array($t."_variant",$miso_terms_coding))	$skip = false;
					}
					if($skip)	continue;
					
					// only keep one transcript
					if(!empty($cgi_transcript) && $eff[1]!=$cgi_transcript)	continue;
					if(!empty($tmp_row))	continue;
					
					// filter duplicates
					$tmp = $eff[0]." ".$eff[5]." ".$eff[6];
					if(in_array($tmp, $tmp_row)) continue;

					$tmp_row[] = $tmp." (".(empty($cgi_transcript)?$eff[1]:$cgi_transcript).")";					
				}				

			}
			//generate string
			$coding = array_shift($tmp_row);
			if(count($tmp_row) > 0)
			{
				$coding .= " (identisch zu ".implode(", ", $tmp_row).")";					
			}

			//report line
			$report_row = array($row[$idx_chr].":".$row[$idx_start]."-".$row[$idx_end], $row[$idx_ref], $row[$idx_obs], number_format($row[$idx_tvf],3)."/".$row[$idx_td], $coding);
			if(!$single_sample)	$report_row = array($row[$idx_chr].":".$row[$idx_start]."-".$row[$idx_end], $row[$idx_ref], $row[$idx_obs], number_format($row[$idx_tvf],3)."/".$row[$idx_td], number_format($row[$nvf_idx],3)."/".$row[$nd_idx], $coding);

			$snv_report->addRow($report_row);
		}
		
		sort($caption_genes);
		$caption = "";
		foreach($caption_genes as $cp)
		{
			$caption .= ", ".$cp;
		}
		$snv_report->addComment(trim($caption,', '));
	}	
	return $snv_report;
}

function report_CNV_table($tsv, $single_sample, $min_regions, $min_zscore, $keep = array(), $filter = array())
{
	$germline = false;
	if(!empty($filter))	$germline = true;

	$cnvs = Matrix::fromTSV($tsv);
	
	$cnv_report = new Matrix();
	$caption_genes = array();
	$caption_cn = array();
	if ($cnvs->rows()!=0)
	{
		$idx_chr = $cnvs->getColumnIndex("chr");
		$idx_start = $cnvs->getColumnIndex("start");
		$idx_end = $cnvs->getColumnIndex("end");
		$idx_size = $cnvs->getColumnIndex("size");
		$idx_zscores = $cnvs->getColumnIndex("region_zscores");
		$idx_copy = $cnvs->getColumnIndex("region_copy_numbers");
		$idx_genes = $cnvs->getColumnIndex("genes");
		$idx_cgi1 = $cnvs->getColumnIndex("CGI_drug_assoc");

		$headers = array('Position','Größe','Typ','Gene');
		$cnv_report->setHeaders($headers);
		
		for($i=0;$i<$cnvs->rows();++$i)
		{
			$line = $cnvs->getRow($i);
			
			//generate list of amplified/deleted genes
			$zscores = explode(",",$line[$idx_zscores]);
				
			// filter CNVs according to best practice filter criteria
			$skip = false;
			$tmp_zscores = array();
			foreach($zscores as $z)
			{
				$tmp_zscores[] = abs($z);
			}
			if(max($tmp_zscores)<$min_zscore)	$skip = true;
			if(count($zscores)<$min_regions)	$skip = true;
			foreach($keep as $k)
			{
				$tmp = explode(",",$line[$idx_genes]);
				if(in_array($k,$tmp) && (max($tmp_zscores)>=$min_zscore))	$skip = false;
			}
			if($skip)	continue;
			
			// filter CNVs for gene array
			if(!empty($filter))
			{
				$skip = true;
				$tmp_genes = array();
				foreach($filter as $f)
				{
					$tmp = explode(",",$line[$idx_genes]);
					if(in_array($f,$tmp))
					{
						$skip = false;
						$tmp_genes[]= $f;
					}
				}
				$line[$idx_genes] = implode(",",$tmp_genes);
			}
			if($skip)	continue;
			
			$copy_median = median(explode(",",$line[$idx_copy]));

			$cgi_drug_assoc = explode(",",$line[$idx_genes]);
			if(!$germline)
			{
				$genes = explode(",", $line[$idx_genes]);
				$cgi_drug_assoc = explode(", ", $line[$idx_cgi1]);
				foreach($genes as $g)
				{
					if(!in_array($g,$cgi_drug_assoc))
					{
						$caption_genes[] = $g;
						if(!isset($caption_cn[$g]))	$caption_cn[$g] = array();
						$caption_cn[$g][] = $copy_median;
					}
				}

				if(empty($line[$idx_cgi1]))	continue;
			}

			$zscore_sum = array_sum(explode(",",$line[$idx_zscores]));
			if ($zscore_sum>0)
			{
				$cnv_report->addRow(array($line[$idx_chr].":".$line[$idx_start]."-".$line[$idx_end],$line[$idx_size],"AMP",$copy_median,implode(", ", $cgi_drug_assoc)));
			}
			else
			{
				$cnv_report->addRow(array($line[$idx_chr].":".$line[$idx_start]."-".$line[$idx_end],$line[$idx_size],"LOSS",$copy_median,implode(", ", $cgi_drug_assoc)));
			}
		}
		
		sort($caption_genes);
		$caption = "";
		foreach($caption_genes as $cp)
		{
			$caption .= "; ".$cp." (CN ".implode(", CN", $caption_cn[$cp]).")";
		}
		$cnv_report->addComment(trim($caption,'; '));	
	}
	
	return $cnv_report;
}
?>
