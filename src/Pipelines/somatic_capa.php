<?php
/**
	@page capa_diagnostic
	@todo submit to queue => shell script
	@todo clean up code (merged analysis_report)
	@todo check single sample mode (multiple Error messages)
	@todo combine BedInfo and lowCoverage from tumor and normal sample
*/

$basedir = dirname($_SERVER['SCRIPT_FILENAME'])."/../";
require_once($basedir."Common/all.php");

//parse command line arguments
$parser = new ToolBase("somatic_capa", "\$Rev: 909 $", "Differential analysis of tumor/reference.");
$parser->addString("pf",  "Project folder.", false);
$parser->addString("tn",  "Tumor sample processed sample number.", false);
$parser->addString("o_folder", "Folder where output will be generated.", false);
//optional
$parser->addInfile("tn_sys",  "Tumor sample processing system INI file (determined from 'tn' by default).", true);
$parser->addString("nn",  "Reference sample GS-number.", true, "");
$parser->addInfile("nn_sys",  "Reference sample processing system INI file (determined from 'nn' by default).", true);
$parser->addInt("td",  "Min-depth for tumor low-coverage / reports.", true, 50);
$parser->addInt("nd",  "Min-depth for normal low-coverage / reports.", true, 20);
$parser->addFlag("abra", "Turn on ABRA realignment.");
$parser->addFlag("amplicon", "Turn on amplicon mode.");
$parser->addFlag("no_db", "No DB import.");
$parser->addFlag("nsc", "Skip sample correlation check (only in pair mode).");
$parser->addFlag("all_variants", "Do not use strelka filter.", true);
extract($parser->parse($argv));

//choose single sample or sample pair mode
$s_txt = "";
$s_tsv = "";
$t_bam = $pf."/Sample_".$tn."/".$tn.".bam";
$n_bam = $pf."/Sample_".$nn."/".$nn.".bam";
if(empty($nn))
{
	$parser->log("Single sample mode.");
	$s_tsv = $o_folder."/".$tn.".tsv";
	$s_txt = $o_folder."/".$tn."_report.txt";
	$extras = "";
	if($amplicon)	$extras .= "-amplicon ";
	if($no_db)	$extras .= "-steps ma,vc,an ";
	if($nsc)	$extras .= "-nsc ";
	$extras .= "-filter_set somatic_diag_capa ";
	if (isset($tn_sys)) $extras .= "-sys_tum $tn_sys ";
	$parser->execTool("php $basedir/Pipelines/somatic_dna.php", "-p_folder . -t_id $tn -n_id na -o_folder $o_folder $extras");	
}
else
{
	$parser->log("Paired Sample mode.");
	$extras = "";
	$s_tsv = $o_folder."/".$tn."-".$nn.".GSvar";
	$s_txt = $o_folder."/".$tn."-".$nn."_report.txt";
	$tmp = "";
	if($abra)	$extras .= "-abra ";
	if($amplicon)	$extras .= "-amplicon ";
	if($no_db)	$extras .= "-steps ma,vc,an,ci ";
	if($nsc)	$extras .= "-nsc ";
	if(!$all_variants)	$extras .= "-reduce_variants_strelka ";
	$extras .= "-filter_set somatic_diag_capa ";
	if (isset($tn_sys)) $extras .= "-sys_tum $tn_sys ";
	if (isset($nn_sys)) $extras .= "-sys_nor $nn_sys ";
	$parser->execTool("php $basedir/Pipelines/somatic_dna.php", "-p_folder . -t_id $tn -n_id $nn -o_folder $o_folder $extras");
}

$system_t = load_system($tn_sys, $tn);
if(empty($system_t['target_file']))	
{
	trigger_error("Tumor target file empty; no report generation.", E_USER_WARNING);
}
else
{
	$var = Matrix::fromTSV($s_tsv);
	// determine target region statistics
	$target = $system_t['target_file'];
	$target_name = $system_t['name_short'];
	list($stdout) = $parser->exec(get_path("ngs-bits")."BedInfo", "-in $target", false);
	$bases_overall = explode(":", $stdout[1]);
	$bases_overall = trim($bases_overall[1]);
	$regions_overall = explode(":", $stdout[0]);
	$regions_overall = trim($regions_overall[1]);

	//determine indices of important variant columns
	$columns = array(0, 1, 2, 3, 4, 5);
	if ($var->rows()!=0)
	{
		$fi_idx = $var->getColumnIndex("filter");
		$vf_idx = $var->getColumnIndex("tumor_af");
		$d_idx = $var->getColumnIndex("tumor_dp");
		$co_idx = $var->getColumnIndex("coding_and_splicing");
		$re_idx = $var->getColumnIndex("variant_type");
		$ge_idx = $var->getColumnIndex("gene");
		$cm_idx = $var->getColumnIndex("COSMIC", true);
		$ma_idx = $var->getColumnIndex("1000g", true);
		$sn_idx = $var->getColumnIndex("dbSNP", true);
		$al1_idx = $var->getColumnIndex("ihdb_allsys_hom");
		$al2_idx = $var->getColumnIndex("ihdb_allsys_het");
		$vu_idx = $var->getColumnIndex("classification");
	}

	//filter variants
	$var_report = new Matrix();
	$var_count_target = 0;
	$statistics_target = get_path("data_folder")."/enrichment/tsCaPa_2013_09_14.bed";

	//get low_cov_statistics
	$target_merged = $parser->tempFile("_merged.bed");
	$parser->exec(get_path("ngs-bits")."BedMerge", "-in $statistics_target -out $target_merged", false);
	//calculate low-coverage regions
	$t_low_cov = $parser->tempFile("_tlowcov.bed");
	$parser->exec(get_path("ngs-bits")."BedLowCoverage", "-in $target_merged -bam $t_bam -out $t_low_cov -cutoff ".$td, false);
	$n_low_cov = $parser->tempFile("_nlowcov.bed");
	$parser->exec(get_path("ngs-bits")."BedLowCoverage", "-in $target_merged -bam $n_bam -out $n_low_cov -cutoff ".$nd, false);
	$low_cov = $parser->tempFile("_lowcov.bed");
	$parser->exec("cat", "$t_low_cov $n_low_cov > $low_cov", false);
	$i_low_cov = $parser->tempFile("_ilowcov.bed");
	$parser->exec(get_path("ngs-bits")."BedMerge", "-in $low_cov -out $i_low_cov", false);
	//annotate gene names (with extended to annotate splicing regions as well)
	$low_cov = $o_folder."/".$tn."-".$nn."_stat_lowcov.bed";
	$parser->exec(get_path("ngs-bits")."BedAnnotateGenes", "-in $i_low_cov -extend 25 -out $low_cov", true);

	if ($var->rows()!=0)
	{
		$var_report->setHeaders(array('Position', 'W', 'M', 'QUALITY', 'COSMIC', 'AA'));
		for($i=0; $i<$var->rows(); ++$i)
		{
			$row = $var->getRow($i);

			//skip filtered variants
			if(!empty($row[$fi_idx]))	continue;

			//format: Gen c.DNA p.AA
			$tmp_row = array();
			$coding = $row[$co_idx];
			if($coding == "" || empty($coding))
			{
				continue;
			}
			else
			{
				foreach(explode(",", $row[$co_idx]) as $c)
				{
					if(empty($c))	continue;
					$eff = explode(":", $c);
					$tmp = $eff[0]." ".$eff[5]." ".$eff[6];
					if(in_array($tmp, $tmp_row)) continue;
					$tmp_row[] = $tmp;					
				}				

			}
			//generate string
			$coding = array_shift($tmp_row);
			if(count($tmp_row) > 0)
			{
				$coding .= " (identisch zu ".implode(", ", $tmp_row).")";					
			}

			//extract cosmic id
			$cosmic = $row[$cm_idx];

			//report line
			$var_report->addRow(array($row[0].":".$row[1]."-".$row[2], $row[3], $row[4], $row[$vf_idx]."/".$row[$d_idx], $cosmic, $coding));
		}				
	}

	//get external sample name
	$tumor_name = basename($t_bam, ".bam");
	list($tname) = explode("_", $tumor_name);
	$tex_name = get_external_sample_name($tname, false);

	$nex_name = "n/a";
	if(!empty($nn))
	{
		$normal_name = basename($n_bam, ".bam");
		list($nname) = explode("_", $normal_name);		
		$nex_name = get_external_sample_name($nname, false);
	}


	// report - produce output
	$report = array();
	$report[] = "Report zur bioinformatischen Analyse:";
	// meta data
	$report[] = "";
	$report[] = "Tumor: $tumor_name (Externer Name: $tex_name)";
	$report[] = "Normal: $normal_name (Externer Name: $nex_name)";
	$report[] = "Datum: ".date("d.m.Y");
	$report[] = "Revision der Analysepipeline: ".get_svn_rev();
	$report[] = "Analysesystem: $target_name";
	$report[] = "Zielregionen: $regions_overall ($bases_overall Basen gesamt)";
	//qc
	$report[] = "";
	$t_qcml = $pf."/Sample_".$tn."/".$tn."_stats_map.qcML";
	$n_qcml = $pf."/Sample_".$nn."/".$nn."_stats_map.qcML";
	$report[] = "QC:";
	$report[] = "  Coverage Tumor 50x: ".get_qc_from_qcml($t_qcml, "QC:2000029", "target region 50x percentage");
	$report[] = "  Durchschnittl. Tiefe Tumor: ".get_qc_from_qcml($t_qcml, "QC:2000025", "target region read depth");
	$report[] = "  Coverage Normal 50x: ".get_qc_from_qcml($n_qcml, "QC:2000029", "target region 50x percentage");
	$report[] = "  Durchschnittl. Tiefe Normal: ".get_qc_from_qcml($n_qcml, "QC:2000025", "target region read depth");
	// variants
	$report[] = "";
	$report[] = "Varianten:";
	$report[] = "  Gefundene Varianten: ".$var->rows();
	$report[] = "  Varianten nach Filtern: ".$var_report->rows();
	if($var_report->rows() > 0)
	{
		$report[] = "  Details:";
		$report[] = "    #".implode("\t", $var_report->getHeaders());
		for ($i=0; $i<$var_report->rows(); ++$i)
		{
			$report[] = "    ".implode("\t",$var_report->getRow($i));
		}
	}
	$report[] = "";

	// low coverage
	$report[] = "Abdeckung:";
	$report[] = "  Target: tsCaPa ($statistics_target)";
	$report[] = "  min. Tiefe Tumor: {$td}x";
	$report[] = "  min. Tiefe Normal: {$nd}x";
	$report[] = "    #gene	no. bases	chr	region(s)";
	$regions = Matrix::fromTSV($low_cov);
	$genes = array();
	$bases_lowcov = 0;
	for ($i=0; $i<$regions->rows(); ++$i)
	{
		list($chr, $start, $end, $gene) = $regions->getRow($i);
		if ($gene=="") $gene = "n/a";

		$count = $end - $start + 1;
		$bases_lowcov += $count;
		if (!isset($reg_by_gene[$gene]))
		{
			$reg_by_gene[$gene] = "$start-$end";
			$count_by_gene[$gene] = $count;
			$chr_by_gene[$gene] = $chr;
		}
		else
		{
			$reg_by_gene[$gene] .= ", $start-$end";
			$count_by_gene[$gene] += $count;
		}
	}
	ksort($reg_by_gene);
	foreach($reg_by_gene as $gene => $regions)
	{
		$bases = $count_by_gene[$gene];
		$chr = $chr_by_gene[$gene];
		$report[] = "    $gene	$bases	$chr	$regions";
	}
	$report[] = "";

	//store output
	$report_file = $o_folder."/".$tn."-".$nn."_report.txt";
	file_put_contents($report_file, implode("\n", $report));
}
