<?php
/**
	@page capa_diagnostic
	@todo remove report generation
	@todo add pharmacogenomics analsis of germline
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
$steps_all = array("ma", "vc", "an", "ci", "db");
$parser->addString("steps", "Comma-separated list of processing steps to perform. Available are: ".implode(",", $steps_all), true, "ma,vc,an,ci,db");
$parser->addFlag("abra", "Turn on ABRA realignment.");
$parser->addFlag("amplicon", "Turn on amplicon mode.");
$parser->addFlag("nsc", "Skip sample correlation check (only in pair mode).");
$parser->addFlag("all_variants", "Do not use strelka filter.", true);
$parser->addFlag("strelka1","",true);
extract($parser->parse($argv));

//init
$single_sample = empty($n_id) || $n_id=="na";
foreach(explode(",", $steps) as $step)
{
	if (!in_array($step, $steps_all)) trigger_error("Unknown processing step '$step'!", E_USER_ERROR);
}

// run somatic pipeline
$s_txt = "";
$s_tsv = "";
$s_vcf = "";
$s_seg = "";
$t_bam = $p_folder."/Sample_".$t_id."/".$t_id.".bam";
$n_bam = $p_folder."/Sample_".$n_id."/".$n_id.".bam";
if($single_sample)
{
	$parser->log("Single sample mode.");
	$s_tsv = $o_folder."/".$t_id.".GSvar";
	$s_vcf = $o_folder."/".$t_id."_var.vcf.gz";
	$s_seg = $o_folder."/".$t_id."_cnvs.seg";
	$s_txt = $o_folder."/".$t_id."_report.txt";
}
else
{
	$parser->log("Paired Sample mode.");
	$s_tsv = $o_folder."/".$t_id."-".$n_id.".GSvar";
	$s_vcf = $o_folder."/".$t_id."-".$n_id."_var.vcf.gz";
	$s_seg = $o_folder."/".$t_id."-".$n_id."_cnvs.seg";
	$s_txt = $o_folder."/".$t_id."-".$n_id."_report.txt";
}

$extras = array("-filter_set non_coding_splicing,off_target,set_somatic_capa", "-steps {$steps}");
if($abra) $extras[] = "-abra";
if($amplicon) $extras[] = "-amplicon";
if($nsc) $extras[] = "-nsc";
if($all_variants) $extras[] = "-keep_all_variants_strelka";
if (isset($t_sys)) $extras[] = "-t_sys $t_sys";
$extras[] = $single_sample ? "-n_id na" : "-n_id $n_id";
if (!$single_sample && isset($n_sys)) $extras[] = "-n_sys $n_sys";
$parser->execTool("Pipelines/somatic_dna.php", "-p_folder $p_folder -t_id $t_id -o_folder $o_folder ".implode(" ", $extras));

$system_t = load_system($t_sys, $t_id);
if(empty($system_t['target_file']))	
{
	trigger_error("Tumor target file empty; no report generation.", E_USER_WARNING);
}
else
{
	// determine target region statistics
	$target = $system_t['target_file'];
	$target_name = $system_t['name_short'];
	list($stdout) = $parser->exec(get_path("ngs-bits")."BedInfo", "-in $target", false);
	$bases_overall = explode(":", $stdout[1]);
	$bases_overall = trim($bases_overall[1]);
	$regions_overall = explode(":", $stdout[0]);
	$regions_overall = trim($regions_overall[1]);

	//determine indices of important variant columns
	$var = Matrix::fromTSV($s_tsv);
	if ($var->rows()!=0)
	{
		$fi_idx = $var->getColumnIndex("filter");
		$qua_idx = $var->getColumnIndex("snp_q");
		$tvf_idx = $var->getColumnIndex("tumor_af");
		$td_idx = $var->getColumnIndex("tumor_dp");
		$co_idx = $var->getColumnIndex("coding_and_splicing");
		$cm_idx = $var->getColumnIndex("COSMIC", true);
		if(!$single_sample)
		{
			$nvf_idx = $var->getColumnIndex("normal_af");
			$nd_idx = $var->getColumnIndex("normal_dp");
		}
	}

	//get low_cov_statistics
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

	$var_report = new Matrix();
	if ($var->rows()!=0)
	{
		$headers = array('Position', 'W', 'M', 'QUALITY', 'F/T_Tumor', 'COSMIC', 'cDNA');
		if(!$single_sample)	$headers = array('Position', 'W', 'M', 'QUALITY', 'F/T_Tumor', 'F/T_Normal', 'COSMIC', 'cDNA');
		$var_report->setHeaders($headers);
		
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
			$report_row = array($row[0].":".$row[1]."-".$row[2], $row[3], $row[4], $row[$qua_idx], number_format($row[$tvf_idx],3)."/".$row[$td_idx], $cosmic, $coding);
			if(!$single_sample)	$report_row = array($row[0].":".$row[1]."-".$row[2], $row[3], $row[4], $row[$qua_idx], number_format($row[$tvf_idx],3)."/".$row[$td_idx], number_format($row[$nvf_idx],3)."/".$row[$nd_idx], $cosmic, $coding);
			$var_report->addRow($report_row);
		}
	}

	//get external sample name
	$tumor_name = basename($t_bam, ".bam");
	list($tname) = explode("_", $tumor_name);
	$tex_name = get_external_sample_name($tname, false);

	$nex_name = "n/a";
	if(!$single_sample)
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
	if(!$single_sample)	$report[] = "Normal: $normal_name (Externer Name: $nex_name)";
	$report[] = "Datum: ".date("d.m.Y");
	$report[] = "Revision der Analysepipeline: ".repository_revision(true);
	$report[] = "Analysesystem: $target_name";
	$report[] = "Zielregionen: $regions_overall ($bases_overall Basen gesamt)";

	//qc
	$report[] = "";
	$t_qcml = $p_folder."/Sample_".$t_id."/".$t_id."_stats_map.qcML";
	if(!$single_sample)	$n_qcml = $p_folder."/Sample_".$n_id."/".$n_id."_stats_map.qcML";
	$report[] = "QC:";
	$report[] = "  Coverage Tumor 100x: ".get_qc_from_qcml($t_qcml, "QC:2000030", "target region 100x percentage");
	$report[] = "  Durchschnittl. Tiefe Tumor: ".get_qc_from_qcml($t_qcml, "QC:2000025", "target region read depth");
	if(!$single_sample)	$report[] = "  Coverage Normal 100x: ".get_qc_from_qcml($n_qcml, "QC:2000030", "target region 100x percentage");
	if(!$single_sample)	$report[] = "  Durchschnittl. Tiefe Normal: ".get_qc_from_qcml($n_qcml, "QC:2000025", "target region read depth");

	// variants
	$report[] = "";
	$report[] = "Varianten:";
	$report[] = "  Gefundene Varianten: ".$var_report->rows();
	if($var_report->rows() > 0)
	{
		$report[] = "  Details:";
		$report[] = "    #".implode("\t", $var_report->getHeaders());
		for ($i=0; $i<$var_report->rows(); ++$i)
		{
			$report[] = "    ".implode("\t",$var_report->getRow($i));
		}
	}

	//CNVs
	$report[] = "";
	$report[] = "CNVs:";
	$cnv_file = $o_folder."/".$t_id.($single_sample ? "" : "-".$n_id)."_cnvs.tsv";
	if(is_file($cnv_file))
	{
		$file = file($cnv_file);
		$genes_amp = array();
		$genes_del = array();
		foreach($file as $line)
		{
			$line = trim($line);
			if ($line=="" || starts_with($line, "##")) continue;
			$report[] = "  ".$line;
			
			//generate list of amplified/deleted genes
			if (!starts_with($line, "#"))
			{
				$parts = explode("\t", $line);
				$zscore = array_sum(explode(",",$parts[7]));
				if ($zscore>0)
				{
					$genes_amp[] = $parts[9];
				}
				else
				{
					$genes_del[] = $parts[9];
				}
			}
		}
		if (count($genes_amp)>0 || count($genes_del)>0)
		{
			$report[] = "";
			$report[] = "  Amplifications: ".implode(",", array_unique($genes_amp));
			$report[] = "  Deletions: ".implode(",", array_unique($genes_del));
		}
		else
		{
			$report[] = "  Keine CNVs gefunden. Siehe $cnv_file fuer QC-Fehler.";
		}
	}
	else
	{
		$report[] = "  Keine CNV-Datei gefunden ($cnv_file).";
	}
	$report[] = "";

	// low coverage
	$report[] = "Lücken:";
	$report[] = "  Target: $target_name ($target)";
	$report[] = "  min. Tiefe Tumor: {$td}x";
	if(!$single_sample)	$report[] = "  min. Tiefe Normal: {$nd}x";
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
	$report_file = $o_folder."/".$t_id."_report.txt";
	if(!$single_sample)	$report_file = $o_folder."/".$t_id."-".$n_id."_report.txt";
	file_put_contents($report_file, implode("\n", $report));


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
	
?>
