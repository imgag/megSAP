<?php
/**
	@page capa_diagnostic
	@todo remove report generation
	@todo add pharmacogenomics analsis of germline
	@todo add CNVs of all variants
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
if($strelka1)	$extras[] = "-strelka1";
$extras[] = $single_sample ? "-n_id na" : "-n_id $n_id";
if (!$single_sample && isset($n_sys)) $extras[] = "-n_sys $n_sys";
$parser->execTool("Pipelines/somatic_dna.php", "-p_folder $p_folder -t_id $t_id -o_folder $o_folder ".implode(" ", $extras));


$system_t = load_system($t_sys, $t_id);
$system_n = array();
if(!$single_sample)	$system_n = load_system($n_sys, $n_id);
if(empty($system_t['target_file']))	
{
	trigger_error("Tumor target file empty; no report generation.", E_USER_WARNING);
}
else
{
	//prepare snvs for report
	$snv = Matrix::fromTSV($s_tsv);
	if ($snv->rows()!=0)
	{
		$fi_idx = $snv->getColumnIndex("filter");
		$tvf_idx = $snv->getColumnIndex("tumor_af");
		$td_idx = $snv->getColumnIndex("tumor_dp");
		$co_idx = $snv->getColumnIndex("coding_and_splicing");
		if(!$single_sample)
		{
			$nvf_idx = $snv->getColumnIndex("normal_af");
			$nd_idx = $snv->getColumnIndex("normal_dp");
		}
	}

	$snv_report = new Matrix();
	if ($snv->rows()!=0)
	{
		$headers = array('Position', 'W', 'M', 'F/T_Tumor', 'cDNA');
		if(!$single_sample)	$headers = array('Position', 'W', 'M', 'F/T_Tumor', 'F/T_Normal', 'cDNA');
		$snv_report->setHeaders($headers);
		
		//extract relevant MISO terms for filtering
		$obo = MISO::fromOBO();
		$miso_terms_coding = array();
		$ids = array("SO:0001580","SO:0001568");
		foreach($ids as $id)
		{
			$miso_terms_coding[] = $obo->getTermNameByID($id);
			$tmp = $obo->getTermChildrenNames($id,true);
			$miso_terms_coding = array_merge($miso_terms_coding,$tmp);
		}
		$miso_terms_coding = array_unique($miso_terms_coding);
		
		for($i=0; $i<$snv->rows(); ++$i)
		{
			$row = $snv->getRow($i);

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
					
					// filter non coding / splicing codings
					$skip = true;
					foreach(explode("&",$eff[2]) as $t)
					{
						if(in_array($t,$miso_terms_coding))	$skip = false;
					}
					if($skip)	continue;
					
					// filter duplicates
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

			//report line
			$report_row = array($row[0].":".$row[1]."-".$row[2], $row[3], $row[4], number_format($row[$tvf_idx],3)."/".$row[$td_idx], $coding);
			if(!$single_sample)	$report_row = array($row[0].":".$row[1]."-".$row[2], $row[3], $row[4], number_format($row[$tvf_idx],3)."/".$row[$td_idx], number_format($row[$nvf_idx],3)."/".$row[$nd_idx], $coding);
			$snv_report->addRow($report_row);
		}
	}

	//prepare cnvs for report	
	$cnv_file = $o_folder."/".$t_id.($single_sample ? "" : "-".$n_id)."_cnvs.tsv";
	$cnvs = Matrix::fromTSV($cnv_file);
	$min_zscore = 5;
	$min_regions = 10;
	$keep = array("MYC","MDM2","MDM4","CDKN2A","CDKN2A-AS1","CDK4","CDK6","PTEN","CCND1","RB1","CCND3","BRAF","KRAS","NRAS");
	
	$cnv_report = new Matrix();
	$genes_amp = array();
	$genes_loss = array();
	if ($cnvs->rows()!=0)
	{
		$headers = array('Position','Größe [bp]','Typ','Gene');
		$cnv_report->setHeaders($headers);
		
		for($i=0;$i<$cnvs->rows();++$i)
		{
			$line = $cnvs->getRow($i);
			$idx_zscores = $cnvs->getColumnIndex("region_zscores");
			$idx_genes = $cnvs->getColumnIndex("genes");
			
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
				if(in_array($k,$tmp) && abs(max($tmp_zscores))>$min_zscore)	$skip = false;
			}
			if($skip)	continue;

			$zscore_sum = array_sum(explode(",",$line[7]));
			if ($zscore_sum>0)
			{
				$genes_amp[] = $line[9];
				$cnv_report->addRow(array($line[0].":".$line[1]."-".$line[2],$line[4],"AMP",$line[$idx_genes]));
			}
			else
			{
				$genes_loss[] = $line[9];
				$cnv_report->addRow(array($line[0].":".$line[1]."-".$line[2],$line[4],"LOSS",$line[$idx_genes]));
			}
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
	$report[] = "Wissenschaftlicher Bericht zu den Proben $t_id ".( !$single_sample ? "/ $n_id" : "");
	// metadata
	$report[] = "";
	$report[] = "Datum: ".date("d.m.Y");
	$report[] = "Revision der Analysepipeline: ".repository_revision(true);
	$tmp = "Tumor: $tumor_name";
	if(!$single_sample)	$tmp .= ", Normal: $normal_name";
	$report[] = $tmp;
	$report[] = "Prozessierungssystem Tumor: ".$system_t['name_manufacturer'];
	$report[] = "Prozessierungssystem Normal: ".$system_n['name_manufacturer'];

	//qc
	$report[] = "";
	$t_qcml = $p_folder."/Sample_".$t_id."/".$t_id."_stats_map.qcML";
	$n_qcml = $p_folder."/Sample_".$n_id."/".$n_id."_stats_map.qcML";
	$s_qcml = $p_folder."/Somatic_".$t_id."-".$n_id."/".$t_id."-".$n_id."_stats_som.qcML";
	$report[] = "QC:";
	$report[] = "  Coverage Tumor 100x: ".get_qc_from_qcml($t_qcml, "QC:2000030", "target region 100x percentage");
	$report[] = "  Durchschnittl. Tiefe Tumor: ".get_qc_from_qcml($t_qcml, "QC:2000025", "target region read depth");
	if(!$single_sample)	$report[] = "  Coverage Normal 100x: ".get_qc_from_qcml($n_qcml, "QC:2000030", "target region 100x percentage");
	if(!$single_sample)	$report[] = "  Durchschnittl. Tiefe Normal: ".get_qc_from_qcml($n_qcml, "QC:2000025", "target region read depth");
	if(!$single_sample && is_file($s_qcml))	$report[] = "  Mutationslast: ".get_qc_from_qcml($s_qcml, "QC:2000053", "somatic variant rate");
	$report[] = "";

	// variants
	$report[] = "Varianten:";
	$report[] = "  Gefundene Varianten: ".$snv_report->rows();
	if($snv_report->rows() > 0)
	{
		$report[] = "  Details:";
		$report[] = "    #".implode("\t", $snv_report->getHeaders());
		for ($i=0; $i<$snv_report->rows(); ++$i)
		{
			$report[] = "    ".implode("\t",$snv_report->getRow($i));
		}
	}
	$report[] = "";

	//CNVs
	$report[] = "CNVs:";
	if(is_file($cnv_file))
	{
		$report[] = "  Gefundene CNVs: ".$cnv_report->rows();
		
		if($cnv_report->rows()>0)
		{
			$report[] = "  Details:";
			$report[] = "    #".implode("\t",$cnv_report->getHeaders());
			for ($i=0; $i<$cnv_report->rows(); ++$i)
			{
				$report[] = "    ".implode("\t",$cnv_report->getRow($i));
			}
			$report[] = "";
			$report[] = "    Amplifizierte Gene: ".implode(",", array_unique($genes_amp));
			$report[] = "    Deletierte Gene: ".implode(",", array_unique($genes_loss));
		}
		else
		{
			$report[] = "  Keine CNVs gefunden. Siehe $cnv_file für QC-Fehler.";
		}
	}
	else
	{
		$report[] = "  Keine CNV-Datei gefunden ($cnv_file).";
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
