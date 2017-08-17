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
$steps_all = array("ma", "vc", "an", "ci", "db","re");
$parser->addString("steps", "Comma-separated list of processing steps to perform. Available are: ".implode(",", $steps_all), true, "ma,vc,an,ci,db");
$parser->addFlag("abra", "Turn on ABRA realignment.");
$parser->addFlag("amplicon", "Turn on amplicon mode.");
$parser->addFlag("nsc", "Skip sample correlation check (only in pair mode).");
$parser->addFlag("all_variants", "Do not use strelka filter.", true);
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

// run somatic_dna
$available_steps = array("ma", "vc", "an", "ci", "db");	//steps that can be used for somatic_dna script
if(count($tmp_steps=array_intersect($available_steps,$steps))>0)
{
	$extras = array("-filter_set non_coding_splicing,off_target,set_somatic_capa", "-steps ".implode(",",$tmp_steps));
	if($abra) $extras[] = "-abra";
	if($amplicon) $extras[] = "-amplicon";
	if($nsc) $extras[] = "-nsc";
	if($all_variants) $extras[] = "-keep_all_variants_strelka";
	if (isset($t_sys)) $extras[] = "-t_sys $t_sys";
	if($strelka1)	$extras[] = "-strelka1";
	$extras[] = $single_sample ? "-n_id na" : "-n_id $n_id";
	if (!$single_sample && isset($n_sys)) $extras[] = "-n_sys $n_sys";
	$parser->execTool("Pipelines/somatic_dna.php", "-p_folder $p_folder -t_id $t_id -o_folder $o_folder ".implode(" ", $extras));
}

if (in_array("re", $steps))
{
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

		$snv_report = new Matrix();
		if ($snv->rows()!=0)
		{
			$idx_chr = $snv->getColumnIndex("chr");
			$idx_start = $snv->getColumnIndex("start");
			$idx_end = $snv->getColumnIndex("end");
			$idx_ref = $snv->getColumnIndex("ref");
			$idx_obs = $snv->getColumnIndex("obs");
			
			$idx_fi = $snv->getColumnIndex("filter");
			$idx_tvf = $snv->getColumnIndex("tumor_af");
			$idx_td = $snv->getColumnIndex("tumor_dp");
			$idx_co = $snv->getColumnIndex("coding_and_splicing");
			if(!$single_sample)
			{
				$nvf_idx = $snv->getColumnIndex("normal_af");
				$nd_idx = $snv->getColumnIndex("normal_dp");
			}

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
				if(!empty($row[$idx_fi]))	continue;

				//format: Gen c.DNA p.AA
				$tmp_row = array();
				$coding = $row[$idx_co];
				if($coding == "" || empty($coding))
				{
					continue;
				}
				else
				{
					foreach(explode(",", $row[$idx_co]) as $c)
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
				$report_row = array($row[$idx_chr].":".$row[$idx_start]."-".$row[$idx_end], $row[$idx_ref], $row[$idx_obs], number_format($row[$idx_tvf],3)."/".$row[$idx_td], $coding);
				if(!$single_sample)	$report_row = array($row[$idx_chr].":".$row[$idx_start]."-".$row[$idx_end], $row[$idx_ref], $row[$idx_obs], number_format($row[$idx_tvf],3)."/".$row[$idx_td], number_format($row[$nvf_idx],3)."/".$row[$nd_idx], $coding);
				$snv_report->addRow($report_row);
			}
		}

		//prepare cnvs for report
		$cnvs = Matrix::fromTSV($s_cnvs);
		$min_zscore = 5;
		$min_regions = 10;
		$keep = array("MYC","MDM2","MDM4","CDKN2A","CDKN2A-AS1","CDK4","CDK6","PTEN","CCND1","RB1","CCND3","BRAF","KRAS","NRAS");
		
		$cnv_report = new Matrix();
		$genes_amp = array();
		$genes_loss = array();
		if ($cnvs->rows()!=0)
		{
			$idx_chr = $cnvs->getColumnIndex("chr");
			$idx_start = $cnvs->getColumnIndex("start");
			$idx_end = $cnvs->getColumnIndex("end");
			$idx_size = $cnvs->getColumnIndex("size");
			$idx_zscores = $cnvs->getColumnIndex("region_zscores");
			$idx_copy = $cnvs->getColumnIndex("region_copy_numbers");
			$idx_genes = $cnvs->getColumnIndex("genes");

			$headers = array('Position','Größe [bp]','Typ','Gene');
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
					if(in_array($k,$tmp) && abs(max($tmp_zscores))>$min_zscore)	$skip = false;
				}
				if($skip)	continue;

				$copy_median = median(explode(",",$line[$idx_copy]));

				$zscore_sum = array_sum(explode(",",$line[$idx_zscores]));
				if ($zscore_sum>0)
				{
					$genes_amp = array_merge($genes_amp,explode(",",$line[$idx_genes]));
					$cnv_report->addRow(array($line[$idx_chr].":".$line[$idx_start]."-".$line[$idx_end],$line[$idx_size],"AMP",$copy_median,implode(", ", explode(",",$line[$idx_genes]))));
				}
				else
				{
					$genes_loss = array_merge($genes_loss, explode(",",$line[$idx_genes]));
					$cnv_report->addRow(array($line[$idx_chr].":".$line[$idx_start]."-".$line[$idx_end],$line[$idx_size],"LOSS",$copy_median,implode(", ", explode(",",$line[$idx_genes]))));
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


		// generate report preformatted using rtf
		$report = array();
		$report[] = "{\rtf1\ansi\deff0 {\fonttbl {\f0 Times New Roman;}}";
		$report[] = "{\footer\pard\qr\fs14 Probe $tumor_name-$normal_name, Seite \chpgn  von {\field{\*\fldinst  NUMPAGES }}\par}";
		$report[] = doc_head("Wissenschaftlicher Bericht zu den Proben $t_id ".( !$single_sample ? "/ $n_id" : ""));
		// metadata
		$report[] = par_head("Metadata:");
		$report[] = "{\pard \fs20 \sa180";
		$report[] = par_table_metadata("Datum:", date("d.m.Y"));
		$report[] = par_table_metadata("Revision der Analysepipeline:", repository_revision(true));
		$report[] = par_table_metadata("Tumor:", $tumor_name);
		if(!$single_sample)	$report[] = par_table_metadata("Normal:", $normal_name);
		$report[] = par_table_metadata("Prozessierungssystem Tumor:", $system_t['name_manufacturer']);
		if(!$single_sample)	$report[] = par_table_metadata("Prozessierungssystem Normal:", $system_n['name_manufacturer']);
		$report[] = "\par}";

		//qc
		$t_qcml = $p_folder."/Sample_".$t_id."/".$t_id."_stats_map.qcML";
		$n_qcml = $p_folder."/Sample_".$n_id."/".$n_id."_stats_map.qcML";
		$s_qcml = $p_folder."/Somatic_".$t_id."-".$n_id."/".$t_id."-".$n_id."_stats_som.qcML";
		$report[] = par_head("Quality parameters:");
		$report[] = "{\pard \fs20 \sa180";
		$report[] = par_table_qc("Coverage Tumor 100x:",get_qc_from_qcml($t_qcml, "QC:2000030", "target region 100x percentage"));
		$report[] = par_table_qc("Durchschnittl. Tiefe Tumor:",get_qc_from_qcml($t_qcml, "QC:2000025", "target region read depth"));
		if(!$single_sample)	$report[] = par_table_qc("Coverage Normal 100x:",get_qc_from_qcml($n_qcml, "QC:2000030", "target region 100x percentage"));
		if(!$single_sample)	$report[] = par_table_qc("Durchschnittl. Tiefe Normal:", get_qc_from_qcml($n_qcml, "QC:2000025", "target region read depth"));
		if(!$single_sample && is_file($s_qcml))	$report[] = par_table_qc("Mutationslast:",get_qc_from_qcml($s_qcml, "QC:2000053", "somatic variant rate"));
		$report[] = "\par}";

		// variants
		$report[] = par_head("SVNs und kleine INDELs:");
		$report[] = "{\pard \fs20 ";
		$report[] = "Gefundene Varianten: ".$snv_report->rows()."\line";
		if($snv_report->rows() > 0)
		{
			$report[] = "Details:\fs8 \line ";
			$report[] = "\par}";
			$report[] = "{\pard \fs20 \sa360";
			
			$report[] = par_table_header_var();
			for ($i=0; $i<$snv_report->rows(); ++$i)
			{
				list($pos,$ref,$obs,$freq_t,$freq_n,$coding) = $snv_report->getRow($i);
				$report[] = par_table_row_var($pos,$ref,$obs,$freq_t,$freq_n,$coding);
			}
		}
		$tmp = "\fs8 \line \fs14 Position - chromosomale Position (".$system_t['build']."); W - Allel Wildtyp; V - Allel Variante; F/T_Tumor - Frequenz und Tiefe im Tumor;";
		$tmp .= "F/T_Normal - Frequenz und Tiefe im Normal; cDNA - cDNA Position und Auswirkung Peptid.\fs20";
		$report[] = $tmp;
		$report[] = "\par}";

		//CNVs
		$report[] = par_head("CNVs:");
		if(is_file($s_cnvs))
		{			
			$report[] = "\fs8 \line \fs14";
			$report [] = "{\pard\fs20\qj";
			$report[] = "Die folgenden Tabellen zeigen das wissenschaftliche Ergebnis der CNV-Analysen, die mit dem CNVHunter Tool durchgeführt wurden. Die Liste der CNVs basiert auf strengen Filterkriterien. Zur Validierung relevanter Veränderungen empfehlen wir eine zweite, unabhängige Methode.";
			$report[] = "\par}";
			$report[] = "\fs8 \line \fs14";
			$report[] = "{\pard \fs20 ";
			$report[] = "Gefundene CNVs: ".$cnv_report->rows()."\line";
			
			if($cnv_report->rows()>0)
			{
				$report[] = "Details:\fs8 \line ";
				$report[] = "\par}";

				$report [] = "{\pard \fs20 \sa360";

				$report[] = par_table_header_cnv();
				for ($i=0; $i<$cnv_report->rows(); ++$i)
				{
					list($pos,$size,$type,$copy_number,$gene) = $cnv_report->getRow($i);
					$report[] = par_table_row_cnv($pos, $size, $type, $copy_number, $gene);
				}
				$tmp = "\fs8 \line \fs14 Position - chromosomale Position (".$system_t['build']."); Größe - CNV-Größe in Basenpaaren; Typ - Verlust (LOSS) oder Amplifikation (AMP);";
				$tmp .= "CN - geschätzte mediane Copy Number in dieser Region;Gene - Gene in dieser Region.\fs20";
				$report[] = $tmp;
				$report[] = "\par}";
				
				$report[] = "{\pard \fs20 ";
				$report[] = par_cnv_genes(array_unique($genes_amp),array_unique($genes_loss));
			}
			else
			{
				$report[] = "Keine CNVs gefunden. Siehe $s_cnvs für QC-Fehler.";
			}
		}
		else
		{
			$report[] = "Keine CNV-Datei gefunden ($s_cnvs).";
		}
		$report[] = "\par}";
		$report[] = "}";

		//store output
		$report_file = $o_folder."/".$t_id."_report.rtf";
		if(!$single_sample)	$report_file = $o_folder."/".$t_id."-".$n_id."_report.rtf";
		file_put_contents($report_file, str_replace(array("\r","\t","\v","\e","\f"), array("\\r","\\t","\\v","\\e","\\f"), implode("\n", $report)));


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

function doc_head($string)
{
	return "{\pard \fs28 \sa180 $string\par}";
}

function par_head($string)
{
	return "{\pard \fs24 $string\par}";
}

function par_table_metadata($name, $value)
{
	return "{\trowd \trgraph180 \fs20 \cellx3000 \cellx6500 \pard\intbl ".$name."\cell \pard\intbl\ql ".$value."\cell\row}";
}

function par_table_qc($name, $value)
{
	return "{\trowd \trgraph180 \fs20 \cellx2500 \cellx4200 \pard\intbl ".$name."\cell \pard\intbl\qr ".$value."\cell\row}";
}

function par_table_header_var()
{
	$string = "{\trowd \trgraph180 \fs20";
	$string .= "\clbrdrt\brdrw18\brdrs \clbrdrl\brdrw18\brdrs \clbrdrb\brdrw18\brdrs \clbrdrr\brdrw18\brdrs \cellx2450";
	$string .= "\clbrdrt\brdrw18\brdrs \clbrdrl\brdrw18\brdrs \clbrdrb\brdrw18\brdrs \clbrdrr\brdrw18\brdrs \cellx2850";
	$string .= "\clbrdrt\brdrw18\brdrs \clbrdrl\brdrw18\brdrs \clbrdrb\brdrw18\brdrs \clbrdrr\brdrw18\brdrs \cellx3250";
	$string .= "\clbrdrt\brdrw18\brdrs \clbrdrl\brdrw18\brdrs \clbrdrb\brdrw18\brdrs \clbrdrr\brdrw18\brdrs \cellx4450";
	$string .= "\clbrdrt\brdrw18\brdrs \clbrdrl\brdrw18\brdrs \clbrdrb\brdrw18\brdrs \clbrdrr\brdrw18\brdrs \cellx5650";
	$string .= "\clbrdrt\brdrw18\brdrs \clbrdrl\brdrw18\brdrs \clbrdrb\brdrw18\brdrs \clbrdrr\brdrw18\brdrs \cellx9500";
	$string .= "\pard\intbl\b\qc Position\cell \pard\intbl\qc W\cell \pard\intbl\qc V\cell \pard\intbl\qc F/T_Tumor\cell \pard\intbl\qc F/T_Normal\cell \pard\intbl\qc cDNA\cell\row}";
	return $string;
}

function par_table_row_var($pos,$ref,$obs,$freq_t,$freq_n,$coding)
{
	$string = "{\trowd \trgraph180 \fs20";
	$string .= "\clbrdrt\brdrw18\brdrs \clbrdrl\brdrw18\brdrs \clbrdrb\brdrw18\brdrs \clbrdrr\brdrw18\brdrs \cellx2450";
	$string .= "\clbrdrt\brdrw18\brdrs \clbrdrl\brdrw18\brdrs \clbrdrb\brdrw18\brdrs \clbrdrr\brdrw18\brdrs \cellx2850";
	$string .= "\clbrdrt\brdrw18\brdrs \clbrdrl\brdrw18\brdrs \clbrdrb\brdrw18\brdrs \clbrdrr\brdrw18\brdrs \cellx3250";
	$string .= "\clbrdrt\brdrw18\brdrs \clbrdrl\brdrw18\brdrs \clbrdrb\brdrw18\brdrs \clbrdrr\brdrw18\brdrs \cellx4450";
	$string .= "\clbrdrt\brdrw18\brdrs \clbrdrl\brdrw18\brdrs \clbrdrb\brdrw18\brdrs \clbrdrr\brdrw18\brdrs \cellx5650";
	$string .= "\clbrdrt\brdrw18\brdrs \clbrdrl\brdrw18\brdrs \clbrdrb\brdrw18\brdrs \clbrdrr\brdrw18\brdrs \cellx9500";
	$string .= "\pard\intbl $pos\cell \pard\intbl\qc $ref\cell \pard\intbl\qc $obs\cell \pard\intbl\qc $freq_t\cell \pard\intbl\qc $freq_n\cell \pard\intbl\qj $coding\cell\row}";
	return $string;
}


function par_table_header_cnv()
{
	$string = "{\trowd \trgraph180 \fs20";
	$string .= "\clbrdrt\brdrw18\brdrs \clbrdrl\brdrw18\brdrs \clbrdrb\brdrw18\brdrs \clbrdrr\brdrw18\brdrs \cellx2450";
	$string .= "\clbrdrt\brdrw18\brdrs \clbrdrl\brdrw18\brdrs \clbrdrb\brdrw18\brdrs \clbrdrr\brdrw18\brdrs \cellx3700";
	$string .= "\clbrdrt\brdrw18\brdrs \clbrdrl\brdrw18\brdrs \clbrdrb\brdrw18\brdrs \clbrdrr\brdrw18\brdrs \cellx4700";
	$string .= "\clbrdrt\brdrw18\brdrs \clbrdrl\brdrw18\brdrs \clbrdrb\brdrw18\brdrs \clbrdrr\brdrw18\brdrs \cellx5200";
	$string .= "\clbrdrt\brdrw18\brdrs \clbrdrl\brdrw18\brdrs \clbrdrb\brdrw18\brdrs \clbrdrr\brdrw18\brdrs \cellx9500";
	$string .= "\pard\intbl\b\qc Position\cell";
	$string .= "\pard\intbl\qc Größe [bp]\cell";
	$string .= "\pard\intbl\qc Typ\cell";
	$string .= "\pard\intbl\qc CN\cell";
	$string .= "\pard\intbl\qc Gene\cell";
	$string .= "\row}";
	return $string;
}

function par_table_row_cnv($pos,$size,$type,$copy_number,$genes)
{
	$string = "{\trowd \trgraph180 \fs20";
	$string .= "\clbrdrt\brdrw18\brdrs \clbrdrl\brdrw18\brdrs \clbrdrb\brdrw18\brdrs \clbrdrr\brdrw18\brdrs \cellx2450";
	$string .= "\clbrdrt\brdrw18\brdrs \clbrdrl\brdrw18\brdrs \clbrdrb\brdrw18\brdrs \clbrdrr\brdrw18\brdrs \cellx3700";
	$string .= "\clbrdrt\brdrw18\brdrs \clbrdrl\brdrw18\brdrs \clbrdrb\brdrw18\brdrs \clbrdrr\brdrw18\brdrs \cellx4700";
	$string .= "\clbrdrt\brdrw18\brdrs \clbrdrl\brdrw18\brdrs \clbrdrb\brdrw18\brdrs \clbrdrr\brdrw18\brdrs \cellx5200";
	$string .= "\clbrdrt\brdrw18\brdrs \clbrdrl\brdrw18\brdrs \clbrdrb\brdrw18\brdrs \clbrdrr\brdrw18\brdrs \cellx9500";
	$string .= "\pard\intbl $pos\cell";
	$string .= "\pard\intbl\qc $size\cell";
	$string .= "\pard\intbl\qc $type\cell";
	$string .= "\pard\intbl\qc $copy_number\cell";
	$string .= "\pard\intbl\qj $genes\cell";
	$string .= "\row}";
	return $string;
}

function par_cnv_genes($array_amplified, $array_loss)
{
	$string = "{\trowd \trgraph180 \fs20";
	$string .= "\cellx1800 \cellx9000";
	$string .= "\pard\intbl Amplifizierte Gene:\cell"; 
	$string .= "\pard\intbl\qj ".implode(", ",$array_amplified)."\cell";
	$string .= "\row}";
	$string .= "{\trowd \trgraph180 \fs20";
	$string .= "\cellx1800 \cellx9000";
	$string .= "\pard\intbl Deletierte Gene:\cell"; 
	$string .= "\pard\intbl\qj ".implode(", ",$array_loss)."\cell";
	$string .= "\row}";	
	return $string;
}

?>
