<?php

/**
	@page somatic_ivac
	
	@todo consider haplotyping for somatic variants and germline variants close to each other
	@todo run DNA and RNA in parallel
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("somatic_ivac", "Differential analysis of tumor/reference exomes. Can combine DNA and RNA data and results are suitable for current vaccination projects.");
$parser->addString("p_folder", "Folder that contains Sample folders.", false);
$parser->addString("t_dna_id",  "Tumor sample DNA processing ID.", false);
$parser->addString("n_dna_id",  "Reference sample DNA processing ID.", false);
$parser->addString("o_folder", "Output folder.", false);
// optional
$parser->addString("t_rna_id",  "Tumor sample RNA processing ID.", true, "na");
$parser->addString("t_rna_fo",  "Tumor sample RNA data folder.", true, "na");
$parser->addString("n_rna_id",  "Normal sample RNA processing ID.", true, "na");
$parser->addString("n_rna_fo",  "Normal sample RNA data folder.", true, "na");
$parser->addInfile("t_dna_sys",  "Tumor DNA processing system INI file (determined from 't_dna_id' by default).", true);
$parser->addInfile("n_dna_sys",  "Normal DNA processing system INI file (determined from 'n_dna_id' by default).", true);
$parser->addInfile("t_rna_sys",  "Tumor RNA processing system INI file (determined from 't_rna_id' by default).", true);
$parser->addInfile("n_rna_sys",  "Normal RNA processing system INI file (determined from 'n_rna_id' by default).", true);
$parser->addFlag("amplicon", "Amplicon library mode (no random distribution of breakpoints).", true);
$parser->addFlag("no_softclip", "Skip soft-clipping of overlapping reads. NOTE: This may increase rate of false positive variants (DNA only).", true);
$steps_all = array("ma", "vc", "fu", "an", "db", "co","im");
$parser->addString("steps", "Comma-separated list of processing steps to perform.", true, implode(",", $steps_all));
$parser->addFlag("rna_only", "Process only RNA.", true);
$parser->addFlag("dna_only", "Process only DNA.", true);
$parser->addFlag("smn_dna", "Skip mapping of DNA normal sample (was done previously).", true);
$parser->addFlag("smt_dna", "Skip mapping of DNA tumor sample (was done previously).", true);
$parser->addFlag("freebayes", "Use freebayes for variant detection (default: strelka).", true);
$parser->addFloat("contamination", "Indicates fraction of tumor cells in normal sample. Freebayes is used for variant calling in contaminated normal samples.", true, 0);
$parser->addFlag("nsc", "Skip sample correlation check.");
extract($parser->parse($argv));

// determine steps to perform
$steps = explode(",", $steps);
foreach($steps as $step)
{
	if (!in_array($step, $steps_all)) trigger_error("Unknown processing step '$step'!", E_USER_ERROR);
}
$o_folder = rtrim($o_folder,'/')."/";

$t_dna_bam = "$p_folder/Sample_$t_dna_id/$t_dna_id.bam";
$n_dna_bam = "$p_folder/Sample_$n_dna_id/$n_dna_id.bam";
$t_rna_bam = "$p_folder/Sample_$t_rna_id/$t_rna_id.bam";
$n_rna_bam = "$p_folder/Sample_$n_rna_id/$n_rna_id.bam";

// get processing systems
$t_dna_sys_file =$t_dna_sys;
$t_dna_sys = load_system($t_dna_sys_file, $t_dna_id);
$n_dna_sys_file = $n_dna_sys;
$n_dna_sys = load_system($n_dna_sys_file, $n_dna_id);
if($t_dna_sys["name_short"] != $n_dna_sys["name_short"]) trigger_error ("System tumor '".$t_dna_sys["name_short"]."' and system reference '".$n_dna_sys["name_short"]."' are different!", E_USER_WARNING);
if($t_dna_sys["build"] != $n_dna_sys["build"]) trigger_error ("System build tumor '".$t_dna_sys["build"]."' and system reference '".$n_dna_sys["build"]."' are different!", E_USER_ERROR);

// (1) run somatic_dna, always pair mode
$o_folder_dna = $o_folder;
$annotated = $o_folder.$t_dna_id."-".$n_dna_id.".GSvar";
$som_vcf = $o_folder.$t_dna_id."-".$n_dna_id."_var_annotated.vcf.gz";
$som_seg = $o_folder.$t_dna_id."-".$n_dna_id."_cnvs.seg";
$tmp_steps = array();
$available_steps = array("ma", "vc", "an", "db");	// available steps that can be used to run somatic_dna script
if(count($tmp_steps=array_intersect($available_steps,$steps))>0 && !$rna_only)
{
	$args = "-t_sys $t_dna_sys_file -n_sys $n_dna_sys_file ";
	$args .= "-steps ".implode(",",$tmp_steps)." -min_af 0.05 ";
	$args .= "-filter_set not-coding-splicing,synonymous ";
	$args .= "--log ".$o_folder_dna."somatic_ivac_dna_".date('YmdHis',mktime()).".log ";
	if($amplicon)	$args .= "-amplicon ";
	if($no_softclip)	$args .= "-no_softclip ";
	if($freebayes)	$args .= "-freebayes ";
	if($smn_dna)	$args .= "-smn ";
	if($smt_dna)	$args .= "-smt ";
	if($nsc)	$args .= "-nsc ";
	if($contamination > 0)	$args .= "-contamination $contamination ";
	if(!is_dir($o_folder_dna))	mkdir($o_folder_dna, 0775, true);
	$parser->execTool("Pipelines/somatic_dna.php", "-p_folder $p_folder -t_id $t_dna_id -n_id $n_dna_id -o_folder $o_folder_dna $args");
}

// (2) run somatic_rna, consider single sample mode
$rna = 0;	// 0 = no RNA, 1 = tumor sample RNA, 2 = tumor and normal sample RNA
if($t_rna_id!="na" && $t_rna_fo!="na" && ($n_rna_id=="na" || $n_rna_fo=="na"))	$rna = 1;
if($t_rna_id!="na" && $t_rna_fo!="na" && $n_rna_id!="na" && $n_rna_fo!="na")	$rna = 2;
$o_folder_rna = $o_folder;

$tmp_steps = array();
$available_steps = array("ma","fu","db");	//steps that can be used for somatic_rna script
if(count($tmp_steps=array_intersect($available_steps,$steps))>0 && !$dna_only)
{
	if($rna>0)
	{
		//somatic_rna pipeline starts with map and ends with db
		if(!is_dir($o_folder_rna))	mkdir($o_folder_rna, 0775, true);

		$args = "";
		$t_rna_sys_file = $t_rna_sys;
		$t_rna_sys = load_system($t_rna_sys_file, $t_rna_id);
		$args .= "-t_sys $t_rna_sys_file ";
		if($rna>1)
		{
			$n_rna_sys_file = $n_rna_sys;
			$n_rna_sys = load_system($n_rna_sys_file, $n_rna_id);
			$args .= "-n_sys $n_rna_sys_file ";
		}
		$args .= "-p_folder $p_folder -t_id $t_rna_id -n_id $n_rna_id -o_folder $o_folder_rna ";
		$args .= "-steps ".implode(",", $tmp_steps)." -t_folder $t_rna_fo -n_folder $n_rna_fo ";
		$args .= "--log ".$o_folder_rna."somatic_ivac_rna_".date('YmdHis',mktime()).".log ";
		$parser->execTool("Pipelines/somatic_rna.php", $args);
	}
}
// compare RNA tumor and DNA tumor (tumor and normal is compared by the pipelines itself
if($rna>0 && !$nsc)
{
	$output = $parser->exec(get_path("ngs-bits")."SampleCorrelation", "-in1 $t_dna_bam -in2 $t_rna_bam -bam -max_snps 4000", true);
	$parts = explode(":", $output[0][1]);
	if ($parts[1]<0.8)	trigger_error("The genotype correlation of samples $t_bam and $n_bam is ".$parts[1]."; it should be above 0.8!", E_USER_ERROR);
}

// (3) combine results of DNA / RNA
if(in_array("co", $steps) && $rna!=0)
{
	$tmp = $parser->tempFile(".vcf");
	$vaf_options = " -depth";
	$tum_rna_bam = $p_folder."/Sample_".$t_rna_id."/".$t_rna_id.".bam";
	
	$parser->exec("bgzip", "-dc $som_vcf > $tmp", false);	// no output logging, because Toolbase::extractVersion() does not return
	$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency", "-in $tmp -bam $tum_rna_bam -out $tmp -name rna_tum $vaf_options", true);
	$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency", "-in $annotated -bam $tum_rna_bam -out $annotated -name rna_tum $vaf_options", true);
	if ($rna==2)
	{
		$ref_rna_bam = $p_folder."/Sample_".$n_rna_id."/".$n_rna_id.".bam";
		$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency", "-in $tmp -bam $ref_rna_bam -out $tmp -name rna_ref $vaf_options", true);
		$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency", "-in $annotated -bam $ref_rna_bam -out $annotated -name rna_ref $vaf_options", true);
	}
	
	$parser->exec("bgzip", "-cf $tmp > $som_vcf", false);	// no output logging, because Toolbase::extractVersion() does not return
	$parser->exec("tabix", "-fp vcf $som_vcf", false);	// no output logging, because Toolbase::extractVersion() does not return
}

// (4) identify if germline variants are close to somatic variants
if(in_array("im", $steps))
{
	$transcript_ids = array();
	$transcript_type = "";
	$file_som = Matrix::fromTSV($annotated);
	$idx_cs = $file_som->getColumnIndex("coding_and_splicing");
	$idx_filter = $file_som->getColumnIndex("filter");
	$known_gene = get_path("data_folder")."/dbs/UCSC/knownGene.txt";
	$no_transcript_id = array();

	// (4a) extract transcripts and their exons for given variants
	for($i=0;$i<$file_som->rows();++$i)	// get all somatic refseq-ids
	{
		$row_som = $file_som->getRow($i);
//		if(!empty($filter))	continue;

		$tmp_ids = array();
		//REFSEQ transcript IDs
		//preg_match('/NM_\d+/',$row_som[$idx_cs],$tmp_ids);
		//ENSEMBL transcript IDs: ENST00000376030
		//preg_match('/NM_\d+/',$row_som[$idx_cs],$tmp_ids);
		foreach(explode(",",$row_som[$idx_cs]) as $t)
		{
			$tmp_id = explode(":",$t)[1];
			if(strpos($tmp_id,"ENST")==0)	$transcript_type = "Ensembl";
			else if(strpos($tmp_id,"NM_")==0)	$transcript_type = "RefSeq";
			else if(!empty($tmp_id))	trigger_error("Unkown transcript information $tmp_id.",E_USER_ERROR);
			$tmp_ids[] = $tmp_id;
		}
		
		if(count($tmp_ids) == 0)
		{
			$no_transcript_id[] = $i;
			continue;
		}
		
		foreach($tmp_ids as $tmp_id)
		{
			if(!isset($transcript_ids[$tmp_id]))
			{
				$transcript_ids[$tmp_id] = array();
				$transcript_ids[$tmp_id]['variant_rows'] = array();
			}
			$transcript_ids[$tmp_id]['variant_rows'][] = $i;
		}
	}
	if(!empty($no_transcript_id))	trigger_error ("No transcript ID was identified in rows: ".implode(", ", $no_transcript_id).". Variants were skipped.",E_USER_NOTICE);
	
	//convert transcript-ids to UCSC-ids
	$transcript_file = "";
	if($transcript_type=="Ensembl")	$transcript_file = get_path("data_folder")."/dbs/UCSC/knownToEnsembl.tsv";
	if($transcript_type=="RefSeq")	$transcript_file = get_path("data_folder")."/dbs/UCSC/knownToRefSeq.tsv";
	
	//$kgxref_file = get_path("data_folder")."/dbs/UCSC/kgXref.txt";
	$ucsc2transcript = array();
	
	$handle = fopen($transcript_file, "r");
	if($handle)
	{
		while(($buffer=fgets($handle)) !== FALSE)
		{
			$row = explode("\t", trim($buffer));

			if(array_key_exists($row[1], $transcript_ids))	$ucsc2transcript[$row[0]] = $row[1];
		}
		fclose($handle);
	}
	else	trigger_error("Could not find file $transcript_file.",E_USER_ERROR);
	
	// extract cDNA information from UCSC-ids
	$handle = fopen($known_gene, "r");
	if($handle)
	{
		while(($buffer=fgets($handle)) !== FALSE)
		{
			$row = explode("\t", trim($buffer));

			if(array_key_exists($row[0], $ucsc2transcript))
			{
				$exon_start = explode(",", $row[8]);
				$exon_end = explode(",", $row[9]);
				$transcript_ids[$ucsc2transcript[$row[0]]]['strand'] = $row[2];
				$transcript_ids[$ucsc2transcript[$row[0]]]['cdsStart'] = $row[5];
				$transcript_ids[$ucsc2transcript[$row[0]]]['cdsEnd'] = $row[6];
				$transcript_ids[$ucsc2transcript[$row[0]]]['exons'] = array();
				for($i=0;$i<count($exon_start);++$i)
				{
					if(empty($exon_start[$i]))	continue;
					$transcript_ids[$ucsc2transcript[$row[0]]]['exons'][] = array('start'=>$exon_start[$i], 'end'=>$exon_end[$i]);
				}
			}
		}
		fclose($handle);
	}
	else trigger_error("Could not open file $known_gene.",E_USER_ERROR);

	// (4b) generate region filter to filter germline variant list
	$max_distance = 40;
	$regions = new Matrix();
	$splicing = array();
	$no_exons = array();
	foreach($transcript_ids as $id => $transcript_id)
	{
		if(!isset($transcript_id['exons']))
		{
			$no_exons[] = $id;
			continue;
		}
		foreach($transcript_id['variant_rows'] as $variant_row)
		{
			$somatic_variant = $file_som->getRow($variant_row);
			$region = array('',-1,-1,-1,-1);
			$region[0] = $somatic_variant[0];
			// get current exon
			$current_exon = -1;
			foreach($transcript_id['exons'] as $exon => $r)
			{
				// exonic variant
				if($somatic_variant[1]>=$r['start'] && $somatic_variant[2]<=$r['end'])	$current_exon = $exon;
			}
			if($current_exon==-1)
			{
				// splicing variant, uncomment following line if region of splicing variants is relevant, TODO: distinguish between splicing, intronic, intergenic, etc.
				// $regions->addRow(array($somatic_variant[0],$somatic_variant[1]-1,$somatic_variant[2],"",$id));
				$splicing[] = $variant_row;
				continue;
			}
			// get start
			$found = false;
			$tmp_start = $somatic_variant[1];
			$tmp_diff = $max_distance;
			$tmp_exon = $current_exon;
			while($found==FALSE)
			{
				$tmp_diff -= ($tmp_start - $transcript_id['exons'][$tmp_exon]['start']);
				if($tmp_diff>0)
				{
					--$tmp_exon;
					if(!isset($transcript_id['exons'][$tmp_exon]['start']))
					{
						trigger_error("Could not determine start point for variant in row $variant_row with starting position ".$somatic_variant[1].", using ".$transcript_id['exons'][$tmp_exon+1]['start']." as start point (first exon).",E_USER_NOTICE);
						$region[1] = $transcript_id['exons'][$tmp_exon+1]['start'];
						$region[3] = $tmp_exon+1;
						$found=TRUE;
						continue;
					}
					$tmp_start = $transcript_id['exons'][$tmp_exon]['end'];
				}
				else if($tmp_diff<=0)
				{
					$region[1] = $transcript_id['exons'][$tmp_exon]['start'] - $tmp_diff - 1;
					$region[3] = $tmp_exon;
					$found=TRUE;
				}
			}
			// get end
			$found = false;
			$tmp_end = $somatic_variant[2];
			$tmp_diff = $max_distance;
			$tmp_exon = $current_exon;
			while($found==FALSE)
			{
				$tmp_diff -= ($transcript_id['exons'][$tmp_exon]['end'] - $tmp_end);
				if($tmp_diff>0)
				{
					++$tmp_exon;
					if(!isset($transcript_id['exons'][$tmp_exon]['start']))
					{
						trigger_error("Could not determine endpoint for variant in row $variant_row with starting position ".$somatic_variant[1].", using last exon as end point.",E_USER_NOTICE);
						$region[2] = $transcript_id['exons'][$tmp_exon-1]['end'];
						$region[4] = $tmp_exon-1;
						$found=TRUE;
						continue;
					}
					$tmp_end = $transcript_id['exons'][$tmp_exon]['start'];
				}
				else if($tmp_diff<=0)
				{
					$region[2] = $transcript_id['exons'][$tmp_exon]['end'] + $tmp_diff;
					$region[4] = $tmp_exon;
					$found=TRUE;
				}
			}
			// extract all exons from start to end
			$variant_id = implode("_", array($somatic_variant[0],$somatic_variant[1],$somatic_variant[2],$somatic_variant[3],$somatic_variant[4]));
			for($exon=$region[3];$exon<=$region[4];++$exon)
			{
				$tmp_start = $transcript_id['exons'][$exon]['start'];
				$tmp_end = $transcript_id['exons'][$exon]['end'];
				if($tmp_start < $region[1])	$tmp_start = $region[1];
				if($tmp_end > $region[2])	$tmp_end = $region[2];
				$regions->addRow(array($region[0],$tmp_start, $tmp_end,$exon+1,$id,$variant_id));
			}
		}
	}
	if(!empty($no_exons))	trigger_error("No exons availabe for ".count($no_exons)." IDs: ".implode(", ",$no_exons),E_USER_NOTICE);
	if(!empty($splicing))	trigger_error("Skipped splicing variants in the following rows: ".implode(", ",array_unique($splicing)).". Variants were skipped.", E_USER_NOTICE);
	$out_reg = $parser->tempFile("_reg.bed");
	$regions->toTSV($out_reg);
	// $regions->toTSV( $o_folder.basename($nor_bam, ".bam").".reg");

	// (4c) find reference variants around coding variants and generate filtered reference variant file
	$tmp_folder = $parser->tempFolder("analyze");
	$nor_vcffile = $tmp_folder."/".$n_dna_id.".vcf.gz";
	$nor_varfile = $tmp_folder."/".$n_dna_id.".GSvar";
	$nor_bam = $p_folder."/Sample_".$n_dna_id."/".$n_dna_id.".bam";
	$extras_vc = "-target $out_reg ";
	$parser->execTool("NGS/vc_freebayes.php", "-bam $nor_bam -out $nor_vcffile -build ". $n_dna_sys['build']." $extras_vc");
	$parser->execTool("Pipelines/annotate.php", "-out_name $n_dna_id -out_folder $tmp_folder -system $t_dna_sys_file -vcf $nor_vcffile");

	// add new col that contains information about surrounding variants
	$new_col = array();
	$ref = Matrix::fromTSV($nor_varfile);
	for($i=0;$i<$file_som->rows();++$i)	// get all somatic refseq-ids
	{
		$row_som = $file_som->getRow($i);
		$variant_id = implode("_", array($row_som[0],$row_som[1],$row_som[2],$row_som[3],$row_som[4]));
		
		//identify regions that belong to this variant
		$target_regions = array();
		for($j=0;$j<$regions->rows();++$j)
		{
			$row_reg = $regions->getRow($j);
			if($row_reg[5]==$variant_id)
			{
				$target_regions[] = array($row_reg[0],$row_reg[1]+1,$row_reg[2]);
			}
		}

		// check if germline variants are in these regions
		$has_germline_variant = false;
		for($j=0;$j<$ref->rows();++$j)
		{
			$row_ref = $ref->getRow($j);
			
			foreach($target_regions as $target_region)
			{
				if($row_ref[0]==$target_region[0])
				{
					if(range_overlap($row_ref[1],$row_ref[2],$target_region[1],$target_region[2]))	$has_germline_variant = true;
				}
			}
		}

		if($has_germline_variant)	$new_col[] = "yes";
		else	$new_col[] = "no";
	}
	$file_som->addCol($new_col, "nearby_germline_variant","There are germline variants surrounding this variant (+/- $max_distance bp, new feature - use with caution, haplotype not considered).");
	$file_som->toTSV($annotated);
	
	// prepare IGV-session for all files
	if(is_dir($p_folder) && is_dir($o_folder))
	{
		$rel_path = relative_path($o_folder, $p_folder);
		$igv_session = array();
		$igv_session[] = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n";
		$igv_session[] = "<Session genome=\"1kg_v37\" hasGeneTrack=\"true\" hasSequenceTrack=\"true\" locus=\"all\" path=\".\" version=\"8\">\n";
		$igv_session[] = "    <Resources>\n";
		if(is_file($som_seg))   $igv_session[] = "\t\t<Resource path=\"".$rel_path."/".$o_folder."/".basename($som_seg)."\"/>\n";
		if(is_file($som_vcf))   $igv_session[] = "\t\t<Resource path=\"".$rel_path."/".$o_folder."/".basename($som_vcf)."\"/>\n";
		if(is_file($t_dna_bam))	$igv_session[] = "\t\t<Resource path=\"".$rel_path."/Sample_".basename($t_dna_bam,".bam")."/".basename($t_dna_bam)."\"/>\n";
		if(is_file($n_dna_bam))	$igv_session[] = "\t\t<Resource path=\"".$rel_path."/Sample_".basename($n_dna_bam,".bam")."/".basename($n_dna_bam)."\"/>\n";
		if(is_file($t_rna_bam))	$igv_session[] = "\t\t<Resource path=\"".$rel_path."/Sample_".basename($t_rna_bam,".bam")."/".basename($t_rna_bam)."\"/>\n";
		if(is_file($n_rna_bam))	$igv_session[] = "\t\t<Resource path=\"".$rel_path."/Sample_".basename($n_rna_bam,".bam")."/".basename($n_rna_bam)."\"/>\n";
		$igv_session[] = "    </Resources>\n";
		$igv_session[] = "    <HiddenAttributes>\n";
		$igv_session[] = "        <Attribute name=\"NAME\"/>\n";
		$igv_session[] = "        <Attribute name=\"DATA FILE\"/>\n";
		$igv_session[] = "        <Attribute name=\"DATA TYPE\"/>\n";
		$igv_session[] = "    </HiddenAttributes>\n";
		$igv_session[] = "</Session>";
		file_put_contents($o_folder."/".$t_dna_id."-".$n_dna_id."_igv.xml",implode("\n",$igv_session));
	}
	else
	{
		trigger_error("IGV-Session-File was not created. Folder $p_folder or $o_folder does not exist.",E_USER_WARNING);
	}
}
?>