<?php

/**
	@page annotate
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("annotate", "Annotate variants called on genome build GRCh37.");
$parser->addString("out_name", "Processed sample ID (e.g. 'GS120001_01').", false);
$parser->addString("out_folder", "Output folder.", false);
//optional
$parser->addInfile("system", "Processing system INI file (automatically determined from NGSD if 'out_name' is a valid processed sample name).", true);
$parser->addString("vcf", "Path to (bgzipped) VCF file (if different from {output_folder}/{out_name}_var.vcf.gz).", true, "");
$parser->addFlag("no_fc", "No format check (vcf/tsv).");
$parser->addFlag("multi", "Enable multi-sample mode.");
$parser->addFlag("somatic", "Enable somatic mode (no variant QC and no GSvar file).", true, "na");
$parser->addFlag("updown", "Don't discard up- or downstream anntations (5000 bases around genes).");
$parser->addInt("threads", "The maximum number of threads used.", true, 2);
extract($parser->parse($argv));

//input file names
$vcffile = $out_folder."/".$out_name."_var.vcf.gz";
if(!empty($vcf))
{
	$vcffile = $vcf;
}

//unzip input VCF if necessary
if (ends_with($vcffile, ".vcf"))
{
	$vcf_unzipped = $vcffile;
}
else
{
	$vcf_unzipped = $parser->tempFile("_unzipped.vcf");
	$parser->exec("bgzip", "-cd $vcffile > $vcf_unzipped", false); //no output logging, because Toolbase::extractVersion() does not return
}

//output file names
$annfile = $parser->tempFile(".vcf");
$annfile_zipped = $out_folder."/".$out_name."_var_annotated.vcf.gz";	
$varfile = $out_folder."/".$out_name.".GSvar";
$stafile = $out_folder."/".$out_name."_stats_vc.qcML";

//get system
$sys = load_system($system, $out_name);
if ($sys['build']!="GRCh37")
{
	trigger_error("Unknown genome build ".$sys['build']." cannot be annotated!", E_USER_ERROR);
}

//annotate VCF
$args = [];
$args[] = "-in ".$vcf_unzipped;
$args[] = "-out ".$annfile;
$args[] = "-build ".$sys['build'];
$args[] = "-threads ".$threads;
$args[] = "-ps_name ".$out_name;
if ($somatic) $args[] = "-somatic";
$parser->execTool("NGS/an_vep.php", implode(" ", $args));

//zip annotated VCF file
$parser->exec("bgzip", "-c $annfile > $annfile_zipped", false); //no output logging, because Toolbase::extractVersion() does not return
$parser->exec("tabix", "-p vcf $annfile_zipped", false); //no output logging, because Toolbase::extractVersion() does not return

//convert to GSvar file
if (!$somatic) //germline only
{
	//calculate variant statistics (after annotation because it needs the ID and ANN fields)
	$parser->exec(get_path("ngs-bits")."VariantQC", "-in $annfile -out $stafile", true);
	
	$args = array("-in $annfile", "-out $varfile", "-blacklist");
	if ($multi) $args[] = "-genotype_mode multi";
	if ($updown) $args[] = "-updown";
	if ($sys['type']=="WGS") $args[] = "-wgs";
	$parser->execTool("NGS/vcf2gsvar.php", implode(" ", $args));
	
	if(db_is_enabled("NGSD"))
	{
		//annotate variant allele frequency and depth from related RNA sample, if available
		//sample id from name
		$db = DB::getInstance("NGSD");
		list($sample_name, $ps_num) = explode("_", $out_name);
		$res = $db->executeQuery("SELECT id, name FROM sample WHERE name=:name", array("name" => $sample_name));
		if (count($res) >= 1)
		{
			$sample_id = $res[0]['id'];

			//related samples
			$res = $db->executeQuery("SELECT * FROM sample_relations WHERE relation='same sample' AND (sample1_id=:sid OR sample2_id=:sid)", array("sid" => $sample_id));

			//collect potential samples (sample IDs)
			$psamples_for_annotation = [];
			foreach ($res as $row)
			{
				$sample_id_annotation = $row['sample1_id'] != $sample_id ? $row['sample1_id'] : $row['sample2_id'];
				$res = $db->executeQuery("SELECT ps.sample_id, ps.process_id, ps.processing_system_id, ps.quality, sys.id, sys.type, CONCAT(s.name, '_', LPAD(ps.process_id, 2, '0')) as psample FROM processed_sample as ps, processing_system as sys, sample as s WHERE ps.sample_id=:sid AND sys.type='RNA' AND ps.processing_system_id=sys.id AND ps.sample_id=s.id AND (NOT ps.quality='bad')", array("sid" => $sample_id_annotation));
				$psamples_for_annotation = array_merge($psamples_for_annotation, array_column($res, 'psample'));
			}
			if (count($psamples_for_annotation) > 0)
			{
				//use BAM file from last entry
				$bam_rna = get_processed_sample_info($db, end($psamples_for_annotation))["ps_bam"];
				$parser->exec(get_path("ngs-bits")."VariantAnnotateFrequency", "-in {$varfile} -out {$varfile} -bam {$bam_rna} -depth -name rna", true);
			}
		}
	}
	
	//check output TSV file (not for somatic)
	if(!$no_fc)
	{
		$parser->execTool("NGS/check_tsv.php", "-in $varfile -build ".$sys['build']);
	}
}





?>
