<?php 
/** 
	@page an_somatic_rna
	
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("an_somatic_rna", "RNA annotation for somatic pipeline.");
$parser->addInfile("t_bam", "Tumor sample BAM file.", false);
$parser->addString("full_prefix", "Full filepath prefix for out files", false);

$parser->addInfile("system", "Processing system file used for tumor DNA sample (resolved from NGSD via tumor BAM by default).", true);
$parser->addInfile("t_rna_bam", "Tumor RNA sample BAM file.", true);
$parser->addString("rna_ref_tissue", "Reference data for RNA annotation.", true);
$parser->addInt("steps", "number of steps being executed with the pipeline (e.g. 4 for 'vc,cn,an,an_rna').", true, 1);
$parser->addFlag("skip_correlation", "Skip sample correlation check.");

extract($parser->parse($argv));

//init
$t_id = basename2($t_bam);
$t_basename = dirname($t_bam)."/".$t_id;
$sys = load_system($system, $t_id);
$genome = genome_fasta($sys['build']);

$variants_gsvar = $full_prefix . ".GSvar"; // GSvar variants
$som_clincnv = $full_prefix . "_clincnv.tsv"; //ClinCNV output file

list($t_name) = explode("_", $t_id);

//Determine tumor rna sample from NGSD via sample_relations
$ps_rna_bams = array();

if( db_is_enabled("NGSD") && !isset($t_rna_bam) )
{
	$db = DB::getInstance("NGSD");
	$res = $db->executeQuery("SELECT id, name FROM sample WHERE name=:name", array("name" => $t_name));
	$sample_id = $res[0]["id"];
	
	//get related RNA samples from NGSD
	if (count($res) >= 1)
	{
		$relations = $db->executeQuery("SELECT * FROM sample_relations WHERE relation='same sample' AND (sample1_id=:sid OR sample2_id=:sid)", array("sid" => $sample_id));
		foreach($relations as $row)
		{
			$sample_id_annotation = $row['sample1_id'] != $sample_id ? $row['sample1_id'] : $row['sample2_id'];
			
			$res = $db->executeQuery("SELECT ps.sample_id, ps.process_id, ps.processing_system_id, ps.quality, sys.id, sys.type, CONCAT(s.name, '_', LPAD(ps.process_id, 2, '0')) as psample FROM processed_sample as ps, processing_system as sys, sample as s WHERE ps.sample_id=:sid AND sys.type='RNA' AND ps.processing_system_id=sys.id AND ps.sample_id=s.id", array("sid" => $sample_id_annotation));
			$rna_ids = array_column($res, 'psample');
			foreach($rna_ids as $rna_id)
			{
				$ps_rna_bams[$rna_id] = get_processed_sample_info($db, $rna_id)["ps_bam"];
			}
		}
	}
}
elseif(isset($t_rna_bam))
{
	$ps_rna_bams[basename($t_rna_bam)] = $t_rna_bam; 
}

if(count($ps_rna_bams) < 1)
{
	if ($steps > 1)
	{
		trigger_error("Skipping step an_rna!\nCouldn't find tumor RNA bam file. For annotation step \"an_rna\" tumor RNA bam file must be specified (via paramter -t_rna_bam or determined via sample_relations).", E_USER_WARNING);
		
		// exit an_somatic_rna.php
		return;
	}
	else
	{
		trigger_error("Couldn't find tumor RNA bam file. For annotation step \"an_rna\" tumor RNA bam file must be specified (via paramter -t_rna_bam or determined via sample_relations).", E_USER_ERROR);
	}
}

//RNA annotation
//Determine reference tissue type 1.) from parameter -rna_ref_tissue or 2.) from NGSD (if available)
if(!isset($rna_ref_tissue) && !db_is_enabled("NGSD"))
{
	trigger_error("For annotation step \"an_rna\" a tissue type for RNA reference data has to be specified.", E_USER_ERROR);
}
if(!isset($rna_ref_tissue) && db_is_enabled("NGSD"))
{
	$db = DB::getInstance("NGSD");
	$tumor_db_id = get_processed_sample_id($db, $t_id);
	$s_id = $db->getValue("SELECT sample_id FROM processed_sample where id=$tumor_db_id");
	$res = $db->getValues("SELECT DISTINCT sdi.disease_info FROM sample s LEFT JOIN sample_relations sr ON s.id=sr.sample1_id OR s.id=sr.sample2_id LEFT JOIN sample_disease_info sdi ON sdi.sample_id=sr.sample1_id OR sdi.sample_id=sr.sample2_id WHERE s.id=$s_id AND sdi.type='RNA reference tissue' AND (sr.relation='same sample' OR sr.relation IS NULL)");
	if(count($res) == 1)
	{
		list($rna_ref_tissue) = $res;
	}
	else
	{
		if (count($res) > 1)
		{
			trigger_error("Found multiple RNA reference tissue in NGSD. Aborting...", E_USER_ERROR);
		}
		else
		{
			trigger_error("Found no RNA reference tissue in NGSD. Skipping reference tissue annotation...", E_USER_WARNING);
		}
	}
}

$rna_counts = array(); //file contains transcript counts
foreach($ps_rna_bams as $rna_id => $rna_bam)
{
	$rna_counts_tmp = glob(dirname($rna_bam)."/*_counts.tsv");
	if(count($rna_counts_tmp) != 1)
	{
		trigger_error("Could not find or found multiple RNA count files in sample folder " . dirname($rna_bam), E_USER_ERROR);
	}
	
	$rna_counts[$rna_id] = $rna_counts_tmp[0];
}

//Annotate data from all detected RNA files
foreach($ps_rna_bams as $rna_id => $rna_bam)
{
	check_genome_build($rna_bam, $sys['build']);
	$rna_count = $rna_counts[$rna_id];
	
	//Calculate sample similarity between tumor and RNA
	if ($skip_correlation)
	{
		trigger_error("The correlation calculation between DNA and RNA ({$rna_id}) is skipped!");
		continue;
	}
	$min_corr = 0.85;  //TODO: evaluate value
	
	if (file_exists($t_bam))
	{
		$output = $parser->execApptainer("ngs-bits", "SampleSimilarity", "-in {$rna_bam} {$t_bam} -mode bam -ref {$genome} -build ".ngsbits_build($sys['build']), [$rna_bam, $t_bam]);
		$correlation = explode("\t", $output[0][1])[3];
		if ($correlation < $min_corr && ! $skip_correlation)
		{
			trigger_error("The genotype correlation of DNA and RNA ({$rna_id}) is {$correlation}; it should be above {$min_corr}!", E_USER_ERROR);
		}
		else
		{
			trigger_error("The genotype correlation of DNA and RNA ({$rna_id}) is {$correlation}.", E_USER_NOTICE);
		}
	}
	else
	{
		trigger_error("BAM file does not exist for tumor DNA sample '{$t_bam}'!", E_USER_ERROR);
	}
	
	//SNVs
	$args = [
	"-gsvar_in $variants_gsvar",
	"-out $variants_gsvar",
	"-rna_id $rna_id",
	"-rna_counts $rna_count",
	"-rna_bam $rna_bam",
	"-rna_target NGSD"
	];
	$parser->execTool("Tools/an_somatic_gsvar.php", implode(" ", $args));
	
	//CNVs
	if(file_exists($som_clincnv))
	{
		$args = [
			"-cnv_in $som_clincnv" ,
			"-out $som_clincnv",
			"-rna_counts $rna_count",
			"-rna_id $rna_id"
		];
		
		if (isset($rna_ref_tissue)) $args[] = "-rna_ref_tissue " .str_replace(" ", 0, $rna_ref_tissue);
		
		$parser->execTool("Tools/an_somatic_cnvs.php",  implode(" ", $args));
	}
}

//Reference tissue SNVs
$args = [
	"-gsvar_in $variants_gsvar",
	"-out $variants_gsvar",
];

if (isset($rna_ref_tissue)) $args[] = "-rna_ref_tissue " .str_replace(" ", 0, $rna_ref_tissue); //Replace spaces by 0 because it is diffcult to pass spaces via command line.

$parser->execTool("Tools/an_somatic_gsvar.php", implode(" ", $args));

?>
