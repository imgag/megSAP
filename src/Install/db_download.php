<?php
/** 
	@page db_download
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("db_download", "Downloads databases for various tools.");
$parser->addString("data_folder", "Path to megSAP data folder", true, repository_basedir()."/data");
$all_dbs = array("kraken2", "STAR", "Clair3");
$parser->addString("dbs", "Comma-separated list of dbs to download/install/build.", true, implode(",", $all_dbs));

extract($parser->parse($argv));

//check dbs
$dbs = explode(",", $dbs);
foreach($dbs as $db)
{
	if (!in_array($db, $all_dbs)) trigger_error("Unknown db '$db'!", E_USER_ERROR);
}

if (in_array("kraken2", $dbs))
{
	print "\n";
	print "### kraken2 ###\n";
	
	$kraken_data_dir = "{$data_folder}/dbs/kraken2_filter_hb";
    $hemoglobin_fa = "{$data_folder}/misc/human_hemoglobin_tx.fa";
	print "Downloading kraken2 database to {$kraken_data_dir} ...\n";
	
	//create kraken2 db folder if not existent
	if (!file_exists($kraken_data_dir))
	{
		exec2("mkdir -p $kraken_data_dir");
	}

    if (!file_exists($hemoglobin_fa))
    {
        trigger_error("human_hemoglobin_tx.fa not found at '{$hemoglobin_fa}'. Run 'download_dbs.sh' first.", E_USER_ERROR);
    }

	//build kraken2 database
    print "Building kraken2 database...\n";
	$parser->execApptainer("kraken2", "kraken2-build", "-db $kraken_data_dir --download-taxonomy --skip-map  --use-ftp", [$kraken_data_dir]);
	$parser->execApptainer("kraken2", "kraken2-build", "-db $kraken_data_dir --add-to-library {$hemoglobin_fa}", [$kraken_data_dir], [dirname($hemoglobin_fa)]);
	$parser->execApptainer("kraken2", "kraken2-build", "--build  --threads 5 --db $kraken_data_dir", [$kraken_data_dir]);
    print "Finished building kraken2 database\n";
}

if (in_array("STAR", $dbs))
{
	print "\n";
	print "### STAR ###\n";
	
	$star_genome_dir = "{$data_folder}/genomes/STAR/GRCh38";
    $star_gtf_file = "{$data_folder}/dbs/gene_annotations/GRCh38.gtf";
    $genome = "{$data_folder}/genomes/GRCh38.fa";

    if (!file_exists($genome))
    {
        trigger_error("Genome not found at '$genome'. Run 'download_GRCh38.sh' first.", E_USER_ERROR);
    }

    if (!file_exists($star_gtf_file))
    {
        trigger_error("GRCh38.gtf not found at '$star_gtf_file'. Run 'download_dbs_rna.sh' first.", E_USER_ERROR);
    }

	//create STAR genome folder if not existent
	if (!file_exists($star_genome_dir))
	{
		exec2("mkdir -p $star_genome_dir");
	}

    print "Indexing genome '{$genome}' using STAR. This may take a while.\n";
	$parser->execApptainer("STAR", "STAR", "--runThreadN 20 --runMode genomeGenerate --genomeDir {$star_genome_dir} --genomeFastaFiles {$genome} --sjdbGTFfile {$star_gtf_file}", [$star_genome_dir, $star_gtf_file, $data_folder]);
}

if (in_array("Clair3", $dbs))
{
	print "\n";
	print "### Clair3 ###\n";
	
    $model_path = get_path("clair3_models");

    //create folder for Clair3 models if not existent
    if (!file_exists($model_path))
	{
		exec2("mkdir -p $model_path");
	}

    //download models
    chdir($model_path);

    print "Downloading Clair3 models from https://github.com/nanoporetech/rerio.git\n";
    exec2("git clone https://github.com/nanoporetech/rerio.git");
    execApptainer("python", "python3", "{$model_path}/rerio/download_model.py --clair3", [$model_path]);
    exec2("mv rerio/clair3_models/* .");
    exec2("rm -rf rerio");

    print "Downloading Clair3 model 'r941_prom_sup_g5014' from 'http://www.bio8.cs.hku.hk/clair3/clair3_models/r941_prom_sup_g5014.tar.gz'\n";
    exec2("wget http://www.bio8.cs.hku.hk/clair3/clair3_models/r941_prom_sup_g5014.tar.gz");
    exec2("tar xzf r941_prom_sup_g5014.tar.gz");
    exec2("rm r941_prom_sup_g5014.tar.gz");

    print "Downloading Clair3 model 'r941_prom_sup_g360+g422' from 'http://www.bio8.cs.hku.hk/clair3/clair3_models/r941_prom_hac_g360+g422.tar.gz'\n";
    exec2("wget http://www.bio8.cs.hku.hk/clair3/clair3_models/r941_prom_hac_g360+g422.tar.gz");
    exec2("tar xzf r941_prom_hac_g360+g422.tar.gz");
    exec2("rm r941_prom_hac_g360+g422.tar.gz");
    exec2("mv ont/ r941_prom_hac_g360+g422/");
}
?>