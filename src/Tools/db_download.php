<?php
/** 
	@page db_download
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("db_download", "Downloads databases for various tools.");
$parser->addString("data_folder", "Path to megSAP data folder", true, repository_basedir()."/data");
$all_dbs = array("vep", "kraken2", "STAR");
$parser->addString("dbs", "Comma-separated list of dbs to download/install/build.", true, implode(",", $all_dbs));

extract($parser->parse($argv));

if (in_array("vep", $dbs))
{
    $vep_data_dir = "{$data_folder}/dbs/ensembl-vep-112";

    // create VEP cache data folder
    if (!file_exists($vep_data_dir))
    {
        exec2("mkdir -p $vep_data_dir");
    }

    // ensure that ftp and cache subdirectories exist
    if (!file_exists("$vep_data_dir/ftp"))
    {
        exec2("mkdir -p $vep_data_dir/ftp");
    }

    if (!file_exists("$vep_data_dir/cache"))
    {
        exec2("mkdir -p $vep_data_dir/cache");
    }

    // download VEP cache data
    $url = "https://ftp.ensembl.org/pub/release-112/variation/indexed_vep_cache/homo_sapiens_vep_112_GRCh38.tar.gz";
    $parser->exec("wget", "-O {$vep_data_dir}/ftp/ ".$url);

    // install ensembl-vep
    $parser->execSingularity("vep", get_path("container_vep"), "INSTALL.pl", "--SPECIES homo_sapiens --ASSEMBLY GRCh38 --AUTO ac --NO_UPDATE --NO_BIOPERL --CACHEDIR $vep_data_dir/cache --CACHEURL $vep_data_dir/ftp --NO_TEST --NO_HTSLIB", [$vep_data_dir]);
}

if (in_array("kraken2", $dbs))
{
	$kraken_data_dir = "{$data_folder}/dbs/kraken2_filter_hb";

	//create kraken2 db folder if not existent
	if (!file_exists($kraken_data_dir))
	{
		exec2("mkdir -p $kraken_data_dir");
	}

	//build kraken2 database
	$parser->execSingularity("kraken2", get_path("container_kraken2"), "kraken2-build", "-db $kraken_data_dir --download-taxonomy --skip-map  --use-ftp", [$kraken_data_dir]);
	$parser->execSingularity("kraken2", get_path("container_kraken2"), "kraken2-build", "-db $kraken_data_dir --add-to-library {$data_folder}/misc/human_hemoglobin_tx.fa", [$kraken_data_dir], ["{$data_folder}/misc/human_hemoglobin_tx.fa"]);
	$parser->execSingularity("kraken2", get_path("container_kraken2"), "kraken2-build", "--build  --threads 5 --db $kraken_data_dir", [$kraken_data_dir]);
}

if (in_array("STAR", $dbs))
{
	$star_genome_dir = "{$data_folder}/genomes/STAR/GRCh38";

	//create STAR genome folder if not existent
	if (!file_exists($star_genome_dir))
	{
		exec2("mkdir -p $star_genome_dir");
	}

	$parser->execSingularity("STAR", get_path("container_STAR"), "STAR", "--runThreadN 20 --runMode genomeGenerate --genomeDir {$star_genome_dir} --genomeFastaFiles {$data_folder}/genomes/GRCh38.fa --sjdbGTFfile {$data_folder}/dbs/gene_annotations/GRCh38.gtf");
}

?>