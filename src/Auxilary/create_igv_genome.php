<?php
/** 
	@page create_igv_genome
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("create_igv_genome", "Creates a genome with Ensembl gene annotation for IGV (GRCh38 only).");
$parser->addString("output_folder", "Target to store JSON file and all required data", false);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
$parser->addString("name", "Name used for the genome. ", true, "Human (hg38) masked/Ensembl");
$parser->addString("id", "ID used for the genome. ", true, "hg38_ensembl_masked");
extract($parser->parse($argv));

//init
$data_folder = get_path("data_folder");

//create output folder
if (!file_exists($output_folder))
{
	exec2("mkdir -p $output_folder");
}

//copy genome to output folder
$genome_fasta = genome_fasta($build);
$parser->copyFile($genome_fasta, $output_folder."/".$build.".fa");
$parser->copyFile($genome_fasta.".fai", $output_folder."/".$build.".fa.fai");

//copy from data folder
$parser->copyFile("{$data_folder}/igv/cytoBandIdeo.txt.gz", "{$output_folder}/cytoBandIdeo.txt.gz");
$parser->copyFile("{$data_folder}/igv/hg38_alias.tab", "{$output_folder}/hg38_alias.tab");

//export gene track
print "Export Ensembl gene track from NGSD ...\n";
$all_transcripts_tmp = $parser->tempFile(".txt");
$mane_transcripts_tmp = $parser->tempFile(".txt");
$all_transcripts = "EnsemblGeneTrack.txt.gz";
$mane_transcripts = "EnsemblGeneTrack_Mane+Clinical.txt.gz";
$parser->execApptainer("ngs-bits", "NGSDExportIgvGeneTrack", "-out {$all_transcripts_tmp} -out_mane {$mane_transcripts_tmp}");
$parser->execApptainer("htslib", "bgzip", "-c {$all_transcripts_tmp} > {$output_folder}/{$all_transcripts}", [], [$output_folder]);
$parser->execApptainer("htslib", "bgzip", "-c {$mane_transcripts_tmp} > {$output_folder}/{$mane_transcripts}", [], [$output_folder]);

//create JSON
print "Create JSON ...\n";
$data = array();
$data["id"] = $id;
$data["name"] = $name;
$data["fastaURL"] = $build.".fa";
$data["indexURL"] = $build.".fa.fai";
$data["cytobandURL"] = "cytoBandIdeo.txt.gz";
$data["aliasURL"] = "hg38_alias.tab";
$data["chromosomeOrder"] = array("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrMT");
$data["tracks"] = array();

$ensembl_mane_track = array();
$ensembl_mane_track["name"] = "Ensembl Genes (MANE)";
$ensembl_mane_track["format"] = "refgene";
$ensembl_mane_track["url"] = $mane_transcripts;
$data["tracks"][] = $ensembl_mane_track;

$ensembl_track = array();
$ensembl_track["name"] = "Ensembl Genes (all)";
$ensembl_track["format"] = "refgene";
$ensembl_track["url"] = $all_transcripts;
$data["tracks"][] = $ensembl_track;


file_put_contents($output_folder."/".$id.".json", json_encode($data, JSON_PRETTY_PRINT));

print "... done\n";
?>
