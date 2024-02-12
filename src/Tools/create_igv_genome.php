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

//create output folder
if (!file_exists($output_folder))
{
	exec2("mkdir -p $output_folder");
}

//copy genome to output folder
$genome_fasta = genome_fasta($build);
$parser->copyFile($genome_fasta, $output_folder."/".$build.".fa");
$parser->copyFile($genome_fasta.".fai", $output_folder."/".$build.".fa.fai");

//download additional files
$url_list = array();
$url_list["cytoband"] = "https://s3.amazonaws.com/igv.org.genomes/hg38/annotations/cytoBandIdeo.txt.gz";
$url_list["alias"] = "https://s3.amazonaws.com/igv.org.genomes/hg38/hg38_alias.tab";
// $url_list["Refseq Genes"] = "https://s3.amazonaws.com/igv.org.genomes/hg38/ncbiRefSeq.sorted.txt.gz";
// $url_list["Refseq Genes index"] = "https://s3.amazonaws.com/igv.org.genomes/hg38/ncbiRefSeq.sorted.txt.gz.tbi";

foreach ($url_list as $track_name => $url) 
{
	print "Downloading '{$track_name}'...\n";
	$parser->exec("wget", "-O {$output_folder}/".basename($url)." ".$url);
}

//export gene track
$ngsbits = get_path("ngs-bits");
// $gff_temp = $parser->tempFile(".gff");
// $gff_temp2 = $parser->tempFile(".gff");
// $parser->exec($ngsbits."NGSDExportGff", "-out {$gff_temp}");
// $parser->exec("sort", "-k1,1V -k4,4n -k5,5rn -k3,3r {$gff_temp} > {$gff_temp2}");
// $gff_file = "NGSDEnsemblGeneExport.gff.gz";
// $parser->exec("bgzip", "-c {$gff_temp2} > {$output_folder}/{$gff_file}");
// $parser->exec("tabix", "-p gff {$output_folder}/{$gff_file}");
print "Export Ensembl gene track from NGSD ...\n";
$all_transcripts_tmp = $parser->tempFile(".txt");
$mane_transcripts_tmp = $parser->tempFile(".txt");
$all_transcripts = "EnsemblGeneTrack.txt.gz";
$mane_transcripts = "EnsemblGeneTrack_Mane+Clinical.txt.gz";
$parser->exec($ngsbits."NGSDExportIgvGeneTrack", "-out {$all_transcripts_tmp} -out_mane {$mane_transcripts_tmp}");
$parser->exec("bgzip", "-c {$all_transcripts_tmp} > {$output_folder}/{$all_transcripts}");
$parser->exec("bgzip", "-c {$mane_transcripts_tmp} > {$output_folder}/{$mane_transcripts}");

//export gene track for searching
print "Export Ensembl search track (hidden BED) from NGSD ...\n";
$bed_file = "NGSDEnsemblGeneExport.bed.gz";
$pipeline = array();
$pipeline[] = array($ngsbits."NGSDExportGenes", "");
$pipeline[] = array("cut", "-f1");
$pipeline[] = array($ngsbits."GenesToBed", "-source ensembl -mode gene -fallback");
$pipeline[] = array("bgzip", "-c > {$output_folder}/{$bed_file}");
$parser->execPipeline($pipeline, "create gene BED track");


//create JSON
print "Create JSON ...\n";
$data = array();
$data["id"] = $id;
$data["name"] = $name;
$data["fastaURL"] = $build.".fa";
$data["indexURL"] = $build.".fa.fai";
$data["cytobandURL"] = basename($url_list["cytoband"]);
$data["aliasURL"] = basename($url_list["alias"]);
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

// $ensembl_track = array();
// $ensembl_track["name"] = "Ensembl Genes";
// $ensembl_track["format"] = "gff3";
// $ensembl_track["url"] = $gff_file;
// $ensembl_track["indexURL"] = $gff_file.".tbi";

$ensembl_bed_track = array();
$ensembl_bed_track["name"] = "Ensembl Genes";
$ensembl_bed_track["format"] = "bed";
$ensembl_bed_track["url"] = $bed_file;
$ensembl_bed_track["hidden"] =  "true";
$data["tracks"][] = $ensembl_bed_track;

// $refseq_track = array();
// $refseq_track["name"] = "Refseq Genes";
// $refseq_track["format"] = "refgene";
// $refseq_track["url"] = basename($url_list["Refseq Genes"]);
// $refseq_track["indexURL"] = basename($url_list["Refseq Genes index"]);

file_put_contents($output_folder."/".$id.".json", json_encode($data, JSON_PRETTY_PRINT));

print "... done\n";
?>
