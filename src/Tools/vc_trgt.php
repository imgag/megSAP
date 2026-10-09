<?php

/**
  @page vc_trgt

*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_trgt", "Call repeat expansions with trgt (PacBio only).");
$parser->addInfile("in", "Input BAM file. Note: .bam.bai file is required!", false);
$parser->addOutfile("out", "Output VCF file.", false);
$parser->addInfile("loci", "BED file containing repeat loci.", false);
$parser->addEnum("gender", "Gender of the sample.", false, ["female", "male", "n/a"]);
//optional
$parser->addInt("threads", "The maximum number of threads used.", true, 2);
$parser->addString("sample_name", "Processed sample name (e.g. 'GS120001_01'). If unset BAM file name will be used.", true);
$parser->addString("build", "The genome build to use.", true, "GRCh38");

extract($parser->parse($argv));

// use BAM file name as fallback if no processed sample name is provided
if(!isset($sample_name)) $sample_name = basename2($in);

//init
$out_folder = dirname($out);
$out_prefix = $out_folder."/".basename2($out);
$out_vcf = $out_prefix.".vcf";
$plot_folder = $out_folder."/repeat_expansions/";
$out_bam = "{$plot_folder}{$sample_name}_spanning.bam";
$genome = genome_fasta($build);


//genotype
$tmp_folder = $parser->tempFolder("_trgt");
$args = [];
$args[] = "--genome {$genome}";
$args[] = "--repeats {$loci}";
$args[] = "--reads {$in}";
$args[] = "--output-prefix {$tmp_folder}/{$sample_name}";
if ($gender == "male") $args[] = "--karyotype XY";
if ($gender == "female") $args[] = "--karyotype XX";
$args[] = "--threads {$threads}";
$args[] = "--sample-name {$sample_name}";
$args[] = "--output-type v"; //write unzipped vcf (makes post-processing easier)

$in_files = [$genome, $loci, $in];
$parser->exec("mkdir", "-p {$plot_folder}");
$parser->execApptainer("trgt", "trgt genotype", implode(" ", $args), $in_files);

//adapt vcf header
$tmp_vcf = "{$tmp_folder}/{$sample_name}.vcf";
$comments = vcf_load_header($tmp_vcf)[0];
$trgt_version = "";
foreach ($comments as $line) 
{
  if (starts_with($line, "##trgtVersion"))
  {
    //##trgtVersion=5.1.0-ec66463
    $trgt_version = trim(explode("=", $line)[1]);
    break;
  }
}
$comments[] = "##source=trgt V".$trgt_version;
$comments[] = "##fileDate=".date("Ymd");
$comments = vcf_sort_comments($comments);
vcf_replace_comments($tmp_vcf, $comments);

//annotate ref motif
$loci_buffer = file($loci, FILE_SKIP_EMPTY_LINES | FILE_IGNORE_NEW_LINES);
$ref_motifs = [];
foreach ($loci_buffer as $line) 
{
  $kv_pairs = explode(";", explode("\t", $line)[3]);
  $repeat_id = "";
  $ref_motif = "";
  foreach ($kv_pairs as $kv_pair) 
  {
    if (starts_with($kv_pair, "ID=")) $repeat_id = trim(explode("=", $kv_pair)[1]);
    if (starts_with($kv_pair, "RefMotif=")) $ref_motif = trim(explode("=", $kv_pair)[1]);
    if (($ref_motif != "") && ($repeat_id != ""))
    {
      $ref_motifs[$repeat_id] = $ref_motif;
      break;
    }
  }
  if (!isset($ref_motifs[$repeat_id])) trigger_error("'ID' of 'RefMotif' entry  missing in line '{$line}'!", E_USER_ERROR);
}
$vcf_buffer = file($tmp_vcf, FILE_SKIP_EMPTY_LINES | FILE_IGNORE_NEW_LINES);
$out_buffer = [];
foreach ($vcf_buffer as $line) 
{
  if (starts_with($line, "##contig=<ID=")) continue;
  if ($line[0] == "#") $out_buffer[] = $line;
  else
  {
    $parts = explode("\t", $line);
    $repeat_id = trim(explode("=", explode(";", $parts[7])[0])[1]);
    $parts[7] = $parts[7].";RefMotif=".$ref_motifs[$repeat_id];
    $out_buffer[] = implode("\t", $parts);
  }
}
file_put_contents($tmp_vcf, implode("\n", $out_buffer));

//sort results
$parser->execApptainer("ngs-bits", "VcfSort", "-in {$tmp_vcf} -out {$out_vcf}", [], [dirname($out_vcf)]);
$samtools_tmp = $parser->tempFolder("_samtools");
$parser->execApptainer("samtools", "samtools sort", "-T {$samtools_tmp}/sort -@ {$threads} -m 1G --reference {$genome} -o {$out_bam} {$tmp_folder}/{$sample_name}.spanning.bam", [$genome], [dirname($out_bam)]);
$parser->execApptainer("samtools", "samtools index", "-@ {$threads} {$out_bam}", [], [dirname($out_bam)]);

//generate plots

//prepare parameters
$args = [];
$args[] = "--genome {$genome}";
$args[] = "--repeats {$loci}";
$args[] = "--vcf {$out_vcf}";
$args[] = "--spanning-reads {$out_bam}";
$in_files = [$genome, $loci, $out_vcf, ];


foreach ($loci_buffer as $line) 
{
  $kv_pairs = explode(";", explode("\t", $line)[3]);
  $repeat_id = "";
  foreach ($kv_pairs as $kv_pair) 
  {
    if (!starts_with($kv_pair, "ID=")) continue;
    $repeat_id = trim(explode("=", $kv_pair)[1]);
    break;
  }
  if ($repeat_id == "") trigger_error("No 'ID' entry found in line '{$line}'!", E_USER_ERROR);

  $plot_file_name = "{$plot_folder}{$sample_name}_repeats_{$repeat_id}.svg";  

  $parser->execApptainer("trgt", "trgt plot", implode(" ", $args)." --image {$plot_file_name} --repeat-id {$repeat_id}", $in_files, [$plot_folder], false, true, false);

}


?>