<?php

/**
  @page vc_straglr

*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_straglr", "Call repeat expansions with straglr. Creates an BED file.");
$parser->addInfile("in", "Input BAM file. Note: .bam.bai file is required!", false);
$parser->addOutfile("out", "Output VCF file.", false);
$parser->addInfile("loci", "BED file containing repeat loci.", false);
//optional
$parser->addInt("threads", "The maximum number of threads used.", true, 2);
$parser->addString("pid", "Processed sample name (e.g. 'GS120001_01'). If unset BAM file name will be used.", true);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
$parser->addFlag("test", "Run in test mode. Skips annotatation of NGSD information.");
extract($parser->parse($argv));

// use BAM file name as fallback if no processed sample name is provided
if(!isset($pid)) $pid = basename2($in);

//init
$out_folder = dirname($out);
$out_prefix = $out_folder."/".basename2($out);
$out_bed = $out_prefix.".bed";
$out_vcf = $out_prefix.".vcf";
$out_tsv = $out_prefix.".tsv";
$plot_folder = $out_folder."/repeat_expansions";
$genome = genome_fasta($build);

//get gender
$gender = "n/a";
if (db_is_enabled("NGSD") && !$test)
{
	$db = DB::getInstance("NGSD", false);
	$info = get_processed_sample_info($db, $pid, false);
	if (!is_null($info))
	{
		$gender = $info['gender'];
	}
}

// prepare command
$args = [];
$args[] = "$in";
$args[] = $genome;
$args[] = $out_prefix;
$args[] = "--loci {$loci}";
$args[] = "--nprocs $threads";
$args[] = "--sample {$pid}";

if ($gender != "n/a") $args[] = "--sex ".$gender[0];

//set bind paths for straglr container
$in_files = array();
$out_files = array();
$in_files[] = $in;
$in_files[] = $genome;
$in_files[] = $loci;
$out_files[] = $out_folder;

// run straglr container
$parser->execApptainer("straglr", "straglr.py", implode(" ", $args), $in_files, $out_files);

//get pathogenic ranges from NGSD
$min_pathogenic = array();
if (db_is_enabled("NGSD") && !$test)
{
	$db = DB::getInstance("NGSD", false);
	$result = $db->executeQuery("SELECT region, repeat_unit, min_pathogenic FROM `repeat_expansion`");
	foreach($result as $row)
	{
		if (!is_null($row["min_pathogenic"]) && $row["min_pathogenic"] != "")
		{
			$min_pathogenic[$row["region"]."_".$row["repeat_unit"]] = (int) $row["min_pathogenic"];
		}
	}
}

// annotate repeat names/ref motif
//read and index catalog
$loci_content = file($loci, FILE_IGNORE_NEW_LINES + FILE_SKIP_EMPTY_LINES);
$catalog = array();
foreach ($loci_content as $line) 
{
    if (starts_with($line, "#")) continue;
    list($chr, $start, $end, $motif, $repeat_id, $repeat_type, $ref_size, $ref_motif) = explode("\t", $line);
    $catalog["{$chr}:{$start}-{$end}"] = array($repeat_id, $repeat_type, $ref_size, $ref_motif);
}
//annotate catalog to output file
$bed_content_in = file($out_bed, FILE_IGNORE_NEW_LINES + FILE_SKIP_EMPTY_LINES);
$bed_content_out = array();
foreach ($bed_content_in as $line) 
{
    if (starts_with($line, "##")) $bed_content_out[] = $line;
    elseif (starts_with($line, "#")) $bed_content_out[] = $line."\trepeat_id\trepeat_type\tref_size\tref_motif\tmin_pathogenic";
    else 
    {
        list($chr, $start, $end, $motif, $repeat_id) = explode("\t", $line);
        $annotation = $catalog["{$chr}:{$start}-{$end}"];
        $ref_motif = $annotation[3];
        if (isset($min_pathogenic["{$chr}:{$start}-{$end}_{$ref_motif}"]))
        {
          $annotation[] = $min_pathogenic["{$chr}:{$start}-{$end}_{$ref_motif}"];
        }
        else 
        {
          $annotation[] = "";
        }
        $bed_content_out[] = $line."\t".implode("\t", $annotation);
    }
}
file_put_contents($out_bed, implode("\n", $bed_content_out));

//annotate reference motif to VCF
$vcf_content_in = file($out_vcf, FILE_IGNORE_NEW_LINES + FILE_SKIP_EMPTY_LINES);
$vcf_content_out = array();
foreach ($vcf_content_in as $line) 
{
    if (starts_with($line, "##")) $vcf_content_out[] = $line;
    elseif (starts_with($line, "#")) 
    {
      $vcf_content_out[] = "##INFO=<ID=REF_MOTIF,Number=.,Type=String,Description=\"Reference motif\">";
      $vcf_content_out[] = "##reference=".genome_fasta($build);
      $vcf_content_out[] = $line;
    }
    else 
    {
        $columns = explode("\t", $line);
        $chr = $columns[0];
        $start = $columns[1];
        $info =  explode(";", $columns[7]);
        foreach ($info as $kv_pair)
        {
          if (starts_with($kv_pair, "END=")) 
          {
            $end = explode("=", $kv_pair)[1];
            break;
          }
        }
        if (!isset($end)) trigger_error("No end position found in line '{$line}'!", E_USER_ERROR);
        $info[] = "REF_MOTIF=".$catalog["{$chr}:{$start}-{$end}"][3];
        $columns[7] = implode(";", $info);
        $vcf_content_out[] = implode("\t", $columns);

    }
}
file_put_contents($out_vcf, implode("\n", $vcf_content_out));



//sort output files
$parser->execApptainer("ngs-bits", "BedSort", "-in {$out_bed} -out {$out_bed}", [$out_bed]);
$parser->execApptainer("ngs-bits", "VcfSort", "-in {$out_vcf} -out {$out_vcf}", [$out_vcf]);


//create plots:
if (file_exists($plot_folder)) $parser->exec("rm", "-r {$plot_folder}"); //delete previous plots
mkdir($plot_folder);
$parser->execApptainer("straglrOn", "straglron.py", "{$out_bed} {$out_tsv} {$loci} -o {$plot_folder} --hist --alleles --bam {$in} --genome ".genome_fasta($build), [$out_folder, $loci, genome_fasta($build)]);



?>