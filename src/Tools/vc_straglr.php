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
$parser->addFlag("include_partials", "Detect and report reads only capturing partial repeats when genotyping."); //TODO Marc/Leon: make this the default?!
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
$args[] = "--max_num_clusters 2";

if ($gender != "n/a") $args[] = "--sex ".$gender[0];

if ($include_partials) $args[] = "--include_partials";

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
    $catalog["{$chr}:{$start}"] = array($repeat_id, $repeat_type, $ref_size, $ref_motif, $end);
}
//annotate catalog to output file
$bed_content_in = file($out_bed, FILE_IGNORE_NEW_LINES + FILE_SKIP_EMPTY_LINES);
$bed_content_out = array();
$wildtype_length = array();
foreach ($bed_content_in as $line) 
{
    if (starts_with($line, "##")) $bed_content_out[] = $line;
    elseif (starts_with($line, "#")) $bed_content_out[] = $line."\trepeat_id\trepeat_type\tref_size\tref_motif\tmin_pathogenic";
    else 
    {
        list($chr, $start, $end, $motif) = explode("\t", $line);
        $annotation = $catalog["{$chr}:{$start}"];
        array_pop($annotation);
        $ref_motif = $annotation[3];
        $repeat_id = $annotation[0];
        if (isset($min_pathogenic["{$chr}:{$start}-{$end}_{$ref_motif}"]))
        {
          $annotation[] = $min_pathogenic["{$chr}:{$start}-{$end}_{$ref_motif}"];
        }
        else 
        {
          $annotation[] = "";
        }
        $bed_content_out[] = $line."\t".implode("\t", $annotation);

        //store wildtype length (2nd allele on het and 1st allele on hom wt)
        $tmp_column = explode("\t", $line);
        if ($tmp_column[8] != "-") $wildtype_length[$repeat_id] = $tmp_column[8];
        else $wildtype_length[$repeat_id] = explode("\t", $line)[5];
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
      $vcf_content_out[] = "##INFO=<ID=RUC_WT,Number=.,Type=Float,Description=\"Wild type repeat unit count of corresponding repeat sequence\">";
      $vcf_content_out[] = "##reference=".genome_fasta($build);
      $vcf_content_out[] = $line;
    }
    else 
    {
        $columns = explode("\t", $line);
        $chr = $columns[0];
        $start = $columns[1];
        
        $info =  explode(";", $columns[7]);
        // $info[] = "END=".$catalog["{$chr}:{$start}"][4]; //add end pos
        $ref_motif = $catalog["{$chr}:{$start}"][3];
        $repeat_id = $catalog["{$chr}:{$start}"][0];
        $info[] = "REF_MOTIF=".$ref_motif;

        if (!isset($wildtype_length[$repeat_id])) trigger_error("No wiltype copy number found for repeat_id '{$repeat_id}'!", E_USER_ERROR);
        $wt_cn = $wildtype_length[$repeat_id];
        
        $format_values = explode(":", $columns[9]);
        //fix genotype
        //special case: male on chrX/Y:
        if (($gender == "male") && (($chr == "chrX" || $chr == "chrY")))
        {
          //there should be only one allele
          if (str_contains($format_values[0], "/")) trigger_error("Invalid genotype '".$format_values[0]."' in line '{$line}'!", E_USER_ERROR);
          // add wt copynumber
          if ($format_values[0] == "0") $info[] = "RUC_WT=".$wt_cn;
        }
        else
        {
          if (($format_values[0] == "0") || ($format_values[0] == "0/0"))
          {
            //hom wt:
            $format_values[0] = "0/0";
            $info[] = "RUC_WT=".$wt_cn;
          }
          else if ($format_values[0] == "1")
          {
            //hom alt
            $format_values[0] = "1/1";
          }
          else if ($format_values[0] == "0/1")
          {
            //het alt
            $info[] = "RUC_WT=".$wt_cn;
          }
          else if ($format_values[0] == "1/2")
          {
            //hom alt with 2 alleles
          }
          else trigger_error("Invalid genotype '".$format_values[0]."' in line '{$line}'!", E_USER_ERROR);
        }

        $columns[7] = implode(";", $info);
        $columns[9] = implode(":", $format_values);
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

$parser->execApptainer("straglrOn", "straglron.py", "{$out_bed} {$out_tsv} {$loci} -o {$plot_folder} --hist --alleles --bam {$in} --genome ".genome_fasta($build), [$in, $out_folder, $loci, genome_fasta($build)]);



?>