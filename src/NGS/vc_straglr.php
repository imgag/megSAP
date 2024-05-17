<?php

/**
  @page vc_straglr

*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_straglr", "Call repeat expansions with straglr. Creates an BED file.");
$parser->addInfile("in", "Input BAM file.", false);
$parser->addOutfile("out", "Output BED file.", false);
$parser->addInfile("loci", "BED file containing repeat loci.", false);
//optional
$parser->addInt("threads", "The maximum number of threads used.", true, 2);
$parser->addString("pid", "Processed sample name (e.g. 'GS120001_01'). If unset BAM file name will be used.", true);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
extract($parser->parse($argv));

// use BAM file name as fallback if no processed sample name is provided
if(!isset($pid)) $pid = basename2($in);

//init
$out_prefix = dirname($out)."/".basename2($out);
$out_vcf = $out_prefix.".vcf";
$straglr = get_path("straglr");

//get gender
$gender = "n/a";
if (db_is_enabled("NGSD"))
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
$args[] = $straglr;
$args[] = "--loci {$loci}";
$args[] = "--nprocs $threads";
$args[] = "--sample {$pid}";
$args[] = "$in";
$args[] = genome_fasta($build);
$args[] = $out_prefix;
if ($gender != "n/a") $args[] = "--sex ".$gender[0];
		
// prepare PATH
putenv("PATH=".getenv("PATH").PATH_SEPARATOR.dirname(get_path("trf")).PATH_SEPARATOR.dirname(get_path("blastn")));
print "\n........\n".getenv("PATH")."\n........\n";
// run straglr
$parser->exec(get_path("python3"), implode(" ", $args));

// annotate repeat names
//read and index catalog
$loci_content = file($loci, FILE_IGNORE_NEW_LINES + FILE_SKIP_EMPTY_LINES);
$catalog = array();
foreach ($loci_content as $line) 
{
    if (starts_with($line, "#")) continue;
    list($chr, $start, $end, $motive, $repeat_id, $repeat_type, $ref_size) = explode("\t", $line);
    $catalog["{$chr}:{$start}-{$end}"] = array($repeat_id, $repeat_type, $ref_size);
}
//annotate catalog to output file
$bed_content_in = file($out, FILE_IGNORE_NEW_LINES + FILE_SKIP_EMPTY_LINES);
$bed_content_out = array();
foreach ($bed_content_in as $line) 
{
    if (starts_with($line, "##")) $bed_content_out[] = $line;
    elseif (starts_with($line, "#")) $bed_content_out[] = $line."\trepeat_id\trepeat_type\tref_size";
    else 
    {
        list($chr, $start, $end, $motive, $repeat_id) = explode("\t", $line);
        $bed_content_out[] = $line."\t".implode("\t", $catalog["{$chr}:{$start}-{$end}"]);
    }
}
file_put_contents($out, implode("\n", $bed_content_out));

//sort output files
$parser->exec(get_path("ngs-bits")."BedSort", "-in {$out} -out {$out}");
$parser->exec(get_path("ngs-bits")."VcfSort", "-in {$out_vcf} -out {$out_vcf}");



//TODO: further post-processing?

// //get tool version
// list($stdout, $stderr, $ec) = $parser->exec(get_path("python3"), $straglr." --version");
// $tool_version = "straglr ".trim($stdout[0]);

// //convert to VCF
// $vcf = $out_prefix.".vcf";
// $fh = fopen2($vcf, 'w');

// //write header
// fwrite($fh, "##filedate=".date("Y-m-d")."\n");
// fwrite($fh, "##source={$tool_version}\n");
// fwrite($fh, "##reference=".genome_fasta($build)."\n");

// //add info headers
// fwrite($fh, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">\n");
// fwrite($fh, "##INFO=<ID=REF,Number=1,Type=Integer,Description=\"Reference copy number\">\n");
// fwrite($fh, "##INFO=<ID=REPID,Number=1,Type=String,Description=\"Repeat identifier as specified in the variant catalog\">\n");
// fwrite($fh, "##INFO=<ID=RL,Number=1,Type=Integer,Description=\"Reference length in bp\">\n");
// fwrite($fh, "##INFO=<ID=RU,Number=1,Type=String,Description=\"Repeat unit in the reference orientation\">\n");

// //add format headers
// fwrite($fh, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
// fwrite($fh, "##FORMAT=<ID=REPCN,Number=1,Type=String,Description=\"Number of repeat units spanned by the allele\">\n");
// fwrite($fh, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n");
// fwrite($fh, "##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Number of supporting reads consistent with the allele\">\n");
// //fwrite($fh, "##INFO=<ID=RU,Number=1,Type=String,Description=\"Repeat unit in the reference orientation\">\n");


// fwrite($fh, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{$pid}\n");

// //convert data lines
// $bed_content = Matrix::fromTSV($out);
// $i_chr = $bed_content->getColumnIndex("chrom");
// $i_start = $bed_content->getColumnIndex("start");
// $i_end = $bed_content->getColumnIndex("end");
// $i_repeat_unit = $bed_content->getColumnIndex("repeat_unit");
// $i_a1_size = $bed_content->getColumnIndex("allele1:size");
// $i_a1_cn = $bed_content->getColumnIndex("allele1:copy_number");
// $i_a1_support = $bed_content->getColumnIndex("allele1:support");
// $i_a2_size = $bed_content->getColumnIndex("allele2:size");
// $i_a2_cn = $bed_content->getColumnIndex("allele2:copy_number");
// $i_a2_support = $bed_content->getColumnIndex("allele2:support");
// $i_repeat_id = $bed_content->getColumnIndex("repeat_id");
// $i_ref_size = $bed_content->getColumnIndex("ref_size");
// for ($i=0; $i < $bed_content->rows(); $i++) 
// { 
//     $vcf_line = array();
//     $chr = trim($bed_content->get($i, $i_chr));
//     $vcf_line[] = $chr; //CHROM
//     $start = ((int) $bed_content->get($i, $i_start)) + 1; 
//     $vcf_line[] = $start; //POS
//     $vcf_line[] = "."; //ID
//     $vcf_line[] = get_ref_seq($build, $chr, $start, $start); //REF
//     $vcf_line[] = "."; //ALT
//     $vcf_line[] = "."; //QUAL
//     $vcf_line[] = ""; //FILTER
    
//     //INFO
//     $info_column = array();
//     $info_column[] = "END=".trim($bed_content->get($i, $i_end));
//     $ref_size = trim($bed_content->get($i, $i_ref_size));
//     $info_column[] = "REF=".$ref_size;
//     $info_column[] = "RU=".trim($bed_content->get($i, $i_repeat_unit));
//     $info_column[] = "REPID=".trim($bed_content->get($i, $i_repeat_id));
//     $info_column[] = "RL=".($bed_content->get($i, $i_ref_size) * strlen($bed_content->get($i, $i_repeat_unit)));
//     $vcf_line[] = implode(";", $info_column);

//     //FORMAT/SAMPLE
//     $is_hemizygous = false;
//     $a1_size = $bed_content->get($i, $i_a1_size);
//     $a1_cn = $bed_content->get($i, $i_a1_cn);
//     $a1_support = (int) $bed_content->get($i, $i_a1_size);
//     $a1_gt = (($a1_cn == $ref_size)?0:1);
//     if (trim($bed_content->get($i, $i_a1_size)) == "-")
//     {
//         if (($chr == "chrY") || (($chr == "chrX") && ($gender == "male")))
//         {
//             //hemizygous
//             $is_hemizygous = true;
//             continue;   
//         }
//         //else: homozygous -> share support reads
//         $a2_size = $a1_size;
//         $a2_cn = $a1_cn;
//         $a2_support = $a1_support / 2;
//         $a1_support = $a2_support;
//         $a2_gt = $a1_gt;
//     }
//     else
//     {
//         //heterozygous
//         $a2_size = (float) $bed_content->get($i, $i_a1_size);
//         $a2_cn = (float) $bed_content->get($i, $i_a1_cn);
//         $a2_support = (int) $bed_content->get($i, $i_a1_size);
//         if ($a2_cn == $ref_size) $a2_gt = 0;
//         elseif ($a1_gt == 0) $a2_gt = 1;
//         else $a2_gt = 2;
//     }

//     $format_values = array();
    
//     if ($is_hemizygous)
//     {
//         $format_values["GT"] = $a1_gt; //GT 
//         $format_values["REPCN"] = $a1_cn; //REPCN 
//         $format_values["AD"] = $a1_support; //AD 
//         $format_values["DP"] = $a1_support; //DP
//     } 
//     else 
//     {
//         $format_values["GT"] = $a1_gt."/".$a2_gt; //GT
//         $format_values["REPCN"] = $a1_cn."/".$a2_cn; //REPCN
//         $format_values["AD"] = $a1_support."/".$a2_support; //AD
//         $format_values["DP"] = $a1_support + $a2_support; //DP
//     }
//     $vcf_line[] = implode(":", array_keys($format_values)); //FORMAT
//     $vcf_line[] = implode(":", array_values($format_values)); //SAMPLE

//     //write assembled line
//     fwrite($fh, implode("\t", $vcf_line)."\n");  

// }

?>