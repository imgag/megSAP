<?php
require_once(dirname($_SERVER['SCRIPT_FILENAME']) . "/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("vc_arriba", "Run fusion detection with Arriba.");

$parser->addInfile("bam", "Input BAM file.", false);

$parser->addOutfile("out_fusions", "Fusion report in TSV format.", false);
$parser->addOutfile("out_vcf", "Fusion report in VCF4.3 format.", false);
$parser->addOutfile("out_discarded", "Discarded fusions in TSV format.", true);
$parser->addOutfile("out_pdf", "Fusion report in PDF format.", true);
$parser->addOutfile("out_bam", "Output BAM file with fusion-supporting reads.", true);
$parser->addString("out_pic_dir", "Output directory which contains all fusion pictures as PNGs.", true);
$parser->addInfile("sv", "Optional structural variants from DNA sequencing, in VCF format.", true);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
extract($parser->parse($argv));


if (!in_array($build, ["GRCh37", "GRCh38"]))
{
    trigger_error("Annotation only available for GRCh37/GRCh38!", E_USER_ERROR);
}


//resolve build string used by Arriba
$arriba_str = [
    "GRCh37" => "hg19_hs37d5_GRCh37",
    "GRCh38" => "hg38_GRCh38"
];
$arriba_build = $arriba_str[$build];


//reference files
$genome = genome_fasta($build);
$gtf = get_path("data_folder") . "/dbs/gene_annotations/{$build}.gtf";
$arriba_ref = get_path("arriba") . "/database";


//run Arriba
$arriba_ver = "v2.4.0";
$args = [
    "-x", $bam,
    "-o", $out_fusions,
    "-a", $genome,
    "-g", $gtf,
    "-b", "{$arriba_ref}/blacklist_{$arriba_build}_{$arriba_ver}.tsv.gz",
    "-k", "{$arriba_ref}/known_fusions_{$arriba_build}_{$arriba_ver}.tsv.gz",
    "-t", "{$arriba_ref}/known_fusions_{$arriba_build}_{$arriba_ver}.tsv.gz",
    "-p", "{$arriba_ref}/protein_domains_{$arriba_build}_{$arriba_ver}.gff3",
    "-X",
    "-f", "no_genomic_support,read_through,same_gene,intragenic_exonic"
];

if (isset($sv)) $args[] = "-d {$sv}";
if (isset($out_discarded)) $args[] = "-O {$out_discarded}";


echo "starting arriba!";
//In rare cases throws segmentation fault while writing discarded_fusion file. 
//Allow that crash while catching others. Writing the discarded file is the last thing it does before freeing resources and finishing -> Check it startet writing the discarded file to be safe the regular output file is complete.
exec(get_path("arriba")."/arriba ". implode(" ", $args)." 2>&1", $output, $exit_code);

echo "finished arriba!\n";

$parser->log("arriba stdout + stderr:");
foreach($output as $line)
{
	$parser->log($line);
}

If ($exit_code != 0)
{
	$startet_dsicarded = false;
	foreach($output as $line)
	{
		if (contains($line, "Writing discarded fusions to file"))
		{
			$startet_discarded = true;
			break;
		}
	}
	
	if ($startet_discarded)
	{
		trigger_error("Arriba exit code was not 0, but the output fusions are complete. The discarded fusions are probably not complete and will be removed.", E_USER_NOTICE);
		if (file_exists($out_discarded))
		{
			exec2("rm $out_discarded");
		}
	}
	else
	{
		trigger_error("Arriba exit code was not 0 and the fusions file is probably not done / complete.", E_USER_ERROR);
	}
}


if (isset($out_vcf)) {
    $args = [
        $genome,
        $out_fusions,
        $out_vcf
    ];
    $parser->exec(get_path("arriba") . "/scripts/convert_fusions_to_vcf.sh", implode(" ", $args));
}


//generate plot in PDF format
//limit to top ~20 fusions
$top_fusions = $parser->tempFile("_top_fusions.txt");
$parser->exec("head", "-n 21 {$out_fusions} > {$top_fusions}");
if (isset($out_pdf)) {
    $plot_args = [
        "--annotation={$gtf}",
        "--fusions={$top_fusions}",
        "--output={$out_pdf}",
        "--alignments={$bam}",
        "--cytobands={$arriba_ref}/cytobands_{$arriba_build}_{$arriba_ver}.tsv",
        "--proteinDomains={$arriba_ref}/protein_domains_{$arriba_build}_{$arriba_ver}.gff3",
        "--minConfidenceForCircosPlot=none"
    ];
    $parser->exec(get_path("rscript")." ". get_path("arriba") . "/draw_fusions.R", implode(" ", $plot_args));
	
	if (isset($out_pic_dir))
	{
		if(!file_exists($out_pic_dir)) mkdir($out_pic_dir);
		$parser->exec("gs", "-sDEVICE=png16m -dTextAlphaBits=4 -r300 -o {$out_pic_dir}/%04d.png $out_pdf");
	}
}




//extract fusion-supporting reads in extra BAM file
if (isset($out_bam)) {
    $fusions = Matrix::fromTSV($out_fusions);
    $reads_rows = $fusions->getCol($fusions->getColumnIndex("read_identifiers"));
    //fix read id concatenation, UMI reads use , in read id
    if ((count($reads_rows) !== 0) && (contains($reads_rows[0], ",,")))
    {
        $reads = explode("|", implode("|", str_replace(",,", ",|", $reads_rows)), -1);
    }
    else
    {
        $reads = explode(",", implode(",", $reads_rows), -1);
    }
    if (count($reads) !== 0)
    {
        $read_ids = $parser->tempFile("_readids.txt");
        file_put_contents($read_ids, implode("\n", $reads));
        $bam_header = $parser->tempFile("_header.txt");
        $bam_records = $parser->tempFile("_records.sam");

        $parser->exec(get_path("samtools"), "view -H {$bam} > {$bam_header}");
        $parser->execPipeline([
            [get_path("samtools"), "view -O SAM {$bam}"],
            ["grep", "-F -f {$read_ids} > {$bam_records}"]
        ], "filter BAM");
        $parser->execPipeline([
            ["cat", "{$bam_header} ${bam_records}"],
            [get_path("samtools"), "view -o {$out_bam}"]
        ], "write BAM");
        $parser->indexBam($out_bam, 1);
    }
}
