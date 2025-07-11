<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME']) . "/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("create_methyl_plot", "Create methylation plots.");

$parser->addString("folder", "Analysis data folder.", false);
$parser->addString("name", "Processed sample name.", false);
$parser->addString("out", "Output table (TSV format).", false);
$parser->addString("regions", "The regions/highlights to analyze/plot. ", false);

// optional
$parser->addInfile("local_bam", "Optional (local) BAM used for analysis instead of BAM/CRAM file in folder.", true);
$parser->addFlag("skip_plot", "Disable methylartist plotting of imprinting sites.");
$parser->addFlag("skip_align_plot", "Do not create alginment plot in methylartist (Useful for large plots).");
$parser->addString("build", "The genome build to use. ", true, "GRCh38");
$parser->addInt("threads", "The maximum number of threads to use.", true, 4);

extract($parser->parse($argv));

$ref_genome = genome_fasta($build);
exec2("mkdir -p {$folder}/methylartist");
$gtf = get_path("data_folder")."/dbs/Ensembl/Homo_sapiens.GRCh38.112.gtf.gz";
$cram_input = false;

if ($local_bam != "")
{
    $bam = $local_bam;
    if (!file_exists($bam)) trigger_error("Local BAM file not found!", E_USER_ERROR);
}
else
{
    $bam = $folder."/".$name.".cram";
    if (!file_exists($bam))
    {
        $bam = $folder."/".$name.".bam";
        if (!file_exists($bam)) trigger_error("BAM/CRAM file not found!", E_USER_ERROR);
    }
    else
    {
        $cram_input = true;
    }
}

//create temp folder for bam
if ($cram_input) $bam_tmp_folder = $parser->tempFolder("bam");

$regions_table = Matrix::fromTSV($regions);

//new data
$methyl_all_avg = [];
$methyl_all_std = [];
$methyl_hp1_avg = [];
$methyl_hp1_std = [];
$methyl_hp2_avg = [];
$methyl_hp2_std = [];
$methyl_nohp_avg = [];
$methyl_nohp_std = [];

$jobs_plotting = array();
for($r=0; $r<$regions_table->rows(); ++$r)
{
    $row = $regions_table->getRow($r);
    if ($cram_input)
    {
        //create tmp bam file for modkit if input is cram
        $tmp_bam = "{$bam_tmp_folder}/{$name}.bam";
        $parser->execApptainer("samtools", "samtools view", "-o {$tmp_bam} -T {$ref_genome} -u --use-index -@ {$threads} {$bam} {$row[3]}:{$row[6]}-{$row[7]}", [$ref_genome, $bam]);
        $parser->indexBam($tmp_bam, $threads);
    }
    else
    {
        //use bam directly
        $tmp_bam = $bam;
    }

    if (!$skip_plot)
    {
        $args = [
            "locus",
            "--bams", $tmp_bam,
            "--interval", "{$row[3]}:{$row[4]}-{$row[5]}",
            "--highlight", "{$row[6]}-{$row[7]}",
            "--ref", $ref_genome,
            "--motif", "CG",
            "--phased",
            "--color_by_hp",
            "--labelgenes",
            "--genes", $row[2],
            "--ignore_ps",
            "--primary_only",
            "--mods", "m",
            "--outfile", "{$folder}/methylartist/{$name}_{$row[0]}.png",
            "--gtf", $gtf
        ];
        //optional parameters
        if ($skip_align_plot) $args[] = "--skip_align_plot";
        
        $in_files = array($tmp_bam, $ref_genome, $gtf);
        $out_files = array("{$folder}/methylartist/");
        $jobs_plotting[] = array("Plotting_".$row[0], $parser->execApptainer("methylartist", "methylartist", implode(" ", $args), $in_files, $out_files, true));
        
    }

    // calculate average methylation
    $pileup_out = $parser->tempFolder();
    $args_modkit = [
        "pileup",
        "--ref", $ref_genome,
        "--cpg",
        "--region", "{$row[3]}:{$row[6]}-{$row[7]}",
        "--combine-strands",
        "--ignore", "h",
        "--partition-tag", "HP",
        "--prefix", "haplotyped",
        "--threads", $threads,
        $tmp_bam,
        $pileup_out . "/",
    ];
    $parser->execApptainer("modkit", "modkit", implode(" ", $args_modkit), [$ref_genome, $tmp_bam], [$pileup_out]);
    
    $haplotypes = [1, 2, "ungrouped"];
    $modified = [];
    $aggregated = [];
    foreach ($haplotypes as $hp)
    {
        $modified[$hp] = [];
        $pileup = $pileup_out . "/haplotyped_{$hp}.bed";
        if (file_exists($pileup))
        {
            $handle = fopen2($pileup, "r");
            while(!feof($handle))
            {
                $line = trim(fgets($handle));
                if(empty($line)) continue;
                $parts = explode("\t", $line);
                if ($parts[3] !== "m") continue;
                $percent_modified = (float) $parts[10];
                $modified[$hp][] = $percent_modified;
            }
        }
        
        if (count($modified[$hp]) > 0) $aggregated[$hp] = [ mean($modified[$hp]), stdev($modified[$hp]) ];
        else $aggregated[$hp] = [ NAN, NAN];
    }

    // calculate average depth
    $haplotypes = [1, 2];
    $hap_bams = [];
    foreach ($haplotypes as $hp)
    {
        $hap_bams[$hp] = $parser->tempFile(".bam", "hap");
        $args_samtools_view = [
            "-o", $hap_bams[$hp],
            "-u",
            "--use-index",
            "--tag", "HP:{$hp}",
            "-F", "SECONDARY,SUPPLEMENTARY",
            "-@", $threads,
            $tmp_bam,
            "{$row[3]}:{$row[6]}-{$row[7]}",
			"-T {$ref_genome}"
        ];       
        $parser->execApptainer("samtools", "samtools view", implode(" ", $args_samtools_view), [$ref_genome, $tmp_bam]);
        $parser->indexBam($hap_bams[$hp], $threads);
    }


    $target_bed = $parser->tempFile(".bed", "hap");
    $coord_start = $row[6]-1;
    file_put_contents($target_bed, "{$row[3]}\t{$coord_start}\t{$row[7]}");
    $result = $parser->execApptainer("ngs-bits", "BedCoverage", "-threads {$threads} -random_access -decimals 6 -bam {$hap_bams[1]} {$hap_bams[2]} {$tmp_bam} -in {$target_bed} -clear -ref {$ref_genome}", [$tmp_bam, $ref_genome]);
    $result_parts = explode("\t", $result[0][1]);
    $average_coverage = array_map("floatval", array_slice($result_parts, 3, 3));

    $all_data = array_merge($modified[1], $modified[2], $modified["ungrouped"]);
    $aggregated["all"] = ["nan", "nan"];
	if (count($all_data)>0)
	{
		$aggregated["all"] = [ number_format(mean($all_data),2), number_format(stdev($all_data), 2) ];
	}
    $methyl_all_avg[] = $aggregated["all"][0];
    $methyl_all_std[] = $aggregated["all"][1];
    $methyl_hp1_avg[] = number_format($aggregated[1][0], 2);
    $methyl_hp1_std[] = number_format($aggregated[1][1], 2);
    $methyl_hp2_avg[] = number_format($aggregated[2][0], 2);
    $methyl_hp2_std[] = number_format($aggregated[2][1], 2);
    $methyl_nohp_avg[] = number_format($aggregated["ungrouped"][0], 2);
    $methyl_nohp_std[] = number_format($aggregated["ungrouped"][1], 2);
    $coverage_all[] = number_format($average_coverage[2], 2);
    $coverage_hp1[] = number_format($average_coverage[0], 2);
    $coverage_hp2[] = number_format($average_coverage[1], 2);
    $coverage_nohp[] = number_format($average_coverage[2] - $average_coverage[1] - $average_coverage[0], 2);
}

$regions_table->addCol($methyl_all_avg, "meth_all_avg", "average methylation");
$regions_table->addCol($methyl_all_std, "meth_all_std", "standard deviation of methylation");
$regions_table->addCol($methyl_hp1_avg, "meth_hp1_avg", "average methylation in haplotype 1");
$regions_table->addCol($methyl_hp1_std, "meth_hp1_std", "standard deviation of methylation in haplotype 1");
$regions_table->addCol($methyl_hp2_avg, "meth_hp2_avg", "average methylation in haplotype 2");
$regions_table->addCol($methyl_hp2_std, "meth_hp2_std", "standard deviation of methylation in haplotype 2");
$regions_table->addCol($methyl_nohp_avg, "meth_nohp_avg", "average methylation in non-haplotyped");
$regions_table->addCol($methyl_nohp_std, "meth_nohp_std", "standard deviation of methylation in non-haplotyped");
$regions_table->addCol($coverage_all, "cov_all", "average coverage");
$regions_table->addCol($coverage_hp1, "cov_hp1", "average coverage in haplotype 1");
$regions_table->addCol($coverage_hp2, "cov_hp2", "average coverage in haplotype 2");
$regions_table->addCol($coverage_nohp, "cov_nohp", "average coverage in non-haplotyped");

$regions_table->toTSV($out);

# run plots in parallel (reverse array since plotting of REs takes longer)
if (count($jobs_plotting) > 0) $parser->execParallel(array_reverse($jobs_plotting), $threads); 



?>