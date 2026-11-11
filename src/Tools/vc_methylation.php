<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME']) . "/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("vc_methylation", "Performs methylation analysis on a given set of imprinting sites and creates plots.");

$parser->addString("folder", "Analysis data folder.", false);
$parser->addString("name", "Processed sample name.", false);
$parser->addString("out", "Output table (TSV format).", false);
$parser->addString("regions", "The regions/highlights to analyze/plot. ", false);
$parser->addEnum("gender", "Gender of the sample.", false, ["female", "male", "n/a"]);

// optional
$parser->addInfile("local_bam", "Optional (local) BAM/CRAM used for analysis instead of BAM/CRAM file in folder.", true);
$parser->addFlag("skip_plot", "Disable methylartist plotting of imprinting sites.");
$parser->addFlag("skip_align_plot", "Do not create alignment plot in methylartist (Useful for large plots).");
$parser->addFlag("skip_cohort_annotation", "Skip cohort annotation.");
// $parser->addFlag("export_methylation_data", "Also exports table with all methylation values");
$parser->addInfile("custom_cohort_table", "Custom sample table with path in TSV format (NGSDExportSamples) used as background cohort.", true);
$parser->addString("build", "The genome build to use. ", true, "GRCh38");
$parser->addInt("threads", "The maximum number of threads to use.", true, 4);
$parser->addFlag("keep_plotdata", "Store plotdata file in output folder.");
extract($parser->parse($argv));

function get_cohort(&$db, $cohort_ps_name, $bam, $size=100)
{
    $sample_info = get_processed_sample_info($db, $cohort_ps_name);
    $run_dates = get_run_dates($db);

    //get sample date
    $sample_date = $run_dates[$sample_info["run_name"]];
        
    //get basecall model (from BAM since sample QC may not imported yet)
    $basecall_model = get_basecall_model($bam);
    if (str_contains($basecall_model, "hac")) $basecall_model = "hac";
    elseif (str_contains($basecall_model, "sup")) $basecall_model = "sup";
    else $basecall_model = "";

    //get related samples
    $related_samples = get_related_processed_samples($db, $cohort_ps_name, "parent-child");
    $related_samples = array_merge($related_samples, get_related_processed_samples($db, $cohort_ps_name, "siblings"));
    $related_samples = array_merge($related_samples, get_related_processed_samples($db, $cohort_ps_name, "same sample"));
    $related_samples = array_merge($related_samples, get_related_processed_samples($db, $cohort_ps_name, "same patient"));
    $related_samples = array_merge($related_samples, get_related_processed_samples($db, $cohort_ps_name, "cousins"));
    $related_samples = array_merge($related_samples, get_related_processed_samples($db, $cohort_ps_name, "twins"));
    $related_samples = array_merge($related_samples, get_related_processed_samples($db, $cohort_ps_name, "twins (monozygotic)"));

    //create sample table
    $args = array();
    $tmp_sample_table = temp_file("sample_table.tsv");
    $args[] = "-out {$tmp_sample_table}";
    $args[] = "-run_finished";
    $args[] = "-ps_quality good,medium";
    $args[] = "-system_type lrGS";
    $args[] = "-system ".$sample_info["sys_name_short"];
    $args[] = "-add_qc";
    $args[] = "-add_path SAMPLE_FOLDER";
    execApptainer("ngs-bits", "NGSDExportSamples", implode(" ", $args));

    //filter table
    $sample_table = Matrix::fromTSV($tmp_sample_table);
    $filtered_sample_table = new Matrix();
    $filtered_sample_table->setHeaders($sample_table->getHeaders());
    $i_sample = $sample_table->getColumnIndex("name");
    $i_run_name = $sample_table->getColumnIndex("run_name");
    $i_gender = $sample_table->getColumnIndex("gender");
    $i_basecall_info = $sample_table->getColumnIndex("qc_basecall_info");
    $i_path = $sample_table->getColumnIndex("path");
    $run_date_diff = array();

    for ($row_idx=0; $row_idx < $sample_table->rows(); $row_idx++) 
    { 
        //filter by gender:
        $gender = $sample_table->get($row_idx, $i_gender);
        if (($sample_info["gender"] != "n/a") && ($sample_info["gender"] != $gender)) continue;

        //filter same sample + related
        $sample = $sample_table->get($row_idx, $i_sample);
        if ($sample == $cohort_ps_name) continue;
        if (in_array($sample, $related_samples, true)) continue;

       //filter by basecall model 
       if (($basecall_model != "") && !str_contains($sample_table->get($row_idx, $i_basecall_info), $basecall_model))
       {
            continue;
       }

       //skip if modkit files/phasing track doesn't exist
       if (!file_exists($sample_table->get($row_idx, $i_path)."/{$sample}_modkit_track_hp1.bed.gz")) continue;
       if (!file_exists($sample_table->get($row_idx, $i_path)."/{$sample}_modkit_track_hp2.bed.gz")) continue;
       if (!file_exists($sample_table->get($row_idx, $i_path)."/{$sample}_modkit_track.bed.gz")) continue;
       if (!file_exists($sample_table->get($row_idx, $i_path)."/{$sample}_phasing_track.bed")) continue;

       // sample (line) passed all filters
       $filtered_sample_table->addRow($sample_table->getRow($row_idx));
       
       //add time diff to current run
       $run_date_diff[] = date_diff(date_create($run_dates[$sample_table->get($row_idx, $i_run_name)]), date_create($sample_date), true)->format('%d');
    }

    //add date diff as column
    $filtered_sample_table->addCol($run_date_diff, "run_date_diff");
    $i_date_diff = $filtered_sample_table->getColumnIndex("run_date_diff");

    //sort by date diff
    $filtered_sample_table->sort($i_date_diff, SORT_NUMERIC);

    //slice table
    if ($filtered_sample_table->rows() > $size)
    {
        $filtered_sample_table->resize($size, $filtered_sample_table->cols());
    }

    return $filtered_sample_table;

}

// get end_date of all runs
function get_run_dates(&$db)
{
    $run_dates = array();
    $res = $db->executeQuery("SELECT name, end_date FROM sequencing_run");
    foreach($res as $row)
    {
        $run_dates[$row["name"]] = $row["end_date"];
    }

    return $run_dates;
}

// get haplotype count
function get_haplotype_count($gender, $chr, $start, $end)
{
    if($gender == 'female') return 2;
    if($chr == "chrX")
    {
        //check pseudo-autosomal region
        if (range_overlap($start, $end, 10001, 2781479) || range_overlap($start, $end, 155701383, 156030895)) return 2;
        return 1;
    }
    if($chr == "chrY") 
    {
        //check pseudo-autosomal region
        if (range_overlap($start, $end, 10001, 2781479) || range_overlap($start, $end, 56887903, 57217415)) return 2;
        return 1;
    }
    if($chr == "chrMT") return 1;
    if($chr == "chrM") return 1;

    // all other chromosomes are diploid
    return 2;
}


//init
$ref_genome = genome_fasta($build);
exec2("mkdir -p {$folder}/methylartist");
if ($build != "GRCh38") trigger_error("Only GRCh38 is supported as build!", E_USER_ERROR);
$gtf = get_path("data_folder")."/dbs/Ensembl/Homo_sapiens.GRCh38.115.gtf.gz";

//input files
$modkit_hp1 = "{$folder}/{$name}_modkit_track_hp1.bed.gz";
if (!file_exists($modkit_hp1)) trigger_error("Modkit file '{$modkit_hp1}' not found!", E_USER_ERROR);
$modkit_hp2 = "{$folder}/{$name}_modkit_track_hp2.bed.gz";
if (!file_exists($modkit_hp2)) trigger_error("Modkit file '{$modkit_hp2}' not found!", E_USER_ERROR);
$modkit_all = "{$folder}/{$name}_modkit_track.bed.gz";
if (!file_exists($modkit_all)) trigger_error("Modkit file '{$modkit_all}' not found!", E_USER_ERROR);
//check for phasing_track file
$phasing_track_file = $folder."/".$name."_phasing_track.bed";
if (!file_exists($phasing_track_file))
{
    trigger_error("Phasing track file '".basename($phasing_track_file)."' not found! Cannot evaluate phasing of methylation regions.", E_USER_WARNING);
}

if ($local_bam != "")
{
    $bam = $local_bam;
    if (!file_exists($bam)) trigger_error("Local BAM/CRAM file not found!", E_USER_ERROR);
}
else
{
    $bam = $folder."/".$name.".cram";
    if (!file_exists($bam))
    {
        $bam = $folder."/".$name.".bam";
        if (!file_exists($bam)) trigger_error("BAM/CRAM file not found!", E_USER_ERROR);
    }
}

//Assemble cohort
if (!$skip_cohort_annotation)
{
    if ($custom_cohort_table != "")
    {
        $cohort_table = Matrix::fromTSV($custom_cohort_table);
    }
    else
    {
        if (!db_is_enabled("NGSD")) trigger_error("No connection to NGSD! Cannot generate cohort!", E_USER_ERROR);
        $cohort_size = 100;
        $db = DB::getInstance("NGSD");
        $cohort_table = get_cohort($db, $name, $bam, $cohort_size);
    }

    $cohort_files_hp1 = [];
    $cohort_files_hp2 = [];
    $cohort_files_all = [];
    $cohort_files_phasing_tracks = [];
    //iterate over the cohort
    $i_sample_name = $cohort_table->getColumnIndex("name");
    $i_path = $cohort_table->getColumnIndex("path");
    for ($i=0; $i < $cohort_table->rows(); $i++) 
    { 
        $cohort_ps_name = $cohort_table->get($i, $i_sample_name);
        $cohort_sample_folder = realpath2($cohort_table->get($i, $i_path));

        $cohort_files_hp1[] = "{$cohort_sample_folder}/{$cohort_ps_name}_modkit_track_hp1.bed.gz";
        $cohort_files_hp2[] = "{$cohort_sample_folder}/{$cohort_ps_name}_modkit_track_hp2.bed.gz";
        $cohort_files_all[] = "{$cohort_sample_folder}/{$cohort_ps_name}_modkit_track.bed.gz";
        $cohort_files_phasing_tracks[] = "{$cohort_sample_folder}/{$cohort_ps_name}_phasing_track.bed";
    }
}


//run BedMethylationStats
$infiles = [$modkit_hp1, $modkit_hp2, $modkit_all, $regions];
$args = array();
$args[] = "-in_hp1 {$modkit_hp1}";
$args[] = "-in_hp2 {$modkit_hp2}";
$args[] = "-in_all {$modkit_all}";
$args[] = "-sites {$regions}";
$args[] = "-gender {$gender}";
if (!$skip_cohort_annotation)
{
    $infiles = array_merge($infiles, $cohort_files_hp1, $cohort_files_hp2, $cohort_files_all, $cohort_files_phasing_tracks);
    $args[] = "-cohort_hp1 ".implode(" ", $cohort_files_hp1);
    $args[] = "-cohort_hp2 ".implode(" ", $cohort_files_hp2);
    $args[] = "-cohort_all ".implode(" ", $cohort_files_all);
    $args[] = "-cohort_phasing_tracks ".implode(" ", $cohort_files_phasing_tracks);
}
$args_stats = $args;
$args_stats[] = "-out {$out}";
if (file_exists($phasing_track_file))
{
    $args_stats[] = "-phasing_track {$phasing_track_file}";
    $infiles[] = $phasing_track_file;
}
execApptainer("ngs-bits", "BedMethylationStats", implode(" ", $args_stats), $infiles, [dirname($out)]);



//assemble plot commands
if (!$skip_plot)
{
    //extract plot data
    $tmp_plot_data = $parser->tempFile("plot_data.bed");
    $args_plot_data = $args;
    $args_plot_data[] = "-out {$tmp_plot_data}";
    execApptainer("ngs-bits", "BedExtractMethylationPlotData", implode(" ", $args_plot_data), $infiles);
    if ($keep_plotdata) $parser->copyFile($tmp_plot_data, "{$folder}/methylartist/{$name}_plotdata.tsv");

    $regions_table = Matrix::fromTSV($regions);
    $jobs_plotting = array();
    for($r=0; $r<$regions_table->rows(); ++$r)
    {
        $row = $regions_table->getRow($r);

        //create combined plot for male samples on X/Y
        $only_one_haplotype = (get_haplotype_count($gender, $row[3], $row[6], $row[7]) == 1);
        
        //methylartist plot
        $args = [
            "locus",
            "--bams", $bam,
            "--interval", "{$row[3]}:{$row[4]}-{$row[5]}",
            "--highlight", "{$row[6]}-{$row[7]}",
            "--ref", $ref_genome,
            "--motif", "CG",
            "--color_by_hp",
            "--labelgenes",
            "--genes", $row[2],
            "--ignore_ps",
            "--primary_only",
            "--mods", "m",
            "--outfile", "{$folder}/methylartist/{$name}_{$row[0]}.png",
            "--gtf", $gtf
        ];

        // use phased plots not for chrX/Y on male samples
        if (!$only_one_haplotype) $args[] =  "--phased";

        //optional parameters
        if ($skip_align_plot) $args[] = "--skip_align_plot";
        
        $infiles = array($bam, $ref_genome, $gtf);
        $outfiles = array("{$folder}/methylartist/");
        $jobs_plotting[] = array("Plotting_".$row[0], $parser->execApptainer("methylartist", "methylartist", implode(" ", $args), $infiles, $outfiles, true));

        
        //cohort plot
        $python_script = repository_basedir()."/src/Tools/generate_methylation_cohort_plot.py";

        $args = [];
        $args[] = "--plot_data {$tmp_plot_data}";
        $args[] = "--out {$folder}/methylartist/{$name}_{$row[0]}_cohort.png";
        $args[] = "--site {$row[0]}";
        $args[] = "--sample_name {$name}";
        $args[] = "--chr {$row[3]}";
        $args[] = "--start {$row[4]}";
        $args[] = "--end {$row[5]}";
        $args[] = "--highlight_start {$row[6]}";
        $args[] = "--highlight_end {$row[7]}";
        if ($only_one_haplotype) $args[] =  "--unphased";

        $infiles = [$python_script];
        $outfiles = ["{$folder}/methylartist/"];
        $jobs_plotting[] = array("Plotting_".$row[0]."_cohort", $parser->execApptainer("python", "python3", $python_script." ".implode(" ", $args), $infiles, $outfiles, true));            
    }

    // run plots in parallel (reverse array since plotting of REs takes longer)
    $return = $parser->execParallel(array_reverse($jobs_plotting), $threads, false, false, false, true); 
    //check failed jobs
    foreach ($return["jobs_failed"] as $job_id) 
    {
        $job_info = $return["jobs"][$job_id];

        // ignore failed jobs beacuse of coverage
        if (($job_info["exit_code"] == 1) && (str_contains(implode("\n", $job_info["stderr"]), "insufficient coverage for plot"))) 
        {
            trigger_error("Processing of job {$job_id} failed due to insufficient coverage!", E_USER_WARNING);
            continue;
        }

        trigger_error("Processing of job {$job_id} failed: ".implode("\n\t", $job_info["stderr"]), E_USER_ERROR);
    }

}

?>