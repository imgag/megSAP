<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME']) . "/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("create_methyl_plot", "Create methylation plots.");

$parser->addString("folder", "Analysis data folder.", false);
$parser->addString("name", "Processed sample name.", false);
$parser->addString("out", "Output table (TSV format).", false);
$parser->addString("regions", "The regions/highlights to analyze/plot. ", false);

// optional
$parser->addInfile("local_bam", "Optional (local) BAM/CRAM used for analysis instead of BAM/CRAM file in folder.", true);
$parser->addFlag("skip_plot", "Disable methylartist plotting of imprinting sites.");
$parser->addFlag("skip_align_plot", "Do not create alginment plot in methylartist (Useful for large plots).");
$parser->addFlag("skip_cohort_annotation", "Skip cohort annotation.");
$parser->addFlag("export_methylation_data", "Also exports table with all methylation values");
$parser->addInfile("custom_cohort_table", "Custom sample table with path in TSV format (NGSDExportSamples) used as background cohort.", true);
$parser->addString("build", "The genome build to use. ", true, "GRCh38");
$parser->addInt("threads", "The maximum number of threads to use.", true, 4);
$parser->addFlag("test", "Run in test-mode with hard-coded cohort.");

extract($parser->parse($argv));


//get (tmp) BAM file, if input is CRAM
function get_bam($bam, $ps_name, $chr, $start, $end, $ref_genome, $threads)
{
    global $parser; //get ToolBase
    global $folder;

    // if BAM is provided: use BAM directly
    if (ends_with($bam, ".bam")) return $bam;

    //else: create temporary file over the region
    $tmp_bam = $parser->tempFolder("tmp_bam")."/{$ps_name}.bam";
    //add 100kb in each direction for better modkit normalization 

    $start = max(1, intval($start) - 100000);
    $end = intval($end) + 100000;
    $parser->execApptainer("samtools", "samtools view", "-b -o {$tmp_bam} -T {$ref_genome} -u --use-index -@ {$threads} {$bam} {$chr}:{$start}-{$end}", [$ref_genome, $bam]);
    $parser->indexBam($tmp_bam, $threads);

    //debug
    $file_size = filesize($tmp_bam) / (1024*1024);
    trigger_error("Temp BAM filesize: \t{$file_size}", E_USER_NOTICE);

    return $tmp_bam;
}

// run modkit for a region and returns (tmp) output folder with the results 
function extract_methylation($bam, $chr, $start, $end, $highlight_start, $highlight_end, $name, $ref_genome, $threads, $unphased=false)
{
    global $parser; //get ToolBase
    // create methyl track for cohort plots
    $tmp_out = $parser->tempFolder()."/";
    $prefix = "{$name}";
    $prefix_highlight = "{$name}_highlight";
    $args_modkit = [
        "pileup",
        $bam,
        (($unphased)? "{$tmp_out}{$prefix}_all.bed": $tmp_out), //for phased output provide folder, otherwise output bed file
        "--ref", $ref_genome,
        "--cpg",
        "--region", "{$chr}:{$start}-{$end}",
        "--combine-strands",
        "--ignore", "h",
        "--threads", $threads,
    ];
    if (!$unphased)
    {
        $args_modkit[] = "--partition-tag HP";
        $args_modkit[] = "--prefix {$prefix}";
    } 

    list($stdout, $stderr, $ec) = $parser->execApptainer("modkit", "modkit", implode(" ", $args_modkit), [$ref_genome, $bam], [$tmp_out], false, true, false, true);
    // workaround to allow 0 reads in BAM region (e.g. in test cases)
    if (($ec != 0) && (end($stderr) != "> Error! zero reads found in bam index")) 
    {
        $parser->toStderr($stdout);
		$parser->toStderr($stderr);
		trigger_error("Call of external tool 'modkit' returned error code '$ec'.", E_USER_ERROR);
    }

    //slice BED to target region
    $haplotypes = [1, 2, "ungrouped"];
    if ($unphased) $haplotypes = ["all"];
    foreach ($haplotypes as $hp)
    {
        $input_bed = "{$tmp_out}{$prefix}_{$hp}.bed";
        //create empty file in case input doesn't exists
        if (!file_exists($input_bed)) touch($input_bed);
        $output_bed = "{$tmp_out}/{$prefix_highlight}_{$hp}.bed";
        $pipeline = [
            ["", "echo '{$chr}\t{$highlight_start}\t{$highlight_end}'"],
            ["", $parser->execApptainer("ngs-bits", "BedIntersect", "-mode in2 -in2 {$input_bed} -out {$output_bed}", [$input_bed], [$tmp_out], true)]
        ];
        $parser->execPipeline($pipeline, "extract highlight region", false);   
    }

    return $tmp_out;
}

// get avg methylation from file(s)
function get_avg_methylation($file_path_prefix, $sort_haplotype=true)
{
    $methylation = array();
    $methylation["all"] = array();
    $methylation["all"]["N_mod"] = array();
    $methylation["all"]["N_valid_cov"] = array();
    $methylation["all"]["fraction_modified"] = array();
    $methylation["hp1_hp2_switched"] = false;
    foreach ([1, 2, "ungrouped"] as $hp)
    {
        $methylation[$hp]["N_mod"] = array();
        $methylation[$hp]["N_valid_cov"] = array();
        $methylation[$hp]["fraction_modified"] = array();
        $buffer = file("{$file_path_prefix}_{$hp}.bed", FILE_IGNORE_NEW_LINES + FILE_SKIP_EMPTY_LINES);
        foreach ($buffer as $line)
        {
            $parts = explode("\t", $line);
            if ($parts[3] !== "m") continue; //skip other methylation

            $n_mod = intval($parts[11]);
            $pos = "{$parts[0]}:{$parts[1]}-{$parts[2]}";
            $methylation[$hp]["N_mod"][$pos] = $n_mod;
            $n_valid_cov = intval($parts[9]);
            $methylation[$hp]["N_valid_cov"][$pos] = $n_valid_cov;
            $methylation[$hp]["fraction_modified"][$pos] = floatval($parts[10]);

            //aggregate all columns
            if (isset($methylation["all"]["N_mod"][$pos]))
            {
                $methylation["all"]["N_mod"][$pos] += $n_mod;
                $methylation["all"]["N_valid_cov"][$pos] += $n_valid_cov;
            } 
            else
            {
                $methylation["all"]["N_mod"][$pos] = $n_mod;
                $methylation["all"]["N_valid_cov"][$pos] = $n_valid_cov;
            } 
        }

        // if ($hp === "ungrouped") $hp_str = "nohp";
        // else $hp_str = "hp{$hp}";
        $methylation[$hp]["avg"] = NAN;
        $methylation[$hp]["std"] = NAN;
        $methylation[$hp]["cov"] = NAN;

        if (count($methylation[$hp]["N_mod"]) > 0)
        {
            $methylation[$hp]["avg"] = mean($methylation[$hp]["fraction_modified"]);
            $methylation[$hp]["std"] = stdev($methylation[$hp]["fraction_modified"]);
            $methylation[$hp]["cov"] = mean($methylation[$hp]["N_valid_cov"]);
        }
    }

    //calculate combined stats
    foreach ($methylation["all"]["N_mod"] as $pos => $n_mod) 
    {
        $n_valid_cov = $methylation["all"]["N_valid_cov"][$pos];
        $methylation["all"]["fraction_modified"][$pos] = floatval($n_mod / $n_valid_cov) * 100.0;
    }
    $methylation["all"]["avg"] = NAN;
    $methylation["all"]["std"] = NAN;
    $methylation["all"]["cov"] = NAN;
    if (count($methylation["all"]["fraction_modified"]) > 0)
    {
        $methylation["all"]["avg"] = mean($methylation["all"]["fraction_modified"]);
        $methylation["all"]["std"] = stdev($methylation["all"]["fraction_modified"]);
        $methylation["all"]["cov"] = mean($methylation["all"]["N_valid_cov"]);
    }
    

    if ($sort_haplotype)
    {
        //move higher HP to front
        if ($methylation[1]["avg"] < $methylation[2]["avg"])
        {
            $tmp = $methylation[1];
            $methylation[1] = $methylation[2];
            $methylation[2] = $tmp;
            $methylation["hp1_hp2_switched"] = true;
        }

    }

    /*

    //TODO: restructure and add depth
    //TODO: sort hp by mean methylation
    $modified = [];
    $aggregated = [];
    foreach ([1, 2, "ungrouped"] as $hp)
    {
        $modified[$hp] = [];
        $buffer = file("{$file_path_prefix}_{$hp}.bed", FILE_IGNORE_NEW_LINES + FILE_SKIP_EMPTY_LINES);
        foreach ($buffer as $line)
        {
            $parts = explode("\t", $line);
            if ($parts[3] !== "m") continue;
            $percent_modified = (float) $parts[10];
            $modified[$hp][] = $percent_modified;
        }        
        if (count($modified[$hp]) > 0) $aggregated[$hp] = [ mean($modified[$hp]), stdev($modified[$hp]) ];
        else $aggregated[$hp] = [NAN, NAN];
    }
    return array($modified, $aggregated);
    */
    return $methylation;
}

//get phasing blocks of region
function get_phasing_info($bed_file, $chr, $start, $end, $gender)
{
    global $parser;

    if (!file_exists($bed_file)) return "NoPhasingInfo";

    $pipeline = [
        ["", "echo '{$chr}\t{$start}\t{$end}'"],
        ["", $parser->execApptainer("ngs-bits", "BedIntersect", "-in2 {$bed_file}", [$bed_file], [], true)]
    ];
    list($stdout, $stderr) = $parser->execPipeline($pipeline, "extract highlight region");

    if (count($stdout) > 1) return "MultiplePhasingBlocks"; // region contains multiple phasing block
    if (count($stdout) < 1) // region not phased
    {
        if (get_haplotype_count($gender, $chr, $start, $end) == 1) return "";
        return "UnphasedRegion";
    }  
    if (trim($stdout[0]) != "{$chr}\t{$start}\t{$end}") return "PartlyUnphasedRegion"; // phasing block doesn't cover whole region
    //else: phasing block covers whole area
    return "";
}


function get_cohort(&$db, $cohort_ps_name, $bam, $size=100)
{
    $sample_info = get_processed_sample_info($db, $cohort_ps_name);
    $run_dates = get_run_dates($db);

    //get sample date
    $sample_date = $run_dates[$sample_info["run_name"]];
        
    //get basecall model (from BAM since sample QC may not imported yet)
    /*
    $qc_id = $db->getValue("SELECT id FROM qc_terms WHERE qcml_id = 'QC:2000149'");
    $basecall_model = $db->getValue("SELECT value FROM `processed_sample_qc` WHERE `qc_terms_id` = {$qc_id} AND processed_sample_id=".$sample_info["ps_id"]);
    */
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
    $args[] = "-no_bad_samples";
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
    $i_quality = $sample_table->getColumnIndex("processed_sample_quality");
    $i_gender = $sample_table->getColumnIndex("gender");
    $i_basecall_info = $sample_table->getColumnIndex("qc_basecall_info");
    $i_path = $sample_table->getColumnIndex("path");
    $run_date_diff = array();

    for ($row_idx=0; $row_idx < $sample_table->rows(); $row_idx++) 
    { 
        //filter by quality (only use 'good' and 'medium')
        $quality = $sample_table->get($row_idx, $i_quality);
        if (!($quality == "good") && !($quality == "medium")) continue;

        //filter by gender:
        $gender = $sample_table->get($row_idx, $i_gender);
        if (($sample_info["gender"] != "n/a") && ($sample_info["gender"] != $gender)) continue;



        //filter same sample + related
        $sample = $sample_table->get($row_idx, $i_sample);
        if ($sample == $cohort_ps_name) continue;
        if (array_search($sample, $related_samples)) continue;

       //filter by basecall model 
       if (($basecall_model != "") && !str_contains($sample_table->get($row_idx, $i_basecall_info), $basecall_model))
       {
            continue;
       }

       //skip if methylation file doesn't exist
       if (!file_exists($sample_table->get($row_idx, $i_path)."/{$sample}_var_methylation.tsv")) continue;

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

    //TODO: slice table
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


$gender = "n/a";
if ($test)
{
    // init test run
    trigger_error("Running in test mode!", E_USER_WARNING);
    $sample_info["sys_name_short"] = "ProcessingSystemForTesting";
    $sample_info["ps_folder"] = $folder;
    $gender = "male";
    if ($custom_cohort_table == "") trigger_error("Custom cohort table required for testing!", E_USER_ERROR);
    $cohort_folder = repository_basedir()."/test/data_out/tool_test_create_methyl_plot/cohort/ProcessingSystemForTesting";
}
else
{
    // init NGSD
    if (db_is_enabled("NGSD"))
    {
        $db = DB::getInstance("NGSD");
        $sample_info = get_processed_sample_info($db, $name);
        $gender = $sample_info["gender"];
        $cohort_folder = get_path("methylation_cohorts")."/".$sample_info["sys_name_short"];
    } 
    else trigger_error("No NGSD connection! Cannot perform cohort analysis.", E_USER_WARNING);
}

//check sample folder
$skip_copy = false;
if (!isset($db))
{
    $skip_copy = true;
    trigger_error("NGSD is disabled! Skip copying/replacing of cohort methylation files.", E_USER_WARNING);
}
else if (realpath($folder) != realpath($sample_info["ps_folder"]))
{
    $skip_copy = true;
    trigger_error("Sample is not in correct folder according to NGSD! Skip copying/replacing of cohort methylation files.", E_USER_WARNING);
}


$ref_genome = genome_fasta($build);
exec2("mkdir -p {$folder}/methylartist");
$gtf = get_path("data_folder")."/dbs/Ensembl/Homo_sapiens.GRCh38.115.gtf.gz";

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
$coverage_all = [];
$coverage_hp1 = [];
$coverage_hp2 = [];
$coverage_nohp = []; 

//filter column
$filter_column = [];

//check for phasing_track file
$phasing_track_file = $folder."/".$name."_phasing_track.bed";
if (!file_exists($phasing_track_file))
{
    trigger_error("Phasing track file '".basename($phasing_track_file)."' not found! Cannot evaluate phasing of methylation regions.", E_USER_WARNING);
}


//cohort methylation table
$raw_methylation_data = array();
$cohort_plot_data = array();


$jobs_plotting = array();
$jobs_plotting_cohort = array();
for($r=0; $r<$regions_table->rows(); ++$r)
{
    $row = $regions_table->getRow($r);

    //create combined plot for male samples on X/Y
    $only_one_haplotype = (get_haplotype_count($gender, $row[3], $row[6], $row[7]) == 1);

    $tmp_bam = get_bam($bam, $name, $row[3], $row[6], $row[7], $ref_genome, $threads);

    if (!$skip_plot)
    {
        if (get_read_count($tmp_bam) < 1) 
        {
            trigger_error("No reads in BAM, skipping methylartist plot for {$row[0]}!", E_USER_WARNING);
        }
        else
        {
            $args = [
                "locus",
                "--bams", $tmp_bam,
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
            
            $in_files = array($tmp_bam, $ref_genome, $gtf);
            $out_files = array("{$folder}/methylartist/");
            $jobs_plotting[] = array("Plotting_".$row[0], $parser->execApptainer("methylartist", "methylartist", implode(" ", $args), $in_files, $out_files, true));    
        }
        
    }

    $prefix = "{$name}";
    $prefix_highlight = "{$name}_highlight";
    $modkit_temp_folder = extract_methylation($tmp_bam, $row[3], $row[4], $row[5], $row[6], $row[7], $name, $ref_genome, $threads);
    $haplotypes = [1, 2, "ungrouped"];
    if ($only_one_haplotype)
    {
        $modkit_temp_folder_all = extract_methylation($tmp_bam, $row[3], $row[4], $row[5], $row[6], $row[7], $name, $ref_genome, $threads, true);
    }
    
    if (isset($db))
    {
        //ignore if sample is not in main path
        if (!$skip_copy)
        {
            //copy to cohort folder
            $target_path = "{$cohort_folder}/{$row[0]}_{$row[3]}_{$row[4]}-{$row[5]}/";
            if (!file_exists($target_path)) mkdir($target_path, 0777, true);
            $target_path_highlight = "{$cohort_folder}/{$row[0]}_highlight_{$row[3]}_{$row[6]}-{$row[7]}/";
            if (!file_exists($target_path_highlight)) mkdir($target_path_highlight, 0777, true);

            foreach ($haplotypes as $hp)
            {
                $parser->copyFile("{$modkit_temp_folder}/{$prefix}_{$hp}.bed", "{$target_path}/{$prefix}_{$hp}.bed");
                $parser->copyFile("{$modkit_temp_folder}/{$prefix_highlight}_{$hp}.bed", "{$target_path_highlight}/{$prefix}_{$hp}.bed");
            }
            if ($only_one_haplotype) 
            {
                $parser->copyFile("{$modkit_temp_folder_all}/{$prefix_highlight}_all.bed", "{$target_path_highlight}/{$prefix}_all.bed");
            }
        }
        
    }

    //check phasing
    $filter_column[] = get_phasing_info($phasing_track_file, $row[3], $row[6], $row[7], $gender);

    $methylation = get_avg_methylation("{$modkit_temp_folder}/{$prefix_highlight}");

    $methyl_all_avg[] = number_format($methylation["all"]["avg"], 3);
    $methyl_all_std[] = number_format($methylation["all"]["std"], 3);
    $methyl_hp1_avg[] = number_format($methylation[1]["avg"], 3);
    $methyl_hp1_std[] = number_format($methylation[1]["std"], 3);
    $methyl_hp2_avg[] = number_format($methylation[2]["avg"], 3);
    $methyl_hp2_std[] = number_format($methylation[2]["std"], 3);
    $methyl_nohp_avg[] = number_format($methylation["ungrouped"]["avg"], 3);
    $methyl_nohp_std[] = number_format($methylation["ungrouped"]["std"], 3);
    $coverage_all[] = number_format($methylation["all"]["cov"], 3);
    $coverage_hp1[] = number_format($methylation[1]["cov"], 3);
    $coverage_hp2[] = number_format($methylation[2]["cov"], 3);
    $coverage_nohp[] = number_format($methylation["ungrouped"]["cov"], 3); 

    //store raw data
    $raw_methylation_data[$row[0]] = array();

    foreach ($methylation["all"]["fraction_modified"] as $pos => $all) 
    {
        if (isset($methylation[1]["fraction_modified"][$pos])) $raw_methylation_data[$row[0]][$pos][$name."_hp1"] = $methylation[1]["fraction_modified"][$pos];
        if (isset($methylation[2]["fraction_modified"][$pos])) $raw_methylation_data[$row[0]][$pos][$name."_hp2"] = $methylation[2]["fraction_modified"][$pos];
        $raw_methylation_data[$row[0]][$pos][$name."_all"] = $all;
    }

    //get cohort plot data
    $cohort_plot_data[$row[0]] = array();
    $cohort_plot_data[$row[0]]["hp1_cohort_files"] = array();
    $cohort_plot_data[$row[0]]["hp2_cohort_files"] = array();
    $cohort_plot_data[$row[0]]["hp1_cohort_files_fallback"] = array();
    $cohort_plot_data[$row[0]]["hp2_cohort_files_fallback"] = array();
    $cohort_plot_data[$row[0]]["out"] = "{$folder}/methylartist/{$name}_{$row[0]}_cohort.png";
    $cohort_plot_data[$row[0]]["site"] = $row[0];
    $cohort_plot_data[$row[0]]["sample_name"] = $name;
    $cohort_plot_data[$row[0]]["highlight_start"] = $row[6];
    $cohort_plot_data[$row[0]]["highlight_end"] = $row[7];
    $cohort_plot_data[$row[0]]["unphased"] = $only_one_haplotype;


    if ($only_one_haplotype)
    {
        $cohort_plot_data[$row[0]]["hp1_file"] = "{$modkit_temp_folder_all}/{$prefix}_all.bed";
    }
    else
    {
        if ($methylation["hp1_hp2_switched"])
        {
            $cohort_plot_data[$row[0]]["hp1_file"] = "{$modkit_temp_folder}/{$prefix}_2.bed";
            $cohort_plot_data[$row[0]]["hp2_file"] = "{$modkit_temp_folder}/{$prefix}_1.bed";
        }
        else
        {
            $cohort_plot_data[$row[0]]["hp1_file"] = "{$modkit_temp_folder}/{$prefix}_1.bed";
            $cohort_plot_data[$row[0]]["hp2_file"] = "{$modkit_temp_folder}/{$prefix}_2.bed";
        }
    }
}


$regions_table->addCol($filter_column, "filter", "");

$regions_table->addCol($methyl_hp1_avg, "meth_hp1_avg", "average methylation in haplotype 1 (higher mean methylation)");
$regions_table->addCol($methyl_hp1_std, "meth_hp1_std", "standard deviation of methylation in haplotype 1 (higher mean methylation)");
$regions_table->addCol($coverage_hp1, "cov_hp1", "average CpG coverage in haplotype 1 (higher mean methylation)");
$regions_table->addCol($methyl_hp2_avg, "meth_hp2_avg", "average methylation in haplotype 2 (lower mean methylation)");
$regions_table->addCol($methyl_hp2_std, "meth_hp2_std", "standard deviation of methylation in haplotype 2 (lower mean methylation)");
$regions_table->addCol($coverage_hp2, "cov_hp2", "average CpG coverage in haplotype 2 (lower mean methylation)");
$regions_table->addCol($methyl_nohp_avg, "meth_nohp_avg", "average methylation in non-haplotyped");
$regions_table->addCol($methyl_nohp_std, "meth_nohp_std", "standard deviation of methylation in non-haplotyped");
$regions_table->addCol($coverage_nohp, "cov_nohp", "average CpG coverage in non-haplotyped");
$regions_table->addCol($methyl_all_avg, "meth_all_avg", "average methylation");
$regions_table->addCol($methyl_all_std, "meth_all_std", "standard deviation of methylation");
$regions_table->addCol($coverage_all, "cov_all", "average CpG coverage");


$regions_table->toTSV($out);

//annotate with cohort samples
if (!$skip_cohort_annotation && (db_is_enabled("NGSD") || $test))
{
    if ($custom_cohort_table != "")
    {
        $cohort_table = Matrix::fromTSV($custom_cohort_table);
    }
    else
    {
        if (!isset($db)) trigger_error("No connection to NGSD! Cannot generate cohort!", E_USER_ERROR);
        $cohort_size = 100;
        $cohort_table = get_cohort($db, $name, $bam, $cohort_size);
    }

    //parse current table 
    $meth_table = Matrix::fromTSV($out);
    $i_identifier = $meth_table->getColumnIndex("identifier");
    $i_chr = $meth_table->getColumnIndex("gene chr");
    $i_highlight_start = $meth_table->getColumnIndex("highlight start");
    $i_highlight_end = $meth_table->getColumnIndex("highlight end");
    $i_gene_start = $meth_table->getColumnIndex("gene start");
    $i_gene_end = $meth_table->getColumnIndex("gene end");
    $i_meth_hp1_avg = $meth_table->getColumnIndex("meth_hp1_avg");
    $i_meth_hp2_avg = $meth_table->getColumnIndex("meth_hp2_avg");
    $i_meth_all_avg = $meth_table->getColumnIndex("meth_all_avg");
    
    $highlight_regions = array();
    $cohort_values = array();
    $sample_methylation = array();
    

    for($r=0; $r<$meth_table->rows(); ++$r)
    {
        $identifier = $meth_table->get($r, $i_identifier);
        $highlight_regions[$identifier] = [
            $meth_table->get($r, $i_chr), 
            $meth_table->get($r, $i_highlight_start), 
            $meth_table->get($r, $i_highlight_end), 
            $meth_table->get($r, $i_gene_start), 
            $meth_table->get($r, $i_gene_end)
        ];
        $cohort_values[$identifier]["hp1"] = array();
        $cohort_values[$identifier]["hp2"] = array();
        $cohort_values[$identifier]["all"] = array();
        $sample_methylation[$identifier]["hp1"] = $meth_table->get($r, $i_meth_hp1_avg);
        $sample_methylation[$identifier]["hp2"] = $meth_table->get($r, $i_meth_hp2_avg); 
        $sample_methylation[$identifier]["all"] = $meth_table->get($r, $i_meth_all_avg); 

    }

    //iterate over the cohort
    $i_sample_name = $cohort_table->getColumnIndex("name");
    $i_path = $cohort_table->getColumnIndex("path");
    $cohort_sample_list = array();
    $excluded_samples = array();
    for ($i=0; $i < $cohort_table->rows(); $i++) 
    { 
        $cohort_ps_name = $cohort_table->get($i, $i_sample_name);
        $cohort_sample_folder = $cohort_table->get($i, $i_path);
        //modify cohort folder path for testing:
        if ($test)$cohort_sample_folder = repository_basedir()."/test/data/create_methyl_plot/".$cohort_sample_folder;

        //get files
        //BAM/CRAM
        $cohort_bam_file = "{$cohort_sample_folder}/{$cohort_ps_name}.bam";
        if (!file_exists($cohort_bam_file)) $cohort_bam_file = "{$cohort_sample_folder}/{$cohort_ps_name}.cram";
        if (!file_exists($cohort_bam_file)) trigger_error("No BAM/CRAM file found in sample folder '{$cohort_sample_folder}'!", E_USER_ERROR);
        //phasing track
        $cohort_phasing_track_file = "{$cohort_sample_folder}/{$cohort_ps_name}_phasing_track.bed";
        if (!file_exists($cohort_phasing_track_file))
        {
            trigger_error("Phasing track file for sample {$cohort_ps_name} ({$cohort_phasing_track_file}) not found! Ignoring for cohort.", E_USER_WARNING);
            continue;
        }
        
        $cohort_sample_list[] = $cohort_ps_name;

        foreach ($highlight_regions as $identifier => list($chr, $start, $end, $gene_start, $gene_end)) 
        {
            //check phasing
            if(get_phasing_info($cohort_phasing_track_file, $chr, $start, $end, $gender) != "")
            {
                if (!isset($excluded_samples[$identifier])) $excluded_samples[$identifier] = array();
                $excluded_samples[$identifier][] = $cohort_ps_name;
                continue;
            }
            $only_one_haplotype = (get_haplotype_count($gender, $chr, $start, $end) == 1);
            
            //check for existing file
            if (isset($cohort_folder))
            {
                $modkit_file_highlight_prefix = "{$cohort_folder}/{$identifier}_highlight_{$chr}_{$start}-{$end}/{$cohort_ps_name}";
                $modkit_file_prefix = "{$cohort_folder}/{$identifier}_{$chr}_{$gene_start}-{$gene_end}/{$cohort_ps_name}";
            }
            $tmp_bam = "";
            if (!isset($cohort_folder) || !file_exists($modkit_file_highlight_prefix."_1.bed") || !file_exists($modkit_file_highlight_prefix."_2.bed") || !file_exists($modkit_file_highlight_prefix."_ungrouped.bed") 
                || !file_exists($modkit_file_prefix."_1.bed") || !file_exists($modkit_file_prefix."_2.bed") || !file_exists($modkit_file_prefix."_ungrouped.bed"))
            {
                //re-create modkit files
                $tmp_bam = get_bam($cohort_bam_file, $cohort_ps_name, $chr, $gene_start, $gene_end, $ref_genome, $threads);
                $modkit_temp_folder = extract_methylation($tmp_bam, $chr, $gene_start, $gene_end, $start, $end, $cohort_ps_name, $ref_genome, $threads);

                //copy files to cohort folder
                if (isset($cohort_folder))
                {
                    //create folders if they don't exist:
                    foreach ([dirname($modkit_file_highlight_prefix), dirname($modkit_file_prefix)] as $target_path) 
                    {
                        if (!file_exists($target_path)) mkdir($target_path, 0777, true);
                    }
                    
                    foreach ([1, 2, "ungrouped"] as $hp)
                    {
                        $parser->copyFile("{$modkit_temp_folder}/{$cohort_ps_name}_{$hp}.bed", "{$modkit_file_prefix}_{$hp}.bed");
                        $parser->copyFile("{$modkit_temp_folder}/{$cohort_ps_name}_highlight_{$hp}.bed", "{$modkit_file_highlight_prefix}_{$hp}.bed");
                    }
                }

                //use temp files to generate plots and stats:
                $modkit_file_highlight_prefix = "{$modkit_temp_folder}/{$cohort_ps_name}_highlight";
                $modkit_file_prefix = "{$modkit_temp_folder}/{$cohort_ps_name}";
                
            }

            // for male samples: also generate an unphased methylation file for X/Y
            if ($only_one_haplotype && (!isset($cohort_folder) || !file_exists($modkit_file_highlight_prefix."_all.bed") || !file_exists($modkit_file_prefix."_all.bed")))
            {
                if ($tmp_bam == "") $tmp_bam = get_bam($cohort_bam_file, $cohort_ps_name, $chr, $gene_start, $gene_end, $ref_genome, $threads);
                $modkit_temp_folder = extract_methylation($tmp_bam, $chr, $gene_start, $gene_end, $start, $end, $cohort_ps_name, $ref_genome, $threads, true);
                if (isset($cohort_folder))
                {
                    $parser->copyFile("{$modkit_temp_folder}/{$cohort_ps_name}_all.bed", "{$modkit_file_prefix}_all.bed");
                    $parser->copyFile("{$modkit_temp_folder}/{$cohort_ps_name}_highlight_all.bed", "{$modkit_file_highlight_prefix}_all.bed");
                }
                
            }

            //get avg methylation
            $methylation = get_avg_methylation($modkit_file_highlight_prefix);
            if (!is_nan($methylation[1]["avg"])) $cohort_values[$identifier]["hp1"][] = $methylation[1]["avg"];
            if (!is_nan($methylation[2]["avg"])) $cohort_values[$identifier]["hp2"][] = $methylation[2]["avg"];
            if (!is_nan($methylation["all"]["avg"])) $cohort_values[$identifier]["all"][] = $methylation["all"]["avg"];

            //store raw data
            foreach ($methylation["all"]["fraction_modified"] as $pos => $all) 
            {
                if (isset($methylation[1]["fraction_modified"][$pos])) $raw_methylation_data[$identifier][$pos][$cohort_ps_name."_hp1"] = $methylation[1]["fraction_modified"][$pos];
                if (isset($methylation[2]["fraction_modified"][$pos])) $raw_methylation_data[$identifier][$pos][$cohort_ps_name."_hp2"] = $methylation[2]["fraction_modified"][$pos];
                $raw_methylation_data[$identifier][$pos][$cohort_ps_name."_all"] = $all;
            }

            // get cohort plot data
            if (!$skip_plot)
            {
                if ($only_one_haplotype)
                {
                    if (get_phasing_info($cohort_phasing_track_file, $chr, $gene_start, $gene_end, $gender) == "")
                    {
                        $cohort_plot_data[$identifier]["hp1_cohort_files"][] = $modkit_file_prefix."_all.bed";
                    }
                    else
                    {
                        $cohort_plot_data[$identifier]["hp1_cohort_files_fallback"][] = $modkit_file_prefix."_all.bed";
                    }
                }
                else
                {
                    if (get_phasing_info($cohort_phasing_track_file, $chr, $gene_start, $gene_end, $gender) == "")
                    {
                        if ($methylation["hp1_hp2_switched"])
                        {
                            $cohort_plot_data[$identifier]["hp1_cohort_files"][] = $modkit_file_prefix."_2.bed";
                            $cohort_plot_data[$identifier]["hp2_cohort_files"][] = $modkit_file_prefix."_1.bed";
                        }
                        else
                        {
                            $cohort_plot_data[$identifier]["hp1_cohort_files"][] = $modkit_file_prefix."_1.bed";
                            $cohort_plot_data[$identifier]["hp2_cohort_files"][] = $modkit_file_prefix."_2.bed";
                        }
                    }
                    else
                    {
                        if ($methylation["hp1_hp2_switched"])
                        {
                            $cohort_plot_data[$identifier]["hp1_cohort_files_fallback"][] = $modkit_file_prefix."_2.bed";
                            $cohort_plot_data[$identifier]["hp2_cohort_files_fallback"][] = $modkit_file_prefix."_1.bed";
                        }
                        else
                        {
                            $cohort_plot_data[$identifier]["hp1_cohort_files_fallback"][] = $modkit_file_prefix."_1.bed";
                            $cohort_plot_data[$identifier]["hp2_cohort_files_fallback"][] = $modkit_file_prefix."_2.bed";
                        }
                    }
                }
                
            }
            
        }

        //log excudeded samples:
        if (count($excluded_samples) > 0)
        {
            $buffer = array();
            foreach ($excluded_samples as $identifier => $sample_list) 
            {
                $buffer[] = "{$identifier}:\t".count($sample_list)."\t".implode(",", $sample_list);
            }
            $parser->log("The following samples were excluded for cohort computation due to phasing issues:", $buffer);
        }

    }

    //calcualte stats
    $cohort_mean_hp1 = array();
    $cohort_stdev_hp1 = array();
    $cohort_mean_hp2 = array();
    $cohort_stdev_hp2 = array();
    $cohort_mean_all = array();
    $cohort_stdev_all = array();
    $cohort_size = array();
    $zscore_hp1 = array();
    $zscore_hp2 = array();
    $zscore_all = array();
    foreach ($cohort_values as $identifier => $values) 
    {
        //TODO: filter too small cohort
        $cohort_mean_hp1[$identifier] = (count($values["hp1"]) > 0)?number_format(mean($values["hp1"]), 3):NAN;
        $cohort_stdev_hp1[$identifier] = (count($values["hp1"]) > 0)?number_format(stdev($values["hp1"]), 3):NAN;
        $cohort_mean_hp2[$identifier] = (count($values["hp2"]) > 0)?number_format(mean($values["hp2"]), 3):NAN;
        $cohort_stdev_hp2[$identifier] = (count($values["hp2"]) > 0)?number_format(stdev($values["hp2"]), 3):NAN;
        $cohort_mean_all[$identifier] = (count($values["all"]) > 0)?number_format(mean($values["all"]), 3):NAN;
        $cohort_stdev_all[$identifier] = (count($values["all"]) > 0)?number_format(stdev($values["all"]), 3):NAN;
        $cohort_size[$identifier] = count($values["hp1"]).", ".count($values["hp2"]).", ".count($values["all"]);

        //calculate zscore
        $zscore_hp1[$identifier] = NAN;
        $zscore_hp2[$identifier] = NAN;
        $zscore_all[$identifier] = NAN;
        if (floatval($cohort_stdev_hp1[$identifier]) != 0.0) $zscore_hp1[$identifier] = number_format((floatval($sample_methylation[$identifier]["hp1"]) - floatval($cohort_mean_hp1[$identifier])) / floatval($cohort_stdev_hp1[$identifier]), 3);
        if (floatval($cohort_stdev_hp2[$identifier]) != 0.0) $zscore_hp2[$identifier] = number_format((floatval($sample_methylation[$identifier]["hp2"]) - floatval($cohort_mean_hp2[$identifier])) / floatval($cohort_stdev_hp2[$identifier]), 3);
        if (floatval($cohort_stdev_all[$identifier]) != 0.0) $zscore_all[$identifier] = number_format((floatval($sample_methylation[$identifier]["all"]) - floatval($cohort_mean_all[$identifier])) / floatval($cohort_stdev_all[$identifier]), 3);
    }

    $meth_table->addComment("#Background cohort:\t".implode(",", $cohort_sample_list));
    $meth_table->addCol(array_values($cohort_mean_hp1), "cohort_mean_hp1");
    $meth_table->addCol(array_values($cohort_stdev_hp1), "cohort_stdev_hp1");
    $meth_table->addCol(array_values($zscore_hp1), "zscore_hp1");
    $meth_table->addCol(array_values($cohort_mean_hp2), "cohort_mean_hp2");
    $meth_table->addCol(array_values($cohort_stdev_hp2), "cohort_stdev_hp2");
    $meth_table->addCol(array_values($zscore_hp2), "zscore_hp2");
    $meth_table->addCol(array_values($cohort_mean_all), "cohort_mean_all");
    $meth_table->addCol(array_values($cohort_stdev_all), "cohort_stdev_all");
    $meth_table->addCol(array_values($zscore_all), "zscore_all");
    $meth_table->addCol(array_values($cohort_size), "cohort_size");

    $meth_table->toTSV($out);

    // create raw data table
    if ($export_methylation_data)
    {
        $header = array("#identifier", "region");
        $table = array();
        $sample_columns = array($name."_hp1", $name."_hp2", $name."_all");
        foreach ($cohort_sample_list as $ps_name) 
        {
            foreach (["_hp1", "_hp2", "_all"] as $hp) $ample_columns[] = $ps_name.$hp;
        }
        $header = array_merge($header, $ample_columns);

        $table[] = implode("\t", $header);

        foreach ($raw_methylation_data as $identifier => $methylation) 
        {
            $positions = array_keys($methylation);
            sort($positions);
            foreach ($positions as $pos) 
            {
                $line = array($identifier, $pos);
                foreach ($ample_columns as $col) 
                {
                    if (isset($methylation[$pos][$col])) $line[] = $methylation[$pos][$col];
                    else $line[] = "";
                }
                $table[] = implode("\t", $line);
            }

        }

        file_put_contents(dirname($out)."/".basename2($out)."_raw_data.tsv", implode("\n", $table));
    }

    //create cohort plotting commands
    if (!$skip_plot)
    {
        foreach ($cohort_plot_data as $identifier => $plot_data) 
        {

            $python_script = repository_basedir()."/src/Tools/generate_methylation_cohort_plot.py";

            $args = array();
            $args[] = "--hp1_file ".$plot_data["hp1_file"];
            if (count($plot_data["hp1_cohort_files"]) > 0) $args[] = "--hp1_cohort_files ".implode(" ", $plot_data["hp1_cohort_files"]);
            if (count($plot_data["hp1_cohort_files_fallback"]) > 0) $args[] = "--hp1_cohort_files_fallback ".implode(" ", $plot_data["hp1_cohort_files_fallback"]);
            
            if ($plot_data["unphased"])
            {
                $args[] = "--unphased";
            }
            else
            {
                $args[] = "--hp2_file ".$plot_data["hp2_file"];
                if (count($plot_data["hp2_cohort_files"]) > 0) $args[] = "--hp2_cohort_files ".implode(" ", $plot_data["hp2_cohort_files"]);
                if (count($plot_data["hp2_cohort_files_fallback"]) > 0) $args[] = "--hp2_cohort_files_fallback ".implode(" ", $plot_data["hp2_cohort_files_fallback"]);
            }
            $args[] = "--out ".$plot_data["out"];
            $args[] = "--site ".$plot_data["site"];
            $args[] = "--sample_name ".$plot_data["sample_name"];
            $args[] = "--highlight_start ".$plot_data["highlight_start"];
            $args[] = "--highlight_end ".$plot_data["highlight_end"];

            $in_files = array();
            $in_files[] = $python_script;
            $in_files[] = $plot_data["hp1_file"];
            if (isset($plot_data["hp2_file"])) $in_files[] = $plot_data["hp2_file"];
            $in_files = array_merge($in_files, $plot_data["hp1_cohort_files"]);
            $in_files = array_merge($in_files, $plot_data["hp2_cohort_files"]);
            $in_files = array_merge($in_files, $plot_data["hp1_cohort_files_fallback"]);
            $in_files = array_merge($in_files, $plot_data["hp2_cohort_files_fallback"]);
            $out_files = array("{$folder}/methylartist/");

            $jobs_plotting_cohort[] = array("Plotting_".$plot_data["site"]."_cohort", 
                $parser->execApptainer("python", "python3", repository_basedir()."/src/Tools/generate_methylation_cohort_plot.py ".implode(" ", $args), $in_files, $out_files, true));

        }
    }
}


# run plots in parallel (reverse array since plotting of REs takes longer)
if (count($jobs_plotting) > 0) 
{
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

if (count($jobs_plotting_cohort) > 0) $parser->execParallel(array_reverse($jobs_plotting_cohort), $threads, true, true, true, true); 

?>