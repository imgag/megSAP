<?php

/** 
	@page create_circos_plot
 */
require_once(dirname($_SERVER['SCRIPT_FILENAME']) . "/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("create_circos_plot", "Creates a Circos plot with CN, CNVs, ROHs and BAFs of the current sample.");

$parser->addString("folder", "Analysis data folder.", false);
$parser->addString("name", "Base file name, typically the processed sample ID (e.g. 'GS120001_01').", false);

// optional
$parser->addString("build", "The genome build to use. The genome must be indexed for BWA!", true, "GRCh38");
$parser->addFloat("cn_offset", "Copy number between 2-offset and 2+offset will not be plotted to reduce the runtime of the script. Filter is deactivated if offset is 0.", true, 0.3);
$parser->addInt("cnv_min_ll", "Minimum loglikelyhood for CNVs to be shown in the plot.", true, 20);
$parser->addInt("cnv_min_nor", "Minimum number of regions for CNVs to be shown in the plot.", true, 3);
$parser->addFloat("cnv_max_af", "Maximum overlap with CNP regions (af_genomes_imgag column) for CNVs.", true, 0.5);
$parser->addFloat("cnv_min_length", "Minimimal CNV length (in kb) for CNVs to be shown in the plot.", true, 15.0);
$parser->addFloat("roh_min_length", "Minimimal ROH length (in kb) for ROHs to be shown in the plot.", true, 300.0);

extract($parser->parse($argv));

//init
$cn_offset = abs($cn_offset);

// get Cicos config and genome files 
$karyotype_file = repository_basedir() . "/data/misc/circos/karyotype.human.{$build}.txt"; 
if (!file_exists($karyotype_file)) 
{
    trigger_error("No karyotype file for build {$build} found!", E_USER_ERROR);
}
$circos_template_file = repository_basedir() . "/data/misc/circos/config_template.conf";
if (!file_exists($circos_template_file)) 
{
    trigger_error("No Circos template file found at \"$circos_template_file\"!", E_USER_ERROR);
}
$chr_file = repository_basedir() . "/data/misc/circos/chr_region.{$build}.txt";
if (!file_exists($chr_file)) 
{
    trigger_error("No chr region file found at \"$chr_file\"!", E_USER_ERROR);
}
$circos_housekeeping_file = repository_basedir() . "/data/misc/circos/housekeeping.conf";
if (!file_exists($circos_housekeeping_file)) 
{
    trigger_error("No Circos config file found at \"$circos_housekeeping_file\"!", E_USER_ERROR);
}

// get telomere file
$telomere_file = repository_basedir() . "/data/misc/centromer_telomer.bed";
if (!file_exists($telomere_file)) 
{
    trigger_error("No telomere file found at \"$telomere_file\"!", E_USER_ERROR);
}


// preprocess data
// parse CN file and generate temp file for circos plot

// create temporary working folder
$temp_folder = $parser->tempFolder("circos");

// clinCNV seg file:
$seg_file = "$folder/${name}_cnvs_clincnv.seg";
if (!file_exists($seg_file)) 
{
    // CnvHunter seg file:
    $seg_file = "$folder/${name}_cnvs.seg";
    if (!file_exists($seg_file)) 
    {
        trigger_error("No SEG file found!", E_USER_ERROR);
    }
}

$cn_temp_file = $temp_folder."/cn.seg";

$input_fh = fopen2($seg_file, "r");
$output_fh = fopen2($cn_temp_file, "w");
$n_copy_numbers = 0;

if ($input_fh) 
{
    while (($buffer = fgets($input_fh)) !== FALSE) 
    {
        // skip comments and header
        if (starts_with($buffer, "#")) continue;
        if (starts_with($buffer, "ID\tchr\tstart\tend")) continue;
        $row = explode("\t", $buffer);

        // skip undefined segments 
        if ((float) $row[5] < 0.0) continue;
        if (((float) $row[5] > (2 - $cn_offset)) && ((float) $row[5] < (2 + $cn_offset))) continue;
        if ($row[5] == "QC failed") continue;

        $cn = min(4.0, (float) $row[5]);

        //write modified seg file to temp
        fwrite($output_fh, $row[1] . "\t" . $row[2] . "\t" . $row[3] . "\t$cn\n");
        $n_copy_numbers++;
    }
    fclose($input_fh);
    fclose($output_fh);
}
else 
{
    if ($error)    trigger_error("Could not open file $seg_file.", E_USER_ERROR);
    return  "Could not open file $seg_file.";
}

// parse CNVs file and generate temp file for circos plot


// clinCNV CNV file:
$cnv_hunter = false;
$cnv_file = "$folder/${name}_cnvs_clincnv.tsv";
if (!file_exists($cnv_file)) 
{
    // CnvHunter CNV file:
    $cnv_file = "$folder/${name}_cnvs.tsv";
    $cnv_hunter = true;
    if (!file_exists($cnv_file)) 
    {
        trigger_error("No CNV file found!", E_USER_ERROR);
    }
}

$cnv_temp_file_del = $temp_folder."/cnv_del.tsv";
$cnv_temp_file_dup = $temp_folder."/cnv_dup.tsv";
$n_cnvs = 0;

// load CNV file as matrix
$cnv_matrix = Matrix::fromTSV($cnv_file);
// get column indices
$chr_idx = $cnv_matrix->getColumnIndex("chr");
$start_idx = $cnv_matrix->getColumnIndex("start");
$end_idx = $cnv_matrix->getColumnIndex("end");
if ($cnv_hunter)
{
    $cn_idx = $cnv_matrix->getColumnIndex("region_copy_numbers");
    $nor_idx = $cnv_matrix->getColumnIndex("region_count");
    $length_idx = $cnv_matrix->getColumnIndex("size");
}
else
{
    $cn_idx = $cnv_matrix->getColumnIndex("CN_change");
    $ll_idx = $cnv_matrix->getColumnIndex("loglikelihood");
    $nor_idx = $cnv_matrix->getColumnIndex("no_of_regions");
    $overlap_af_idx = $cnv_matrix->getColumnIndex("overlap af_genomes_imgag");
    $length_idx = $cnv_matrix->getColumnIndex("length_KB");
}

$output_fh_del = fopen2($cnv_temp_file_del, "w");
$output_fh_dup = fopen2($cnv_temp_file_dup, "w");
for ($row_idx=0; $row_idx < $cnv_matrix->rows(); $row_idx++) 
{ 
    // filter CNV
    // no_of_regions
    if ((int) $cnv_matrix->get($row_idx, $nor_idx) < $cnv_min_nor) continue;
    // CnvHunter:
    if ($cnv_hunter)
    {
        // CNV size
        if (((float) $cnv_matrix->get($row_idx, $length_idx) / 1000.0) < $cnv_min_length) continue;
    }
    else
    // ClinCNV
    {
        // loglikelihood
        if ((int) $cnv_matrix->get($row_idx, $ll_idx) < $cnv_min_ll) continue;
        // overlap with CNPs
        if ((float) $cnv_matrix->get($row_idx, $overlap_af_idx) > $cnv_max_af) continue;
        // CNV size
        if ((float) $cnv_matrix->get($row_idx, $length_idx) < $cnv_min_length) continue;
    }
    


    // write modified cnv files to temp
    if ((int) $cnv_matrix->get($row_idx, $cn_idx) < 2) 
    {
        fwrite($output_fh_del, $cnv_matrix->get($row_idx, $chr_idx)."\t".$cnv_matrix->get($row_idx, $start_idx)."\t".$cnv_matrix->get($row_idx, $end_idx)."\n");
    } 
    else 
    {
        fwrite($output_fh_dup, $cnv_matrix->get($row_idx, $chr_idx)."\t".$cnv_matrix->get($row_idx, $start_idx)."\t".$cnv_matrix->get($row_idx, $end_idx)."\n");
    }
    $n_cnvs++;
}

fclose($output_fh_del);
fclose($output_fh_dup);

// parse ROHs file and generate temp file for circos plot
$roh_file = "$folder/${name}_rohs.tsv";

$roh_temp_file = $temp_folder."/roh.tsv";
$n_rohs = 0;

if (!file_exists($roh_file)) 
{
    trigger_error("WARNING: No ROH file found!", E_USER_WARNING);
    
    // create empty temp file
    touch($roh_temp_file);
}
else
{
    $input_fh = fopen2($roh_file, "r");
    $output_fh = fopen2($roh_temp_file, "w");

    if ($input_fh) 
    {
        while (($buffer = fgets($input_fh)) !== FALSE) 
        {
            // skip comments and header
            if (starts_with($buffer, "#")) continue;
            $row = explode("\t", $buffer);

            // skip small ROHs 
            if ((float) $row[5] < $roh_min_length) continue;

            $cn = min(4.0, (float) $row[5]);

            // write modified roh file to temp
            fwrite($output_fh, $row[0] . "\t" . $row[1] . "\t" . $row[2] . "\n");
            $n_rohs++;
        }
        fclose($input_fh);
        fclose($output_fh);
    } 
    else 
    {
        if ($error) trigger_error("Could not open file $roh_file.", E_USER_ERROR);
        return  "Could not open file $roh_file.";
    }
}

// parse BAFs file and generate temp file for circos plot
$baf_file = "$folder/${name}_bafs.igv";

$baf_temp_file = $temp_folder."/bafs.tsv";
$n_bafs = 0;

if (!file_exists($baf_file)) 
{
    trigger_error("WARNING: No BAF file found!", E_USER_WARNING);

    // create empty temp file
    touch($baf_temp_file);
}
else
{
    $input_fh = fopen2($baf_file, "r");
    $output_fh = fopen2($baf_temp_file, "w");

    if ($input_fh) 
    {
        while (($buffer = fgets($input_fh)) !== FALSE) 
        {
            // skip comments and header
            if (starts_with($buffer, "#")) continue;
            if (starts_with($buffer, "Chromosome\tStart\tEnd")) continue;
            $row = explode("\t", $buffer);

            // write modified baf file to temp
            fwrite($output_fh, $row[0] . "\t" . $row[1] . "\t" . $row[2] . "\t" . $row[4] . "\n");
            $n_bafs++;
        }
        fclose($input_fh);
        fclose($output_fh);
    } 
    else 
    {
        if ($error) trigger_error("Could not open file $baf_file.", E_USER_ERROR);
        return  "Could not open file $baf_file.";
    }
}

// parse SV file and extract high-quality translocations
$sv_file = "$folder/${name}_manta_var_structural.bedpe";
$sv_filter = repository_basedir() . "/data/misc/circos/sv_filter.ini";

$sv_temp_file = $temp_folder."/sv.tsv";
$n_bnds = 0;

if (!file_exists($sv_file)) 
{
    trigger_error("WARNING: No SV file found!", E_USER_WARNING);

    // create empty temp file
    touch($sv_temp_file);
}
else
{
    // filter SV file
    list($stdout, $stderr, $return_code) = $parser->exec(get_path("ngs-bits")."SvFilterAnnotations", "-in $sv_file -out $sv_temp_file -filters $sv_filter", true, false, true);

    // abort if filter fails
    if($return_code != 0)
    {
        trigger_error("WARNING: SV file filtering failed! Skipping SV break points in circos plot.", E_USER_WARNING);
        // create empty temp file
        touch($sv_temp_file);
    }

    // cleanup filtered SV file
    $svs = file($sv_temp_file);
    $sv_positions = array();
    foreach ($svs as $sv_line) 
    {
        if (starts_with($sv_line, "#")) continue;
        $pos = array_slice(explode("\t", $sv_line), 0, 6);
        $sv_positions[] = implode("\t", $pos)."\n";
        $n_bnds++;
    }
    file_put_contents($sv_temp_file, $sv_positions);
}



// create dummy file to add sample name to plot
$sample_label = $temp_folder."/sample_label.txt";
$output_fh = fopen2($sample_label, "w");
fwrite($output_fh, "chr1\t1\t1000000\t$name\n");
fclose($output_fh);

// create modified Circos config file
$file_names = array();
$file_names["[OUTPUT_FOLDER]"] = $temp_folder;
$file_names["[PNG_OUTPUT]"] = "${name}_circos.png";
$file_names["[KARYOTYPE_FILE]"] = $karyotype_file;
$file_names["[CHR_FILE]"] = $chr_file;
$file_names["[TELOMERE_FILE]"] = $telomere_file;
$file_names["[BAF_FILE]"] = $baf_temp_file;
$file_names["[CN_FILE]"] = $cn_temp_file;
$file_names["[CNV_DUP_FILE]"] = $cnv_temp_file_dup;
$file_names["[CNV_DEL_FILE]"] = $cnv_temp_file_del;
$file_names["[ROH_FILE]"] = $roh_temp_file;
$file_names["[SV_BND_FILE]"] = $sv_temp_file;
$file_names["[SAMPLE_LABEL]"] = $sample_label;
$file_names["[LABEL_FOLDER]"] = repository_basedir()."/data/misc/circos";
$file_names["[HOUSEKEEPING_FILE]"] = $circos_housekeeping_file;

// parse circos template file and replace file names
$circos_config_file = $temp_folder."/circos_config.conf";

$input_fh = fopen2($circos_template_file, "r");
$output_fh = fopen2($circos_config_file, "w");

if (!$output_fh) 
{
    if ($error)    trigger_error("Could not open file $circos_config_file.", E_USER_ERROR);
    return  "Could not open file $circos_config_file.";
}
if ($input_fh) 
{
    while (($buffer = fgets($input_fh)) !== FALSE) 
    {
        // replace file names
        $modified_line = strtr($buffer, $file_names);

        // write modified baf file to temp
        fwrite($output_fh, $modified_line);
    }
    fclose($input_fh);
    fclose($output_fh);
} 
else 
{
    if ($error) trigger_error("Could not open file $circos_template_file.", E_USER_ERROR);
    return  "Could not open file $circos_template_file.";
}

// report preprocessing stats
trigger_error("Preprocessing finished.", E_USER_NOTICE);
trigger_error("Remaining CN entries: \t$n_copy_numbers", E_USER_NOTICE);
trigger_error("Remaining CNV entries: \t$n_cnvs", E_USER_NOTICE);
trigger_error("Remaining ROH entries: \t$n_rohs", E_USER_NOTICE);
trigger_error("Remaining BAF entries: \t$n_bafs", E_USER_NOTICE);
trigger_error("Remaining BND entries: \t$n_bnds", E_USER_NOTICE);


// create Circos plot
$perl_cpan = get_path("perl_cpan");
putenv("PERL5LIB=".$perl_cpan."/lib/perl5/:".getenv("PERL5LIB"));
$circos_bin = get_path("circos");
$parser->exec($circos_bin, "-nosvg -conf $circos_config_file", true);

// copy PNG to sample folder
$parser->moveFile("$temp_folder/${name}_circos.png", "$folder/${name}_circos.png");