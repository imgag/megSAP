<?php
/** 
	@page create_qcml
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("create_qcml", "Collects QC terms and creates qcml files for somatic pipeline.");
$parser->addString("full_prefix", "Full filepath prefix for out files", false);
$parser->addString("t_basename", "Basename of tumor sample file including filepath", false);
$parser->addFlag("tumor_only", "Annotation for tumor-only analysis.");

extract($parser->parse($argv));

//Collect QC terms if necessary
$terms = array();
$sources = array();

//CNVs:
$som_clincnv = $full_prefix . "_clincnv.tsv"; //ClinCNV output file
if (file_exists($som_clincnv))
{
    $cnv_count_hq = 0;
    $cnv_count_hq_autosomes = 0;
    $cnv_count_loss = 0;
    $cnv_count_gain = 0;
    $h = fopen2($som_clincnv, 'r');
    while(!feof($h))
    {
        $line = trim(fgets($h));
        if ($line=="") continue;
        
        if ($line[0]!="#")
        {
            $parts = explode("\t", $line);
            $ll = $parts[11];
            if ($ll>=20) // TODO find somatic specific cut-off
            {
                ++$cnv_count_hq;
                
                $chr = $parts[0];
                if (is_numeric(strtr($chr, ["chr"=>""])))
                {
                    ++$cnv_count_hq_autosomes;
                    $cn = $parts[5];
                    if ($cn<2) ++$cnv_count_loss;
                    if ($cn>2) ++$cnv_count_gain;
                }
            }
        }
    }
    fclose($h);
    
    //counts (all, loss, gain)
    $terms[] = "QC:2000044\t{$cnv_count_hq}"; // somatic CNVs count
    if ($cnv_count_hq_autosomes>0)
    {
        $terms[] = "QC:2000118\t".number_format(100.0*$cnv_count_loss/$cnv_count_hq_autosomes, 2); // percentage losses
        $terms[] = "QC:2000119\t".number_format(100.0*$cnv_count_gain/$cnv_count_hq_autosomes, 2); // percentage gains
    }
    $sources[] = $som_clincnv;
}

// HRD score:
$hrd_file = $full_prefix."_HRDresults.txt";
if (file_exists($hrd_file))
{
    foreach (file($hrd_file) as $line)
    {
        if (starts_with($line , '""')) continue;
        list($sample, $loh, $tai, $lst, $hrd) = explode("\t", trim($line));
        $terms[] = "QC:2000062\t{$loh}";
        $terms[] = "QC:2000063\t{$tai}";
        $terms[] = "QC:2000064\t{$lst}";
        $terms[] = "QC:2000126\t{$hrd}";
    }
    $sources[] = $hrd_file;
}

//virus:
$viral = "{$t_basename}_viral.tsv"; // viral sequences results
if (file_exists($viral))
{
    $detected_viruses = [];
    foreach(file($viral) as $line)
    {
        if (starts_with($line, "#")) continue;
        
        list($chr, $start, $end, $v_name, $reads, $coverage, $coverage_rel, $mismatches, $ident) = explode("\t", $line);
        $v_base_name = explode("_", $v_name)[0];
        if (intval($reads) > 100 && ! in_array($v_base_name, $detected_viruses))
        {
            $detected_viruses[] = $v_base_name;
        }
    }
    
    $value = "None";
    if (count($detected_viruses) > 0)
    {
        $value = implode(", ", $detected_viruses);
    }
    
    $terms[] = "QC:2000130\t{$value}";
    $sources[] = $viral;
}

//MSI-status
$msi_o_file = $full_prefix . "_msi.tsv"; //MSI
if (!$tumor_only && file_exists($msi_o_file))
{
    foreach(file($msi_o_file) as $line)
    {
        if (trim($line) == "" || $line[0] == "#") continue;
        
        list($total, $somatic, $percent) = explode("\t", $line);
    }
    
    $terms[] = "QC:2000141\t{$percent}";
    $sources[] = $msi_o_file;
    
}

//HLA: TODO
// if (file_exists($hla_file_tumor))
// {
    
    // $sources[] = $hla_file_tumor;
// }

// if (!$single_sample && file_exists($hla_file_normal))
// {
    
    // $sources[] = $hla_file_normal;
// }

//TODO mutational signatures

//create qcML file
$qc_other = $full_prefix."_stats_other.qcML";
if (count($sources) > 0) 
{
    $tmp = $parser->tempFile("qc.tsv");
    file_put_contents($tmp, implode("\n", $terms));
    $in_files = $sources;
    $parser->execApptainer("ngs-bits", "TsvToQC", "-in $tmp -out $qc_other -sources ".implode(" ", $sources), $in_files, [dirname($qc_other)]);
}

?>