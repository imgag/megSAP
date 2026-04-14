<?php
/** 
	@page create_qcml
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("create_qcml", "Collects QC terms and creates qcml files for somatic pipeline.");
$parser->addString("full_prefix", "Path prefix for out files (absolute path plus base file name).", false);
$parser->addString("viral_tsv", "Virus detection TSV file.", false);
$parser->addInfile("roi", "Processing system target region (needed for TMB calculation for tumor-only.", true);
$parser->addFlag("tumor_only", "Special handling for tumor-only: no MSI, etc.");
extract($parser->parse($argv));

function count_variants($gsvar_filtered, $roi)
{
	global $parser;
	$output = 0;
	
	list($stdout) = $parser->execApptainer("ngs-bits", "VariantFilterRegions", "-in {$gsvar_filtered} -reg {$roi}  -mode gsvar", [$roi], [], false, false);
	foreach($stdout as $line)
	{
		$line = nl_trim($line);
		if ($line=="" || $line[0]=='#') continue;
		
		$parts = explode("\t", $line);
		
		//filter
		if (contains($parts[7], "freq-tum")) continue;
		if (contains($parts[7], "GERMLINE")) continue;
		if (contains($parts[7], "RefCall")) continue;
		
		++$output;
	}
	
	return $output;
}

function calculate_tmb($gsvar_filtered, $roi)
{
	//no ROI > abort
	if ($roi=="") return "";
	
	//init
	global $parser;
	$exons = repository_basedir()."/data/gene_lists/gene_exons.bed";
	$blacklist = repository_basedir()."/data/gene_lists/somatic_tmb_blacklist.bed";
	$tmb_tsg = repository_basedir()."/data/gene_lists/somatic_tmb_tsg.bed";
	
	//determine usable target region
	$tmp = $parser->tempFile(".bed");
	$parser->execApptainer("ngs-bits", "BedIntersect", "-in {$roi} -in2 {$exons} -out {$tmp}", [$roi, $exons]);
	$roi_usable = $parser->tempFile(".bed");
	$parser->execApptainer("ngs-bits", "BedSubtract", "-in {$tmp} -in2 {$blacklist} -out {$roi_usable}", [$blacklist]);
	$roi_usable_size_mb = bed_size($roi_usable) / 1000000.0;
	if ($roi_usable_size_mb==0) return "";

	//determine usable TSG region
	$tmb_tsg_usable = $parser->tempFile(".bed");	
	$parser->execApptainer("ngs-bits", "BedIntersect", "-in {$tmb_tsg} -in2 {$roi_usable} -out {$tmb_tsg_usable}", [$tmb_tsg]);

	//determine exome sizes
	$exome_size_mb = bed_size($exons) / 1000000.0;
	
	//calculate mutation burden
	$somatic_var_count = count_variants($gsvar_filtered, $roi_usable);
	$somatic_var_count_tsg = count_variants($gsvar_filtered, $tmb_tsg_usable);
	$tmb =  ( ($somatic_var_count - $somatic_var_count_tsg) * $exome_size_mb / $roi_usable_size_mb + $somatic_var_count_tsg ) / $exome_size_mb;
	return number_format($tmb, 2);
}

//Collect QC terms if necessary
$terms = [];
$sources = [];

//CNVs
$som_clincnv = $full_prefix . "_clincnv.tsv"; //ClinCNV output file
if (file_exists($som_clincnv))
{
    $cnv_count_hq = 0;
    $cnv_count_hq_autosomes = 0;
    $cnv_count_loss = 0;
    $cnv_count_gain = 0;
    $h = fopen2($som_clincnv, 'r');
	$max_clonality = -1;
    while(!feof($h))
    {
        $line = trim(fgets($h));
        if ($line=="") continue;
        
        if ($line[0]!="#")
        {
            $parts = explode("\t", $line);
            $ll = $parts[11];
            if ($ll>=20) //TODO Alexander/Marc: find somatic specific cut-off for tumor-normal and tumor-only
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
			
			#tumor content estimate - max clonality
			$clonality = $parts[9];
			if (is_numeric($clonality) && floatval($clonality) > $max_clonality) $max_clonality = floatval($clonality);
        }
    }
    fclose($h);
    
	//tumor estimate:
	if ($max_clonality > 0) $terms[] = "QC:2000151\t{$max_clonality}";
	
    //counts (all, loss, gain)
    $terms[] = "QC:2000044\t{$cnv_count_hq}"; // somatic CNVs count
    if ($cnv_count_hq_autosomes>0)
    {
        $terms[] = "QC:2000118\t".number_format(100.0*$cnv_count_loss/$cnv_count_hq_autosomes, 2); // percentage losses
        $terms[] = "QC:2000119\t".number_format(100.0*$cnv_count_gain/$cnv_count_hq_autosomes, 2); // percentage gains
    }
    $sources[] = $som_clincnv;
}

//HRD score
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

//virus detection
if (file_exists($viral_tsv))
{
    $detected_viruses = [];
    foreach(file($viral_tsv) as $line)
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
    $sources[] = $viral_tsv;
}

//MSI status
$msi_o_file = $full_prefix . "_msi.tsv";
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

//TODO HLA + other mutational signatures

//tumor-only: small variants QC
if ($tumor_only)
{
	//filter GSvar for somatic variants
	$gsvar = "{$full_prefix}.GSvar";
	$filter = repository_basedir()."/data/misc/gsvar_filters/tumor_only.txt";
    $gsvar_filtered = $parser->tempFile(".GSvar");
	$parser->execApptainer("ngs-bits", "VariantFilterAnnotations", "-in {$gsvar} -filters {$filter} -out {$gsvar_filtered}", [$gsvar, $filter]);
	
	//parses filtered GSvar
	$c_gnomad = -1;
	$count_som_vars = 0;
	$somatic_indel_perc = 0;
	$known_som_vars_perc = 0;
	$ti_count = 0;
	$tv_count = 0;
	$h = fopen2($gsvar_filtered, 'r');
	while(!feof($h))
	{
		$line = nl_trim(fgets($h));
		if ($line=="") continue;

		$parts = explode("\t", $line);

		if ($line[0]=='#')
		{
			if (!starts_with($line, "##")) $c_gnomad = array_search("gnomAD", $parts);
			continue;
		}
		
		++$count_som_vars;
		
		//indel and TV ratio
		$ref = $parts[3];
		$alt = $parts[4];
		if (strlen($ref)>1 || strlen($alt)>1 || $ref=="-" || $alt=="-")
		{
			++$somatic_indel_perc;
		}
		else if (($alt=="A" && $ref=="G") || ($alt=="G" && $ref=="A") || ($alt=="T" && $ref=="C") || ($alt=="C" && $ref=="T"))
		{
			++$ti_count;
		}
		else
		{
			++$tv_count;
		}
		
		//known
		$af = $parts[$c_gnomad];
		if (is_numeric($af) && $af>0)
		{
			++$known_som_vars_perc;
		}
	}
	fclose($h);
	if ($c_gnomad==-1) trigger_error("Tumor-only GSvar file does not contain 'gnomAD' column!", E_USER_ERROR);
	
	$terms[] = "QC:2000041\t{$count_som_vars}";
	$known_som_vars_perc = number_format(100*$known_som_vars_perc/$count_som_vars, 2);
	$terms[] = "QC:2000045\t{$known_som_vars_perc}";
	$somatic_indel_perc = number_format(100*$somatic_indel_perc/$count_som_vars, 2);
	$terms[] = "QC:2000042\t{$somatic_indel_perc}";
	$terms[] = "QC:2000043\t".number_format($ti_count/$tv_count, 2);
	$tmb = calculate_tmb($gsvar_filtered, $roi);
	print "TMB: $tmb\n";
	if ($tmb!="") $terms[] = "QC:2000053\t{$tmb}";
}

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