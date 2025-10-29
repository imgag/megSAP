<?php
/** 
	@page create_baf_file
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("create_baf_file", "Creates BAF file from VCF file.");
$parser->addString("out_folder", "Target folder to store BAF file.", false);
//optional
$parser->addInfile("vcf", "Input VCF file.", true);
$parser->addInfile("sample_list", "TXT File with multiple sample names to create baf file for", true);
$parser->addInfile("t_bam", "Input tumor bam/cram.", true);
$parser->addInfile("n_gsvar", "Input normal GSvar file for tumor sample baf creation.", true);
$parser->addString("sample_name", "Processed sample name. Inferred from vcf file if none is given", true);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
$parser->addString("old_baf_folder", "Original baf folder", true);
$parser->addInt("mq_thresh", "Sets a mapping quality threshold for variants to be flagged as n/a in BAF", true, 1);
$parser->addFlag("tumor", "flag for baf creation of tumor samples. Needs normal GSvar file and tumor bam!");
$parser->addFlag("test", "test flag to skip out_file existence check");
extract($parser->parse($argv));

function get_mq($cols, $vcf) 
{
    $info = explode(";", $cols[7]);

    $mqm = null;
    foreach ($info as $field) 
    {
        if (strpos($field, "MQM=") === 0) 
        {
            $mqm = floatval(substr($field, 4));
            break;
        }
    }
    if ($mqm === null) 
    {
        foreach ($info as $field) {
            if (strpos($field, "MQ=") === 0) 
            {
                $mqm = floatval(substr($field, 3));
                break;
            }
        }
    }

	return $mqm;
}

function create_baf($vcf, $out_file, $s_col)
{
	global $mq_thresh;

	$in_handle  = gzopen2($vcf,"r");
	$out_handle = fopen2($out_file,"w");
	//process VCF header (extract variant caller and sample index)
	$var_caller = false;
	$source_lines = 0;
	$sample_idx = false;
	while(!feof($in_handle))
	{
		$line = nl_trim(fgets($in_handle));
		if($line=="") continue;
		
		//variant caller
		if(starts_with($line, "##source="))	
		{
			++$source_lines;
			$line = strtolower($line);
			if(contains($line, "freebayes")) $var_caller = "freebayes";
			if(contains($line, "strelka")) $var_caller = "strelka";
			if(contains($line, "umivar2")) $var_caller = "umivar2";
			if(contains($line, "dragen"))
			{
				if ($var_caller != "dragen") $var_caller = "dragen";
				else --$source_lines;
			}
			if(contains($line, "deepvariant")) $var_caller = "deepvariant";
			if(contains($line, "varscan2")) $var_caller = "varscan2";

		}

		//workaround for dragen vcfs with no source line
		if(!$var_caller)
		{
			if(contains($line, "DRAGENCommandLine")) 
			{
				$var_caller = "dragen";
				++$source_lines;
			}
		}
		
		//header line
		if(starts_with($line, "#CHROM")) 
		{
			$parts = explode("\t", $line);
			$sample_idx = array_search($s_col, $parts);
			break;
		}
	}
	if($source_lines!=1) trigger_error("Found ".($source_lines==0 ? "no" : "several")." 'source' entries in VCF header of vcf '{$vcf}' (needed to identify the variant caller).", E_USER_ERROR);
	if($var_caller===false) trigger_error("Unknown variant caller from 'source' entry in VCF header.", E_USER_ERROR);
	if($sample_idx===false) trigger_error("Could not identify sample column '$s_col' in VCF header.", E_USER_ERROR);
	
	//process variants
	while(!feof($in_handle))
	{
		$line = trim(fgets($in_handle));
		if(starts_with($line,"#")) continue;
		if(empty($line)) continue;
		
		$cols = explode("\t", $line);
		if (count($cols)<9) trigger_error("VCF file line contains less than 10 columns: $line\n", E_USER_ERROR);
		list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format) = $cols;
		
		if(strlen($ref) != 1 || strlen($alt) != 1 ) continue; //Skip INDELs
		$mq = get_mq($cols, $vcf);

		//calculate DP/AF
		if ($mq !== null && $mq < $mq_thresh)
		{
			$dp = 0;
			$af = "n/a";
		}
		else if ($chr=="chrMT") //special handling of mito
		{
			$formats = explode(":",$format);
			if (in_array("AO", $formats)) list($dp, $af) = vcf_freebayes($format, $cols[$sample_idx], $pos);
			else if (in_array("AF", $formats)) list($dp, $af) = vcf_dragen_var($format, $cols[$sample_idx], $pos);
		}
		else if (str_contains($filter, "low_mappability")) //special handling of low_mappability variants (might be called with freebayes or deepvariant)
		{
			$formats = explode(":",$format);
			if (in_array("AO", $formats)) list($dp, $af) = vcf_freebayes($format, $cols[$sample_idx], $pos);
			else if (in_array("AF", $formats) || in_array("VAF", $formats)) list($dp, $af) = vcf_deepvariant($format, $cols[$sample_idx]);
		}
		else if($var_caller=="freebayes" || str_contains($filter, "mosaic")) //mosaic is always called with freebayes
		{
			list($dp, $af) = vcf_freebayes($format, $cols[$sample_idx], $pos);
		}
		else if($var_caller=="strelka")
		{
			list($dp, $af) = vcf_strelka_snv($format, $cols[$sample_idx], $alt);
		}
		else if ($var_caller == "dragen")
		{
			list($dp, $af) = vcf_dragen_var($format, $cols[$sample_idx], $pos);
		}
		else if($var_caller == "umivar2")
		{
			list($dp, $af, $p_value, $alt_count, $strand, $seq_context, $homopolymer, $m_af, $m_ref, $m_alt) = vcf_umivar2($filter, $info, $format, $cols[$sample_idx]);
		}
		else if ($var_caller =="deepvariant")
		{
			list($dp, $af) = vcf_deepvariant($format, $cols[$sample_idx]);
		}
		else if ($var_caller == "varscan2")
		{
			list($dp, $af) = vcf_varscan2($format, $cols[$sample_idx]);
		}

		fputs($out_handle,"{$chr}\t{$pos}\t{$pos}\t{$chr}_{$pos}\t{$af}\t{$dp}\n");
	}
	gzclose($in_handle);
	fclose($out_handle);
}

//get genome
$genome = genome_fasta($build);

if (!$tumor)
{
	if (isset($vcf))
	{
		//check that vcf OR sample_list is given as input
		if (isset($sample_list)) trigger_error("You can either provide a single 'vcf' or a 'sample_list' as input, not both", E_USER_ERROR);

		//check that vcf file exists
		if (!file_exists($vcf))
		{
			trigger_error("VCF file '{$vcf}' not found. Cannot generate BAF file for this sample!", E_USER_WARNING);
			return;
		} 

		//get sample name from input vcf
		is_null($sample_name) ? $s_col = basename($vcf, "_var.vcf.gz") : $s_col = $sample_name;
		$out_file = $out_folder."/{$s_col}.tsv";	

		//Abort if out_file exists to prevent interference with other jobs
		if(file_exists($out_file) && !$test) return;

		//generate baf file
		create_baf($vcf, $out_file, $s_col);
	}
	elseif (isset($sample_list))
	{
		if (!file_exists($sample_list)) trigger_error("Given list of samples is not readable: $sample_list", E_USER_ERROR);

		$in_list_handle = fopen($sample_list, "r");
		while(!feof($in_list_handle))
		{
			$line = trim(fgets($in_list_handle));
			$out_file = $out_folder."/{$line}.tsv";

			//Abort if out_file exists to prevent interference with other jobs
			if(file_exists($out_file) && !$test) continue;

			if(empty($line)) continue;
			if(file_exists($out_file) && !$test) continue;

			list($stdout, $stderr) = execApptainer("ngs-bits", "SamplePath", "-ps $line");
			$ps_folder = trim(implode("", $stdout));
			$vcf = $ps_folder."/{$line}_var.vcf.gz";

			if (!file_exists($vcf)) 
			{
				if (isset($old_baf_folder))	$parser->exec("cp","{$old_baf_folder}/{$line}.tsv {$out_file}");
				else trigger_error("Vcf '{$vcf}' not found for sample '{$line}' in sample folder '{$ps_folder}', skipping sample", E_USER_WARNING);
			}
			else create_baf($vcf, $out_file, $line);
		}
	}
	else
	{
		trigger_error("Either 'vcf' or 'sample_list' must be given as parameter", E_USER_ERROR);
	}
}
else
{
	$tmp_out = $parser->tempFile(".tsv");
	$parser->execApptainer("ngs-bits", "VariantAnnotateFrequency", "-in {$n_gsvar} -bam {$t_bam} -depth -out {$tmp_out} -ref {$genome}", [$n_gsvar, $t_bam, $genome]);

	$in_handle  = fopen2($tmp_out,"r");
	$out_handle = fopen2($out_file,"w");

	while(!feof($in_handle))
	{
		$line = trim(fgets($in_handle));
		if(starts_with($line,"#")) continue;
		if(empty($line)) continue;
		
		$parts = explode("\t",$line);
		list($chr,$start,$end,$ref,$alt) = $parts;
		if(strlen($ref) != 1 || $ref == "-" || $alt == "-" || strlen($alt) != 1 ) continue; //Skip INDELs

		$freq = $parts[count($parts)-2]; //freq annotation is 2nd last col
		$depth= $parts[count($parts)-1]; //depth annotation is last col
		fputs($out_handle,"{$chr}\t{$start}\t{$end}\t{$chr}_{$start}\t{$freq}\t{$depth}\n");
	}
	fclose($in_handle);
	fclose($out_handle);
}
?>