<?php 
/** 
	@page an_spliceai
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

$parser = new ToolBase("an_spliceai", "Scores and annotates the SpliceAI prediction.");
$parser->addInfile("in",  "Input file in VCF format.", false);
$parser->addOutfile("out", "Output file in VCF format.", false);
//optional

$parser->addFloat("max_af", "Maximum gnomAD AF of variants that are scored.", true, 0.01);
$parser->addInt("max_vars", "The maximum number of variants to score.", true, 15000);
$parser->addString("build", "The genome build to use.", true, "GRCh38");
$parser->addInt("threads", "The maximum number of threads used.", true, 1);
extract($parser->parse($argv));

//returns the number of variants in a VCF file
function vcf_variant_count($vcf, $lines_containing=null)
{
	$count = 0;
	
	$h = fopen2($vcf, 'r');
	while(!feof($h))
	{
		$line = trim(fgets($h));
		if ($line=="" || $line[0]=='#') continue;
		
		if (is_null($lines_containing) || contains($line, $lines_containing)) ++$count;
	}
	fclose($h);
	
	return $count;
}

//filters VCF by AF and INFO key
function filter_by_annotation($vcf_in, $vcf_out, $max_af)
{
	//init
	$csq_gnomad_idx = -1;
	$c_written = 0;
	
	$h = fopen2($vcf_in, 'r');
	$h_o = fopen2($vcf_out, 'w');
	while(!feof($h))
	{
		$line = nl_trim(fgets($h));
		if ($line=="") continue;
		
		//header
		if(starts_with($line, "#"))
		{
			fwrite($h_o, $line."\n");
			
			//extract VEP entry indices
			if(starts_with($line, '##INFO=<ID=CSQ,'))
			{
                preg_match('/.*Description=\"(.*)\".*/', $line, $match);
				if(count($match) < 2) continue;
                $entries = explode("|", $match[1]);
				for($i=0; $i<count($entries); ++$i)
				{
					$entry = trim($entries[$i]);
					if($entry=='gnomAD_AF')
					{
						$csq_gnomad_idx = $i;
					}
				}
			}	
		}
		else //content
		{
			$parts = explode("\t", $line);
			if(count($parts) < 8) trigger_error("Invalid VCF content line in {$vcf_in}: Missing INFO column.", E_USER_ERROR);
			
			//skip already annotated variants
			$skip = false;
            $info = explode(';', $parts[7]);
			foreach($info as $info_entry)
			{
				$info_entry = trim($info_entry);
				if(starts_with($info_entry, 'SpliceAI='))
				{
					$skip = true;
				}
			}
			if ($skip) continue;
			
			//skip variant with too high AF
			$af = 0.0;
			foreach($info as $info_entry)
			{
				$info_entry = trim($info_entry);
				if(starts_with($info_entry, 'CSQ='))
				{
					$info_entry = explode('=', $info_entry, 2);
					$csq_annotations = explode('|', $info_entry[1]);

					//gnomAD AF (exome)
					$value = $csq_annotations[$csq_gnomad_idx];
					if (is_numeric($value))
					{
						$af = max($af, $value);
					}
				}
				//gnomAD AF (genome)
				else if(starts_with($info_entry, 'gnomADg_AF='))
				{
					$info_entry = explode('=', $info_entry, 2);
					$value = $info_entry[1];
					if (is_numeric($value))
					{
						$af = max($af, $value);
					}
				}
				//gnomAD AF (mito)
				else if(starts_with($info_entry, 'gnomADm_AF_hom='))
				{
					$info_entry = explode('=', $info_entry, 2);
					$value = $info_entry[1];
					if (is_numeric($value))
					{
						$af = max($af, $value);
					}
				}
			}
			if($af>$max_af) continue;
			
			//write
			fwrite($h_o, $line."\n");
			++$c_written;
		}
	}
	fclose($h);
	fclose($h_o);
	
	return $c_written;
}

//checks if contig lines are present. If not, adds contig lines. They are needed by SpliceAI.
function add_contigs_if_missing($vcf)
{
	global $parser;
	global $build;
	
	//check if contig lines are contained
	$contains_contig_lines = false;
	$line_below_reference_info = 1; //write contigs in second line if no reference genome line is given

	$h = fopen2($vcf, 'r');
	$line_nr = 0;
	while(!feof($h))
	{
		++$line_nr;
		$line = trim(fgets($h));
		if ($line=="") continue;
		if (!starts_with($line, "#")) break;

		if(starts_with($line, "##reference"))
		{
			$line_below_reference_info = $line_nr;
		}
		if (starts_with($line, "##contig"))
		{
			$contains_contig_lines = true;
		}
	}
	fclose($h);
	if($contains_contig_lines) return;
	
	$parser->log("Contig lines are missing. Adding them from genome FAI file!");
	
	//parse contig data from genome FAI file
	$new_contigs = array();
	$fai = file(genome_fasta($build).".fai");
	foreach($fai as $line)
	{
		$line = trim($line);
		if ($line=="") continue;
		
		list($chr, $len) = explode("\t", $line);
		if (contains($chr, "chrUn") || contains($chr, "_random") || contains($chr, "_decoy")) continue;
		$new_contigs[] = "##contig=<ID={$chr},length={$len}>\n";
	}
	
	//add contiig lines to VCF
	$lines = file($vcf);
	array_splice($lines, $line_below_reference_info, 0, $new_contigs);  
	file_put_contents($vcf, implode("", $lines));
}

function annotate_spliceai_scores($in, $vcf_filtered, $out)
{
	global $parser;
	global $build;
	global $threads;
	$ngsbits = get_path("ngs-bits");
	$splice_env = get_path("Splicing", true);
	
	//check if SpliceAI INFO header is already present in the input VCF
	$header_old = "";
	foreach(file($in) as $line)
	{
		$line = trim($line);
		if ($line=="") continue;
		if (!starts_with($line, "#")) break;
		
		if (starts_with($line, "##INFO=<ID=SpliceAI,"))
		{
			$header_old = $line;
		}
	}
	
	//call SpliceAI
	$tmp1 = $parser->tempFile("_spliceai_new_annotations.vcf");
	$args = array();
	$args[] = "-I {$vcf_filtered}";
	$args[] = "-O {$tmp1}";
	$args[] = "-R ".genome_fasta($build);
	$args[] = "-A ".strtolower($build);
	putenv("PYTHONPATH");
	$parser->exec("OMP_NUM_THREADS={$threads} {$splice_env}/splice_env/bin/python3 {$splice_env}/splice_env/lib/python3.6/site-packages/spliceai", implode(" ", $args), true);
	
	//no variants scored => copy input to output
	$var_count = vcf_variant_count($tmp1, "SpliceAI=");
	if($var_count==0)
	{
		$parser->copyFile($in, $out);
		return;
	}
	
	//sort, zip and index scored variants to make them usable with VcfAnnotateFromVcf
	$tmp2 = $parser->tempFile("_spliceai_new_annotations.vcf.gz");
	$parser->exec("{$ngsbits}/VcfSort", "-in {$tmp1} -out {$tmp1}");
	$parser->exec("bgzip", "-c $tmp1 > $tmp2");
	$parser->exec("tabix", "-f -p vcf $tmp2");
	
	//annotate variants in input with scored variants
	$tmp3 = $in;
	if ($header_old!="") //remove old header. Otherwise there will be two in the output
	{
		$tmp3 = $parser->tempFile("_input_header_fixed.vcf.gz");
		$parser->exec("grep", "-v '##INFO=<ID=SpliceAI,' {$in} > {$tmp3}");
	}
	$parser->exec("{$ngsbits}/VcfAnnotateFromVcf", "-in {$tmp3} -annotation_file {$tmp2} -info_ids SpliceAI -out {$out} -threads {$threads}");
	if ($header_old!="") //replace new by old header
	{		
		$header_old = str_replace("\"", "\\\"", $header_old); //SpliceAI header has doubles quotes which need to be escaped
		exec2("sed -i 's/##INFO=<ID=SpliceAI.*/{$header_old}/g' {$out}");
	}
}

//filter for private variants (low AF and not annotated in the step before)
$tmp1 = $parser->tempFile("_spliceai_filtered_annotations.vcf");
$var_count = filter_by_annotation($in, $tmp1, $max_af);
$parser->log("Variants after annotation filtering (AF too high or already annotated): $var_count");

//abort if no variants must be scored
if($var_count==0)
{
	$parser->copyFile($in, $out);
	return;
}

//filter for variants in SpliceAI transcript regions
$spliceai_regions = $parser->tempFile("spliceai_scoring_regions.bed");
$splice_env = get_path("Splicing", true);
exec2("cut -f 2,4,5 -d'\t' {$splice_env}/splice_env/lib/python3.6/site-packages/spliceai/annotations/".strtolower($build).".txt | sed 's/^/chr/' | sed '1d' > {$spliceai_regions}");
$tmp2 = $parser->tempFile("_spliceai_filtered_regions.vcf");
$parser->exec(get_path("ngs-bits")."/VcfFilter", "-reg {$spliceai_regions} -in {$tmp1} -out {$tmp2}", true);
$var_count = vcf_variant_count($tmp2);
$parser->log("Variants after SpliceAI transcript regions filter: {$var_count}");

//abort if no variants or too many variants must be scored
if($var_count==0 || $var_count > $max_vars)
{
	if ($var_count > $max_vars) trigger_error("SpliceAI annotation skipped: Number of private variants ({$var_count}) is above max_vars ({$max_vars})!", E_USER_NOTICE);
	$parser->copyFile($in, $out);
	return;
}

//perform splicing annotations
add_contigs_if_missing($tmp2);
annotate_spliceai_scores($in, $tmp2, $out);

?>