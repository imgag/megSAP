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
				//gnomAD AF (genome)
				if(starts_with($info_entry, 'gnomADg_AF='))
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
	
	//add contig lines to VCF
	$lines = file($vcf);
	array_splice($lines, $line_below_reference_info, 0, $new_contigs);  
	file_put_contents($vcf, implode("", $lines));
}

function annotate_spliceai_scores($in, $vcf_filtered, $out)
{
	global $parser;
	global $build;
	global $threads;
	global $transcript_annotations;
	
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
	$args[] = "-A ".$transcript_annotations;
	$args[] = "-M 1"; //enable masked scores
	$args[] = "-D 50";//TODO

	//set bind paths for container execution
	$in_files = [genome_fasta($build), $vcf_filtered, $transcript_annotations];

	//run spliceai container
	$spliceai_command = $parser->execApptainer("spliceai", "spliceai", implode(" ", $args), $in_files, [], true);
	$prefix = container_platform()=='apptainer' ? "APPTAINERENV" : "SINGULARITYENV";
	$parser->exec("{$prefix}_OMP_NUM_THREADS={$threads} {$spliceai_command}", "");

	//no variants scored => copy input to output
	$var_count = vcf_variant_count($tmp1, "SpliceAI=");
	if($var_count==0)
	{
		$parser->copyFile($in, $out);
		return;
	}
	
	//sort, zip and index scored variants to make them usable with VcfAnnotateFromVcf
	$tmp2 = $parser->tempFile("_spliceai_new_annotations.vcf.gz");
	$parser->execApptainer("ngs-bits", "VcfSort", "-in {$tmp1} -out {$tmp1}");
	$parser->execApptainer("htslib", "bgzip", "-c $tmp1 > $tmp2");
	$parser->execApptainer("htslib", "tabix", "-f -p vcf $tmp2");
	
	//annotate variants in input with scored variants
	$tmp3 = $in;
	if ($header_old!="") //remove old header. Otherwise there will be two in the output
	{
		$tmp3 = $parser->tempFile("_input_header_fixed.vcf.gz");
		$parser->exec("grep", "-v '##INFO=<ID=SpliceAI,' {$in} > {$tmp3}");
	}
	$parser->execApptainer("ngs-bits", "VcfAnnotateFromVcf", "-in {$tmp3} -source {$tmp2} -info_keys SpliceAI -out {$out} -threads {$threads}", [$tmp3], [dirname($out)]);
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

//create SpliceAI transcript regions BED file
$transcript_annotations = repository_basedir()."/data/misc/spliceai_gene_annotations/{$build}_original.txt";//TODO
$spliceai_regions = $parser->tempFile("spliceai_scoring_regions.bed");
$pipeline = [];
$pipeline[] = array("cut", "-f 2,4,5 -d'\t' {$transcript_annotations}");
$pipeline[] = array("sed", "'1d' > {$spliceai_regions}");
$parser->execPipeline($pipeline, "Dragen small variants post processing");

//filter for variants in SpliceAI transcript regions
$tmp2 = $parser->tempFile("_spliceai_filtered_regions.vcf");
$parser->execApptainer("ngs-bits", "VcfFilter", "-reg {$spliceai_regions} -in {$tmp1} -out {$tmp2} -ref ".genome_fasta($build), [genome_fasta($build)]);
$var_count = vcf_variant_count($tmp2);
$parser->log("Variants after SpliceAI transcript regions filter: {$var_count}");

//abort if no variants or too many variants must be scored
if($var_count==0 || $var_count > $max_vars)
{
	if ($var_count > $max_vars) trigger_error("SpliceAI annotation skipped: Number of private variants ({$var_count}) is zero or above max_vars ({$max_vars})!", E_USER_NOTICE);
	$parser->copyFile($in, $out);
	return;
}

//perform splicing annotations
add_contigs_if_missing($tmp2);
annotate_spliceai_scores($in, $tmp2, $out);

?>