<?php

/** 
	@page mip_generator
	@todo move calculation folder to temp and copy only files that are necessary
	@todo step 3: further filtering of MIPs, check if target region is covered sufficiently and update low-coverage statistics
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");
error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//
$parser = new ToolBase("mip_generator", "");
$parser->addInfile("target", "Target file (bed format).", false);
$parser->addString("project_name", "Prefix for files generated.", false);
$parser->addString("out_folder", "Out folder.", false);
//optional
$parser->addEnum("mode", "Select parameter set for MIP generation.", true,array("blood","ffpe","cfdna"),"blood");
$parser->addInt("min_capture_size", "Min. capture size.", true, 120);
$parser->addInt("max_capture_size", "Max. capture size.", true, 140);
$parser->addFloat("min_svr", "Minimum SVR score.", true, 1.4);
$parser->addFlag("no_svr", "Do not filter for svr score.");
extract($parser->parse($argv));

$temp_folder = $parser->tempFolder("MipGenerator_");

// 1. extract common SNPs
$output = array();
$tg_db = get_path('data_folder')."/dbs/1000G/1000g_v5b.vcf.gz";
$bed = file($target);
foreach($bed as $reg)
{
	list($c, $s, $e) = explode("\t", $reg);
	$c = substr($c, 3);
	$s -= 1000;
	$e += 1000;
	list($snps,) = $parser->exec("tabix", "$tg_db $c:$s-$e", false); //no output logging, because Toolbase::extractVersion() does not return
	foreach($snps as $snp)
	{
		$skip_variant = false;
		
		list(, , , , $obs, , , $info) = explode("\t", $snp);
		$info = explode(";", $info);
		unset($af);
		
		foreach($info as $entry)
		{
			if(starts_with($entry, "AF="))
			{
				$af = substr($entry, 3);
			}
		}
		
		if($af<0.01)
		{
			$skip_variant = true;
		}

		$obs = explode(",",$obs);
		foreach($obs as $o)
		{		
			if (starts_with($o,'<'))
			{
				$skip_variant = true;
			}
		}
		
		if(!$skip_variant)
		{
			$output[] = $snp;
		}
	}
}
//remove duplicates
$output = array_unique($output);
//store with headers
list($headers) = exec2("zcat $tg_db | head -1000 | egrep \"^#\"", false);	
$unsorted = $temp_folder."/unsorted.vcf";
file_put_contents($unsorted, implode("\n", $headers)."\n".implode("\n", $output));
//sort
$common_snps = $temp_folder."/common_snps.vcf";
$parser->exec(get_path('ngs-bits')."VcfStreamSort", "-in $unsorted -out $common_snps", true);
//index
$parser->exec("bgzip", "-f $common_snps", false); //no output logging, because Toolbase::extractVersion() does not return
$parser->exec("tabix", "$common_snps.gz", false); //no output logging, because Toolbase::extractVersion() does not return

// 2. generate MIPs
putenv("PATH=".dirname(get_path("samtools")).":".dirname(get_path("bwa")).":".getenv("PATH"));
$args = "-regions_to_scan $target ";
$args .= "-project_name $temp_folder/$project_name ";
$args .= "-bwa_genome_index ".get_path('data_folder')."/genomes/GRCh37.fa -bwa_threads 2 ";
$args .= "-snp_file $common_snps ";
$args .= "-arm_length_sums 40,41,42 -tag_sizes 8,0 -score_method mixed ";
$args .= "-min_capture_size $min_capture_size -max_capture_size $max_capture_size ";
if($mode=="ffpe")	$args .= "-half_seal_both_strands off -double_tile_strands_separately on ";
if($mode=="cfdna")	$args .= "-half_seal_both_strands off -double_tile_strands_separately on ";
$parser->exec(get_path('mipgen'), $args, true);
$p_mips = $out_folder."/".$project_name."_mips_picked.txt";
$parser->moveFile("$temp_folder/$project_name.picked_mips.txt", $p_mips);

// 3. filter MIPs and generate low-cov-statistics
$tmp_p_cov = $temp_folder."/tmp_mips_cov.bed";
$mips = Matrix::fromTSV($p_mips,"\t",">");
$i_svr = $mips->getColumnIndex("svr_score");
$i_chr = $mips->getColumnIndex("chr");
$i_start = $mips->getColumnIndex("mip_scan_start_position");
$i_end = $mips->getColumnIndex("mip_scan_stop_position");
$i_strand = $mips->getColumnIndex("probe_strand");
$i_name = $mips->getColumnIndex("mip_name");
$c_mips = 0;
$p_mips_new = $out_folder."/".$project_name."_mips_picked.tsv";
$mips_new = new Matrix();
$mips_new->setHeaders($mips->getHeaders());
//scores
$bed = new Matrix();
for($i=0;$i<$mips->rows();++$i)
{
	$row = $mips->getRow($i);
	if($row[$i_svr]<$min_svr && !$no_svr)	continue;
	$mips_new->addRow($row);
	$bed->addRow(array("chr".$row[$i_chr],$row[$i_start]-1,$row[$i_end],$row[$i_name],0,$row[$i_strand]));
	$c_mips += $row[$i_end]-($row[$i_start]-1);
}
$mips_new->toTSV($p_mips_new);
$p_cov = $out_folder."/".$project_name."_target_covered.bed";
$bed->toTSV($p_cov);
// statistics
$c_cov = 0;
$cov = Matrix::fromTSV($p_cov);
for($i=0;$i<$cov->rows();++$i)
{
	$c_cov += $cov->getRow($i)[2]-$cov->getRow($i)[1];
}

// 4. missing regions
$p_miss = $out_folder."/".$project_name."_target_missing.bed";
$c_miss = 0;
$parser->exec(get_path("ngs-bits")."/BedSubtract","-in $target -in2 $p_cov -out $p_miss", true);
$miss = Matrix::fromTSV($p_miss);
for($i=0;$i<$miss->rows();++$i)
{
	$c_miss += $miss->getRow($i)[2]-$miss->getRow($i)[1];
}
$parser->exec(get_path("ngs-bits")."/BedAnnotateGenes","-in $p_miss -out $p_miss", true);

$parser->log("Pre-Filter used for MIPs: ".($no_svr?"none ":"Minimum 'svr=$min_svr', ").". Used mode '$mode'");
$parser->log("Total number of MIPs:\t".($mips_new->rows()).".");
$parser->log("Total length of given target region:\t".($c_cov+$c_miss)." bp.");
$parser->log("Total length of (filtered) MIPs:\t$c_mips bp.");
$parser->log("Target region missed by (filtered) MIPs:\t$c_miss bp (".number_format($c_miss/($c_cov+$c_miss)*100,2)." %).");
$parser->log("Target region covered by (filtered) MIPs:\t$c_cov bp (".number_format($c_cov/($c_cov+$c_miss)*100,2)." %).");
?>