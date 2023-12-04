<?php

/** 
	@page mip_generator
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
$parser->addInt("overlap", "MIP overlap size (start value).", true, 0);
$parser->addString("snp_db", "Common variants database, which contains 'AF' info annoation.", true, get_path("data_folder")."/dbs/gnomAD/gnomAD_genome_v3.1.2_GRCh38.vcf.gz"); 
extract($parser->parse($argv));

$temp_folder = $parser->tempFolder("MipGenerator_");

//check db file
if (!file_exists($snp_db)) trigger_error("VCF given for parameter '-snp_db' ('${snp_db}') not found!", E_USER_ERROR);
if (!file_exists($snp_db.".tbi") && !file_exists($snp_db.".csi")) trigger_error("VCF index file for parameter '-snp_db' is missing!", E_USER_ERROR);

// 1. extract common SNPs
$output = array();
$bed = file($target);
foreach($bed as $reg)
{
	list($c, $s, $e) = explode("\t", $reg);
	$e = trim($e);
	//GRCh38 requires chr prefix: do not cut
	//$c = substr($c, 3);
	$s -= 1000;
	$e += 1000;

	list($snps,) = $parser->exec("tabix", "$snp_db $c:$s-$e", false); //no output logging, because Toolbase::extractVersion() does not return

	foreach($snps as $snp)
	{
		list(, , , , $obs, , , $info) = explode("\t", $snp);
		
		//skip low AF variants
		$af = 0.0;
		$info = explode(";", $info);
		foreach($info as $entry)
		{
			if(starts_with($entry, "AF="))
			{
				$af = substr($entry, 3);
			}
		}
		if($af<0.01)
		{
			continue;
		}
		
		//skip CNVs etc.
		$obs = explode(",", $obs);
		foreach($obs as $o)
		{		
			if (starts_with($o,'<'))
			{
				continue;
			}
		}
		
		$output[] = $snp;
	}
}
//remove duplicates
$output = array_unique($output);
//store with headers
list($headers) = exec2("zcat $snp_db | head -1000 | egrep \"^#\"", false);	
$unsorted = $temp_folder."/unsorted.vcf";
file_put_contents($unsorted, implode("\n", $headers)."\n".implode("\n", $output));
//sort
$common_snps = $temp_folder."/common_snps.vcf";
$parser->exec(get_path('ngs-bits')."VcfStreamSort", "-in $unsorted -out $common_snps", true);
//index
$parser->exec("bgzip", "-f $common_snps", false); //no output logging, because Toolbase::extractVersion() does not return
$parser->exec("tabix", "$common_snps.gz", false); //no output logging, because Toolbase::extractVersion() does not return
$common_snps = "${common_snps}.gz";

// 2. generate MIPs
putenv("PATH=".dirname(get_path("samtools")).":".dirname(get_path("bwa")).":".getenv("PATH"));
$args = "-regions_to_scan $target ";
$args .= "-project_name $temp_folder/$project_name ";
$args .= "-bwa_genome_index ".genome_fasta("GRCh38")." -bwa_threads 2 ";
$args .= "-snp_file $common_snps ";
$args .= "-starting_mip_overlap $overlap ";
$args .= "-arm_length_sums 40,41,42 -tag_sizes 8,0 -score_method mixed ";
$args .= "-min_capture_size $min_capture_size -max_capture_size $max_capture_size ";
$args .= "-capture_increment 1 " ;
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
$i_e_copies = $mips->getColumnIndex("ext_probe_copy");
$i_l_copies = $mips->getColumnIndex("lig_probe_copy");
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
	$info = sprintf("Name=%s;SVR=%.2f;ext_probe_copy=%d;lig_probe_copy=%d", $row[$i_name], $row[$i_svr], $row[$i_e_copies], $row[$i_l_copies]);
	$color = "22,111,255";
	if($row[$i_svr]<$min_svr) $color = "255,53,22";
	$bed->addRow(array("chr".$row[$i_chr],$row[$i_start]-1,$row[$i_end],$info,$row[$i_svr],$row[$i_strand],$row[$i_start]-1,$row[$i_start]-1,$color));
	$c_mips += $row[$i_end]-($row[$i_start]-1);
}
$mips_new->toTSV($p_mips_new);
$p_cov = $out_folder."/".$project_name."_target_covered.bed";
$bed->addComment("gffTags");
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