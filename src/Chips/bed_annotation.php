<?php	

/** 
	@page bed_annotation
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("bed_annotation", "Annotates a BED file with gene and OMIM information.");
$parser->addInfile("in", "BED input file.", false);
$parser->addOutfile("out",  "TSV output file.", false);
extract($parser->parse($argv));

//annotate gene info
$temp = temp_file("_anno.bed");
$parser->exec(get_path("ngs-bits")."BedAnnotateGenes","-in $in -out $temp", true);

//load OMIM file
$omim = array();
$file = file(get_path("data_folder")."/dbs/OMIM/omim.bed");
foreach ($file as $line)
{
	$line = nl_trim($line);
	if ($line=="" || $line[0]=="#") continue;
	
	list($chr, $start, $end, $anno) = explode("\t", $line);
	$anno = strtr($anno, array('_'=>' ', '['=>'', ']'=>''));
	list($mim, $gene, $details) = explode(' ', $anno, 3);
	$omim[$chr][] = array($start, $end, strtoupper($gene), $mim, $details);
}

//annotate OMIM info and write as TSV
$output = new Matrix();
$output->setHeaders(array("locus", "locus_size_mb", "gene", "OMIM_number", "OMIM_details"));
$file = file($temp);
foreach ($file as $line)
{
	$line = nl_trim($line);
	if ($line=="" || $line[0]=='#' || starts_with($line, "track") || starts_with($line, "browser")) continue;
	
	list($chr, $start, $end, $genes) = explode("\t", $line);
	
	//extract omim info per gene
	$omim_info = array();
	foreach($omim[$chr] as list($o_start, $o_end, $gene, $mim, $details))
	{
		if (range_overlap($start, $end, $o_start, $o_end))
		{
			$omim_info[$gene] = array($mim, $details);
		}
	}
	
	//merge genes with OMIM genes
	$genes = explode(", ", strtoupper($genes));
	$genes = array_merge($genes, array_keys($omim_info));
	$genes = array_unique($genes);
	
	//write one line per gene
	$locus = $chr.":".$start."-".$end;
	$locus_size = number_format(($end-$start)/1000000, 3);
	foreach($genes as $gene)
	{
		@list($mim, $details) = $omim_info[$gene];
		$output->addRow(array($locus, $locus_size, $gene, $mim, $details));
	}
}
$output->sort(1, SORT_NUMERIC, SORT_DESC);
$output->toTSV($out);


?>