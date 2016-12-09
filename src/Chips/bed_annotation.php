<?php	

/** 
	@page bed_annotation
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("bed_annotation", "\$Rev: 67 $", "Annotates a BED file with OMIM and UCSC data.");
$parser->addInfile("in", "BED input file.", false);
$parser->addOutfile("out",  "TSV output file.", false);
extract($parser->parse($argv));

//load required DBs
$genes = parse_ucsc_genes();
$omim = parse_omim_genes();

//load input
$regions = load_bed($in);

//annotated
$output = new Matrix();
$output->setHeaders(array("locus", "locus_size_mb", "gene_name", "gene_coordinates", "OMIM_title", "OMIM_disorders", "inheritance"));
foreach ($regions as $region)
{
	list($chr, $start, $end) = $region;
	$locus = $chr.":".$start."-".$end;
	$locus_size = number_format(($end-$start)/1000000, 3);
	$gene_list = array();
	$gene_positions = array();
	foreach($genes[$chr] as $gene)
	{
		list($gene_start, $gene_end, $gene_name) = $gene;
		
		if (range_overlap($start, $end, $gene_start, $gene_end))
		{
			$gene_list[] = $gene_name;
			$gene_positions[$gene_name] = $chr.":".$gene_start."-".$gene_end;
		}
	}
	
	$gene_list = array_unique($gene_list);
	
	foreach ($gene_list as $gene)
	{
		if (isset($omim[$chr][$gene]))
		{
			list ($title, $disorders) = $omim[$chr][$gene];
			
			if (count($disorders)==0)
			{
				$disorders[] = array("", "");
			}
			
			for ($i=0; $i<count($disorders); ++$i)
			{
				$output->addRow(array($locus, $locus_size, $gene, $gene_positions[$gene], $title, $disorders[$i][0], $disorders[$i][1]));
			}
		}
		else
		{
			$output->addRow(array($locus, $locus_size, $gene, $gene_positions[$gene], "", "", ""));
		}
	}
}

$output->sort(1, SORT_NUMERIC, SORT_DESC);

$output->toTSV($out);

?>