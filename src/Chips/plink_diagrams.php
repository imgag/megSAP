<?php

/** 
	@page plink_diagrams
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

/**
	@brief Shows an interactive 2D plot using gnuplot
	
	@param labels An associative array of the labels used for the plot. The following keys are supported:
				- title
				- xlabel
				- ylabel
	@param data_labels An array of lables for the plot data (number must match the number of labels).
	@param data An array of one-dimensional data arrays to plot. 
	@param options An associative array of options The following keys are supported:
				- file If set a png file with this name is created instead of a interactive plot.
				- xrange X-axis range in format '0:100'.
				- yrange Y-axis range in format '0:100'.

	@note Several datasets can be plotted by handing an array of data arrays as the $data argument
	
	@ingroup plotting
*/
function plot($labels, $data_labels, $data, $options = null)
{
	//init labels if not set (to avoid warnings)
	if (!isset($labels["title"])) $labels["title"]="";
	if (!isset($labels["xlabel"])) $labels["xlabel"]="";
	if (!isset($labels["ylabel"])) $labels["ylabel"]="";
	
	//check input
	if (!is_array($data[0]))
	{
		$data = array($data);
	}
	
	if (!is_array($labels))
	{
		warning("Invalid first argument. Must be an array!", __FILE__, __LINE__, __FUNCTION__);
	}
	if (!is_array($data_labels))
	{
		warning("Invalid second argument. Must be an array!", __FILE__, __LINE__, __FUNCTION__);
	}
	if (count($data_labels)!=count($data))
	{
		warning("Wrong label count", __FILE__, __LINE__, __FUNCTION__);
	}
	
	//store datasets
	$tmp_file = temp_file();
	
	$i=0;
	foreach ($data as $dataset)
	{
		file_put_contents($tmp_file."_".$i.".txt", implode("\n", $dataset));
		++$i;
	}
	
	//handle options
	if (isset($options["file"]))
	{
		$file = $options["file"];
		if (!ends_with($file, ".png"))
		{
			$file .= ".png";
		}
	}
	
	//create gnuplot script
	$script = array();
	if (isset($file))
	{
		$script[] = "set terminal png size 1024,768";
		$script[] = "set output '$file'";
	}
	$script[] = "set title '".$labels["title"]."'";
	$script[] = "set xlabel '".$labels["xlabel"]."'";
	$script[] = "set ylabel '".$labels["ylabel"]."'";
	if (isset($options["xrange"]))
	{
		$script[] = "set xrange [".$options["xrange"]."]";
	}
	if (isset($options["yrange"]))
	{
		$script[] = "set yrange [".$options["yrange"]."]";
	}

	$line = "plot ";
	$i=0;
	foreach ($data as $dataset)
	{
		if ($i!=0) $line .= ", ";
		$line .= "'".$tmp_file."_".$i.".txt' title '";
		if (isset($data_labels[$i]))
		{
		  $line .= $data_labels[$i];
		}
		$line .= "' with lines";
		++$i;
	}
	$script[] = $line."";
	
	//store script
	file_put_contents($tmp_file.".gnuplot", implode("\n",$script));
	
	//execute script
	exec("gnuplot ".$tmp_file.".gnuplot");
}

// parse command line arguments
$parser = new ToolBase("plink_diagrams", "Creates diagrams for all chromosomes showing homocygosity regions of affected and non-affected indivisuals.");
$parser->addInfile("in",  "PLINK homocygosity mapping output file.", false);
$parser->addInt("bin_size",  "Bin size in bases (one bin is displayed as one pixel).", true, 100000);
$parser->addOutfile("out",  "Diagram output file base name. Default is family ID if unset.", true);
extract($parser->parse($argv));

$file = file($in);

// determine family ID
$parts = explode("\t", $file[2]);
$fam = $parts[0];

// determine output file base name
if (!isset($out))
{
	$out = $fam;
}

//count number of affected and not affected
$a = array();
$na = array();
$seen = array();
foreach($file as $line)
{
	$parts = explode("\t", $line);
	$iid = $parts[1];
	if (!in_array($iid, $seen))
	{
		if ($parts[2]==2)
		{
			$a[] = $iid;
		}
		else if ($parts[2]==1)
		{
			$na[] = $iid;
		}
		$seen[] = $iid;
	}
}
$all = array_merge($a, $na);
$count_a = count($a);
$count_na = count($na);

//create diagrams
$chrs = chr_list();
foreach($chrs as $chr)
{
	$bps = chr_info($chr);
	
	$bins = array();
	foreach($all as $iid)
	{
		$bins[$iid] = array_fill(0, $bps/$bin_size+1, 0);
	}
	
	foreach($file as $line)
	{
		$parts = explode("\t", $line);
		if ($parts[3]==$chr)
		{
			$iid = $parts[1];
			$start = $parts[6];
			$end = $parts[7];
			
			for($i=floor($start/$bin_size); $i<=floor($end/$bin_size); ++$i)
			{
				$bins[$iid][$i] = 1;
			}			
		}
	}
	
	// sum up affected and not affected (normalize with a/na counts)
	$bins_a = array_fill(0, $bps/$bin_size+1, 0);
	$bins_na = array_fill(0, $bps/$bin_size+1, 0);
	for($i=0; $i<count($bins_a); ++$i)
	{
		$a_sum = 0;
		$na_sum = 0;
		
		foreach($all as $iid)
		{
			if (in_array($iid, $a))
			{
				$a_sum += $bins[$iid][$i];
			}
			else if (in_array($iid, $na))
			{
				$na_sum += $bins[$iid][$i];
			}
		}
		
		$bins_a[$i]  = 0.5 + 99 * $a_sum  / $count_a;
		$bins_na[$i] = 0;
		if ($count_na!=0)
		{
			$bins_na[$i] = 0.5 + 99 * $na_sum / $count_na;
		}
	}
	
	// plot
	$options = array(
		"file" => $out."_chr".$chr,
		"yrange" => "0:100"
		);
	
	$labels = array(
		"title"=>"$fam chr$chr",
		"xlabel"=>"coordinate [x$bin_size bases]",
		"ylabel"=>"%"
		);
		
	plot($labels, array("affected (".count($a).")", "not affected (".count($na).")"), array($bins_a, $bins_na), $options);
}

?>