<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("converter_bgi2tsv", "Compares two tsv files.");
$parser->addInfile("gff",  "Input file / folder containing variant files (if folder: files are supposed to end with '*.gff.gz').", false);
$parser->addOutfile("out",  "Output file in tsv-format.", false);
extract($parser->parse($argv));

//get gzipped files in folder
$files = array($gff);
if(is_dir($gff))
{
	exec("find $gff -name '*.gff.gz'", $files);
}

//convert bgi to tsv (multiple gzipped files)
$headers = array('chr' => 'chr', 'start' => 'start', 'end' => 'end', 
				 'ref' => 'ref', 'alleles' => 'obs', 'mutType' => 'genotype', 
				 'score' => 'snp_q', 'type' => 'type', 'ID' => 'ID', 'state' => 'state', 
				 'sample' => 'sample', 'function' => 'function',
				 'gene' => 'gene', 'rest' => 'rest');

$handle1 = fopen2($out, "w");
fwrite($handle1, "#".implode("\t", $headers)."\n");
foreach($files as $file)
{
	
	//read file
	$handle2 = gzopen2($file , "r");
	while (!feof($handle2))
	{
		$line = nl_trim(fgets($handle2));
		if ($line=="") continue;
		$cols = explode("\t", $line);
		
		$row = array();
		
		for($i = 0; $i <= 7; ++$i)
		{
			$row[] = array_shift($cols);
		}
		$row[] = implode("", $cols);
		
		$fields = array('chr' => $row[0], 'start' => (int)$row[3], 'end' => (int)$row[4],
						'score' => (int)$row[5], 'type' => $row[2], 'sample' => strstr(basename($file), ".", true));

		
		//extract fields and convert column 8 
		$row[8] = trim($row[8], ';');
		$attrs = explode(";", $row[8]);
		
		for($ii = 0; $ii<count($attrs); ++$ii)
		{
			$attrs[$ii] = trim($attrs[$ii]);
		}
		
		$function = array();
		$gene = array();
		$tmp = true;
		for($ii = 0; $ii<count($attrs); ++$ii)
		{
			$attr = trim($attrs[$ii]);
			list($label, $value) = explode("=", $attr);
			
			if($label == 'name') 
			{
				$rest  = array_slice($attrs, $ii);
	
				for($ii; $ii < count($attrs); ++$ii)
				{
					$attr = trim($attrs[$ii]);
					list($label, $value) = explode("=", $attr);
					
					if($label == 'function') 
					{
						$function[$value] = true;
					}
			
					if($label == 'name') 
					{
						$gene[$value] = true;
					}
				}
				
				if(count($function) > 0) 
				{
					$fields['function'] = implode(",", array_keys($function));
				}
				
				if(count($gene) > 0) 
				{
					$fields['gene'] = implode(",", array_keys($gene));
				}
				
				$fields['rest'] = implode(";", $rest).";";			
				break;
			}

			$fields[$label] = $value;
		}
		
		//resort, reformat and generate output
		$output = array();
		foreach($headers as $column => $label)
		{
			$value = null;
			if(array_key_exists($column, $fields))
			{
				$value = $fields[$column];
			}
						
			$output[$label] = $value;
		}

		
		$output['genotype'] = trim($output['genotype'], 'oe');
		$output['genotype'] = lcfirst($output['genotype']);
		
		if($row[2]=="SNP")
		{
			if($output['genotype'] == "het")
			{
				$output['obs'] = preg_replace('/'.$output['ref'].'/', '', $output['obs']);
				$output['obs'] = trim($output['obs'], '/');
			}
			if($output['genotype'] == "hom")
			{
				$output['obs'] = substr($output['obs'], -1);
			}
		}
		
		if($row[2]=="Indel")
		{
			$output['state'] = $fields['status'];
			$output['ref'] = '-';
			$output['obs'] = $fields['Base'];
				
			if($fields['indelType'] < 0)
			{
				$output['ref'] = $fields['Base'];
				$output['obs'] = '-';
			}
		}
		
		//convert to tsv
		fwrite($handle1, implode("\t", $output)."\n");
	}
	gzclose($handle2);
}
fclose($handle1);