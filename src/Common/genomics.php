<?php

require_once("functions.php");
require_once("Matrix.php");
require_once("db.php");

/**
	@brief Corrects INDEL variants with surrounding reference bases.
	@ingroup genomics
*/
function correct_indel($start, $ref, $obs)
{
    //trigger error if the variant is already normalized (positions of indels would be messed up)
    if($ref=="-" || $ref=="" || $obs=="-" || $obs=="")
    {
		trigger_error("Cannot correct position of already normalized variant $start:$ref>$obs!", E_USER_ERROR);
    }
	
    $ref = strtoupper($ref);
    $obs = strtoupper($obs);
    
    //remove common first base
    if($ref!="" && $obs!="" && $ref[0]==$obs[0])
    {
		$ref = substr($ref,1);
		$obs = substr($obs,1);
		$start +=1;
    }
    
    //remove common suffix
    $suff = strlen(common_suffix($ref, $obs));
    if ($suff>0)
    {
		$ref = substr($ref,0,-$suff);
		$obs = substr($obs,0,-$suff);
    }

    //remove common prefix
    $pref = strlen(common_prefix($ref, $obs));
    if ($pref>0)
    {
		$ref = substr($ref,$pref);
		$obs = substr($obs,$pref);
		$start += $pref;
    }
 
    //determine start and end
    $end = $start;        
    $ref_c = strlen($ref);
    $obs_c = strlen($obs);
	
    if ($obs_c==1 && $ref_c==1) //SNV
	{
		//nothing to do for SNVs
	}
	else if($ref_c == 0) //insertion
    {
        $ref="-";
		
		//change insertions from before the coordinate to are after the coordinate!!!!
		$start -= 1; 
		$end -= 1;
    }
    else if($obs_c == 0) //deletion
    {
        $end = $start + $ref_c -1;        
        $obs="-";
    }
    else if($obs_c>=1 && $ref_c>1) //complex indel
    { 
        $end = $start + $ref_c -1;       
    }
    
    return array($start, $end, $ref, $obs);
}

/**
	@brief Returns the reverse complementary sequence of a given sequence.
	@ingroup genomics
*/
function rev_comp($input)
{
	//check for invalid chars
	if (!preg_match("/^[acgtACGT]*$/", $input))
	{
		trigger_error("The input sequence '$input' contains characters other than 'acgtACGT'.", E_USER_ERROR);
	}
	
	return strtr(strrev($input), array("A"=>"T", "T"=>"A", "G"=>"C", "C"=>"G", "a"=>"t", "t"=>"a", "c"=>"g", "g"=>"c"));
}

/**
	@brief Returns a genomic reference sequence (1-based chromosomal coordinates).	
	@ingroup genomics
*/
function get_ref_seq($chr, $start, $end)
{
	//fix chromosome for GHCh37
	if ($chr=="chrM") $chr = "chrMT";
	
	//get sequence
	$output = array();
	exec(get_path("samtools")." faidx ".get_path("local_data")."/GRCh37.fa $chr:{$start}-$end 2>&1", $output, $ret);
	if ($ret!=0)
	{
		trigger_error("Error in get_ref_seq: ".implode("\n", $output), E_USER_ERROR);
	}
	
	return implode("", array_slice($output, 1));
}

/**
	@brief Returns the two bases (array) from IUPAC one base code.
	@ingroup genomics
*/
function from_IUPAC($in)
{
	$nuc1 = "none";
	$nuc2 = "none";
	//check for invalid chars
	if (strlen($in)!=1 | !preg_match("/^[ACTGKMRYSW]$/i", $in))
	{
		trigger_error("Genotype '$in' not valid for IUPAC-code. Should be within 'ACTGKMRYSW'.", E_USER_ERROR);
	}
		
	if($in == "K") {$nuc1 = "G"; $nuc2 = "T";}
	elseif($in == "M" or $in == "m") {$nuc1 = "A"; $nuc2 = "C";}
	elseif($in == "R" or $in == "r") {$nuc1 = "A"; $nuc2 = "G";}
	elseif($in == "Y" or $in == "y") {$nuc1 = "C"; $nuc2 = "T";}
	elseif($in == "S" or $in == "s") {$nuc1 = "C"; $nuc2 = "G";}
	elseif($in == "W" or $in == "w") {$nuc1 = "A"; $nuc2 = "T";}
	elseif($in == "A" or $in == "a") {$nuc1 = "A"; $nuc2 = "A";}
	elseif($in == "C" or $in == "c") {$nuc1 = "C"; $nuc2 = "C";}
	elseif($in == "T" or $in == "t") {$nuc1 = "T"; $nuc2 = "T";}
	elseif($in == "G" or $in == "g") {$nuc1 = "G"; $nuc2 = "G";}
         	
	return array($nuc1, $nuc2);
}

/**
	@brief Sanitizes the a chromosome string.
	
	The following changes are made:
	- strips the leading 'chr' if it is present.
	- changes M, X and Y to upper-case.
	
	@return The sanitized chromosome string.
	@ingroup genomics
*/
function chr_trim($chr)
{
	$chr = strtoupper($chr);
	if (strlen($chr)>3 && $chr[0]=="C" && $chr[1]=="H" && $chr[2]=="R")
	{
		$chr = substr($chr, 3);
	}
	return $chr;
}

/**
	@brief Checks if a chromosome string is valid: 1, 2, ..., $max, X, Y, M.
	
	@return The sanitized chromosome string.
	@ingroup genomics
*/
function chr_check($chr, $max = 22, $fail_trigger_error = true)
{
	$chr = chr_trim($chr);
	
	if($chr!="X" && $chr!="Y" && $chr!="M" && $chr!="MT" && (!ctype_digit($chr) || $chr<1 || $chr>$max))
	{
		if ($fail_trigger_error)
		{
			trigger_error("The input string '$chr' is not a valid chromosome!", E_USER_ERROR);
		}
		else
		{
			return false;
		}
	}
	
	return $chr;
}


/**
	@brief Returns the base count of chromosomes for different builds.
	@ingroup genomics
*/
function chr_info($chr, $build = "hg19")
{
	if($build=="hg19")
	{
		if ($chr=="all")  return 3095693983;
		$chr = chr_check($chr);
		if ($chr=="M") return 16571;
		if ($chr=="1") return 249250621;
		if ($chr=="2") return 243199373;
		if ($chr=="3") return 198022430;
		if ($chr=="4") return 191154276;
		if ($chr=="5") return 180915260;
		if ($chr=="6") return 171115067;
		if ($chr=="7") return 159138663;
		if ($chr=="8") return 146364022;
		if ($chr=="9") return 141213431;
		if ($chr=="10") return 135534747;
		if ($chr=="11") return 135006516;
		if ($chr=="12") return 133851895;
		if ($chr=="13") return 115169878;
		if ($chr=="14") return 107349540;
		if ($chr=="15") return 102531392;
		if ($chr=="16") return 90354753;
		if ($chr=="17") return 81195210;
		if ($chr=="18") return 78077248;
		if ($chr=="19") return 59128983;
		if ($chr=="20") return 63025520;
		if ($chr=="21") return 48129895;
		if ($chr=="22") return 51304566;
		if ($chr=="X") return 155270560;
		if ($chr=="Y") return 59373566;
	}
	else if($build=="mm9")
	{
		if ($chr=="all")  return 2654911517;
		$chr = chr_check($chr, 19);
		if ($chr=="M") return 16299;
		if ($chr=="1") return 197195432;
		if ($chr=="2") return 181748087;
		if ($chr=="3") return 159599783;
		if ($chr=="4") return 155630120;
		if ($chr=="5") return 152537259;
		if ($chr=="6") return 149517037;
		if ($chr=="7") return 152524553;
		if ($chr=="8") return 131738871;
		if ($chr=="9") return 124076172;
		if ($chr=="10") return 129993255;
		if ($chr=="11") return 121843856;
		if ($chr=="12") return 121257530;
		if ($chr=="13") return 120284312;
		if ($chr=="14") return 125194864;
		if ($chr=="15") return 103494974;
		if ($chr=="16") return 98319150;
		if ($chr=="17") return 95272651;
		if ($chr=="18") return 90772031;
		if ($chr=="19") return 61342430;
		if ($chr=="X") return 166650296;
		if ($chr=="Y") return 15902555;
	}
	else
	{
		trigger_error("chr_info: Unknown build '$build'.", E_USER_ERROR);
	}

	return 0;
}


/**
	@brief Returns the list of all chromosomes.
	@ingroup genomics
*/
function chr_list()
{
	return array_merge(range(1,22), array("X","Y"));
}

/**
	@brief Function to get central organized paths to tools.
	
	@param name name of variable in ini-file.
	@return value of variable in ini-file (e.g. path).
		
	@ingroup helpers
*/
function get_path($name, $throw_on_error=true)
{
	$dir = repository_basedir();
	
	//determine INI file to name
	if (isset($GLOBALS["path_ini"]))
	{
		$ini_file = $GLOBALS["path_ini"];
	}
	else
	{
		if (file_exists($dir."settings.ini"))
		{
			$ini_file = $dir."settings.ini";
		}
		else
		{
			$ini_file = $dir."settings.ini.default";
		}
	}
		
	//parse ini file
	$parsed_ini = parse_ini_file($ini_file);	
	if($parsed_ini===FALSE)
	{
		trigger_error("Could not parse INI file '$ini_file'.",E_USER_ERROR);
	}
	
	//get value
	if (!isset($parsed_ini[$name]) && $throw_on_error)
	{
		trigger_error("Could not find key '$name' in settings file '$ini_file'!", E_USER_ERROR);
	}
	@$value = $parsed_ini[$name];

	//replace [path] by base path
	$value = str_replace("[path]", $dir, $value);
	
	return $value;
}


/**
	@brief Function to get central organized db credentials.
	
	@param database name of database in ini-file.
	@param name name of variable in ini-file.
	@return value of variable in ini-file (e.g. path).
		
	@ingroup helpers
*/
function get_db($db, $name)
{
	$values = get_path($name);
	
	if (!isset($values[$db]))
	{
		trigger_error("get_db could not find value '$name' for DB '$db'!", E_USER_ERROR);
	}
	
	return $values[$db];
}

/**
	@brief Returns the list of all databases from the settings INI file.
	
	@ingroup helpers
*/
function db_names()
{
	return array_keys(get_path('db_name'));
}

/**
	@brief Returns if a database is enabled (all properties set).
	
	@ingroup helpers
*/
function db_is_enabled($name)
{
	$db_properties = array('db_host', 'db_name', 'db_user', 'db_pass');
	foreach($db_properties as $db_prop)
	{
		$entries = get_path($db_prop);
		if (!isset($entries[$name]) || trim($entries[$name])=="")
		{
			return false;
		}
	}
	
	return true;
}

///Loads a processing system INI file. If the file name is empty, the system is determine from the processed sample name, written to a temporary file and the filename is set to that temporary file.
function load_system(&$filename, $ps_name = "")
{	
	//determine system from processed sample name
	if (is_null($filename) || $filename=="")
	{
		//get ID of processed sample
		$db_conn = DB::getInstance("NGSD");
		$ps_id = get_processed_sample_id($db_conn, $ps_name, false);
		if ($ps_id==-1)
		{
			trigger_error("load_system: Cannot determine processing system - processed sample name '$ps_name' is invalid!", E_USER_ERROR);
		}
		
		//store system INI file
		$sys_name = $db_conn->getValue("SELECT sys.name_short FROM processing_system sys, processed_sample ps WHERE ps.processing_system_id=sys.id AND ps.id=$ps_id");
		$filename = temp_file(".ini", "pro_sys_".$sys_name."_");
		store_system($db_conn, $sys_name, $filename);
	}
	
	return parse_ini_file($filename);
}

function store_system(&$db_conn, $name_short, $filename)
{
	//execute query
	$res = $db_conn->executeQuery("SELECT sys.*, g.build FROM processing_system sys, genome g WHERE sys.genome_id=g.id and sys.name_short=:name_short", array("name_short"=>$name_short));
	if(count($res)==0)
	{
		trigger_error("Processing system with short name '$name_short' not found in NGSD!", E_USER_ERROR);
	}
	
	//store output
	$output = array();
	$output[] = "name_short = \"".$name_short."\"";
	$output[] = "name_manufacturer = \"".$res[0]['name_manufacturer']."\"";
	$output[] = "target_file = \"".$res[0]['target_file']."\"";
	$output[] = "adapter1_p5 = \"".$res[0]['adapter1_p5']."\"";
	$output[] = "adapter2_p7 = \"".$res[0]['adapter2_p7']."\"";
	$output[] = "shotgun = ".$res[0]['shotgun'];
	$output[] = "umi_type = \"".$res[0]['umi_type']."\"";
	$output[] = "type = \"".$res[0]['type']."\"";
	$output[] = "build = \"".$res[0]['build']."\"";
	$output[] = "";
	file_put_contents($filename, implode("\n", $output));
}

/**
	@brief Loads the qcML terms for NGS from the OBO file.
	@ingroup helpers
*/
function load_qc_terms()
{
	//do nothing if already loaded
	if (isset($GLOBALS["qcml"])) return;
	
	//load terms
	$terms = array();
	
	$current = array();
	$h = fopen(repository_basedir()."/data/misc/qc-cv.obo", "r");
	while(!feof($h))
	{
		$line = trim(fgets($h));
		if ($line=="")
		{
			continue;
		}
		else if ($line=="[Term]")
		{
			if (isset($terms[$current['id']])) trigger_error("duplicate qcML term id {$current['id']}!", E_USER_ERROR);
			$terms[$current['id']] = $current;
			$current = array();
		}
		else if (starts_with($line, "id:"))
		{
			$current['id'] = trim(substr($line, 3));
		}
		else if (starts_with($line, "name:"))
		{
			$current['name'] = trim(substr($line, 5));
		}
		else if (starts_with($line, "def:"))
		{
			$parts = explode("\"", $line);
			$current['def'] = trim($parts[1]);
		}
		else if (starts_with($line, "xref: value-type:xsd\:"))
		{
			$parts = explode(":", strtr($line, "\"", ":"));
			$current['type'] = trim($parts[3]);
		}
		else if (starts_with($line, "comment:"))
		{
			//nothing to do here
		}
		else if(starts_with($line, "is_obsolete:"))
		{
			//nothing to do here
		}
	}
	if (isset($terms[$current['id']])) trigger_error("duplicate qcML term id '{$current['id']}'!", E_USER_ERROR);
	$terms[$current['id']] = $current;
	
	//remove QC terms we do not need
	foreach($terms as $id => $data)
	{
		//terms that are not in the NGS namespace
		if (!starts_with($id, "QC:2"))
		{
			unset($terms[$id]);
		}
		//no type (parent terms needed to structure the ontology only)
		if (!isset($data['type']))
		{
			unset($terms[$id]);
		}
	}
	ksort($terms);
	$GLOBALS["qcml"] = $terms;
}

/**
	@brief Determines the gender based on a list of genotypes.
	
	@param genotypes The array of genotypes
	@param het The genotype string used for heterocygote.
	@param male Cutoff for male (below this fraction).
	@param female Cutoff for female (above this fraction).
	@return Returns an array with gender ('m' or 'f') and herocygote ratio. Or FALSE (if undetermined).
	
	@ingroup genomics
*/
function gender($genotypes, $het, $male, $female)
{
	$counts = array_count_values($genotypes);
	
	if (!isset($counts[$het]))
	{
		$counts[$het] = 0;
	}

	$ratio = $counts[$het] / count($genotypes);
	
	if ($ratio > $female)
	{
		return array('f', $ratio);
	}
	if ($ratio < $male)
	{
		return array('m', $ratio);
	}
	
	// unsure
	return array(false, $ratio);
}

/**
	@brief Returns the processed sample database ID for a processed sample name, or -1 if the processed sample is not in the database.
 */
function get_processed_sample_id(&$db_conn, $name, $error_if_not_found=true)
{
	//split name
	$parts = explode("_", $name."_99");
	list($sname, $id) = $parts;
	$id = ltrim($id, "0");
	
	//query NGSD
	try 
	{
		$res = $db_conn->executeQuery("SELECT ps.id FROM processed_sample ps, sample s WHERE ps.sample_id=s.id AND s.name=:name AND ps.process_id=:id", array('name' => $sname, "id"=>$id));
	}
	catch(PDOException $e)
	{
		if ($error_if_not_found) throw $e;
		return -1;
	}
	
	//processed sample not found
	if (count($res)!=1)
	{
		if ($error_if_not_found) trigger_error("Could not find processed sample with name '$name' in NGSD!", E_USER_ERROR);
		return -1;
	}
	
	return $res[0]['id'];
}

///Updates normal sample entry for given tumor sample.
function updateNormalSample(&$db_conn, $ps_tumor, $ps_normal, $overwrite = false)
{
	//TUMOR
	//get processed sample ID
	list($s, $p) = explode("_", $ps_tumor);
	$p = (int)$p;
	$res = $db_conn->executeQuery("SELECT ps.id FROM processed_sample ps, sample s WHERE ps.sample_id=s.id AND s.name='$s' and ps.process_id='$p'");
	$ps_tid = $res[0]['id'];
	// get normal_id
	$res = $db_conn->executeQuery("SELECT normal_id FROM processed_sample ps, sample s WHERE ps.sample_id=s.id AND s.name='$s' and ps.process_id='$p'");
	$n = $res[0]['normal_id'];

	//NORMAL
	//get processed sample ID
	list($s, $p) = explode("_", $ps_normal);
	$p = (int)$p;
	$res = $db_conn->executeQuery("SELECT ps.id FROM processed_sample ps, sample s WHERE ps.sample_id=s.id AND s.name='$s' and ps.process_id='$p'");
	$ps_nid = $res[0]['id'];

	if(empty($ps_nid) || empty($ps_tid))	trigger_error("Could not find either tumor ($ps_tumor) or normal ($ps_normal) in NGSD.");

	if(!empty($n) && $ps_nid!=$n)	trigger_error("Different normal sample found in NGSD (NGSD-IDs $ps_nid!=$n) for tumor ($ps_tumor).",E_USER_WARNING);

	if(empty($n) || $overwrite)
	{
		$db_conn->executeStmt("UPDATE processed_sample SET normal_id='$ps_nid' WHERE id='$ps_tid'");
		return true;
	}

	return false;
}

//@TODO implement *-syntax for HGVS compliance (end of coding)
function convert_hgvs2genomic($transcript, $cdna, $error = true)
{
	//extract cDNA position and
	$chr = null;
	$start = null;
	$offset1 = 0;
	$end = null;
	$offset2 = 0;
	$strand = null;
	$ref = null;
	$obs = null;
	$matches = array();	//for preg_match
	$e = null;
	if(preg_match("/^c\.(?<start>\d+)(?<offset1>[\-\+]\d+)?(?<ref>[ACGT])[\>\<](?<obs>[ACGT])$/i",$cdna,$matches)!=0)	//SNV
	{
		$result = convert_coding2genomic($transcript, $matches["start"], $matches["start"],$error);
		if(is_array($result))	list($chr,$start,$end,$strand) = $result;
		else	$e = $result;
		if(!empty($matches["offset1"]))	$offset1 = $matches["offset1"];
		$offset2 = $offset1;
		$ref = $matches["ref"];
		$obs = $matches["obs"];			
	}
	else if(preg_match("/^c\.(?<start>\d+)(?<offset1>[\-\+]\d+)?_?(?<end>\d+)?(?<offset2>[\-\+]\d+)?del(?<ref>[CATG]+)?$/i",$cdna,$matches)!=0)	//Deletion
	{
		if(empty($matches["end"]))	$matches["end"] = $matches["start"];	//if no end position is given
		
		$result = convert_coding2genomic($transcript, $matches["start"], $matches["end"],$error);
		if(is_array($result))	list($chr,$start,$end,$strand) = $result;
		else	$e = $result;
		
		if(!empty($matches["offset1"]))	$offset1 = $matches["offset1"];
		$offset2 = $offset1;
		if(!empty($matches["offset2"]))	$offset2 = $matches["offset2"];
		if(!empty($matches["ref"]))	$ref = $matches["ref"];
		$obs = "-";
	}
	else if(preg_match("/^c\.(?<start>\d+)(?<offset1>[\-\+]\d+)?_?(?<end>\d+)?(?<offset2>[\-\+]\d+)?del(?<ref_count>\d+)?$/i",$cdna,$matches)!=0)	//Deletion, e.g. c.644-12del16
	{
		if(empty($matches["end"]))	$matches["end"] = $matches["start"];	//if no end position is given
		
		$result = convert_coding2genomic($transcript, $matches["start"], $matches["end"],$error);
		if(is_array($result))	list($chr,$start,$end,$strand) = $result;
		else	$e = $result;

		if(!empty($matches["offset1"]))	$offset1 = $matches["offset1"];
		$offset2 = $offset1;
		if(!empty($matches["ref_count"]))	$offset2 += $matches["ref_count"]-1;	//if no end position is given		
		$obs = "-";
	}
	else if(preg_match("/^c\.(?<start>\d+)(?<offset1>[\-\+]\d+)?_?(?<end>\d+)?(?<offset2>[\-\+]\d+)?ins(?<obs>[CATG]+)$/i",$cdna,$matches)!=0)	//Insertion
	{
		//skip end and offset2, since insertion is always next to start (both splicing and coding)
		$result = convert_coding2genomic($transcript, $matches["start"], $matches["start"],$error);
		if(is_array($result))	list($chr,$start,$end,$strand) = $result;
		else	$e = $result;

		//offsets
		if(!empty($matches["offset1"]))	$offset1 = $matches["offset1"];
		if(!empty($matches["offset2"]))	$offset2 = $matches["offset2"];
		if($strand=="+" && $offset1!=0 && $offset2!=0)	$offset1 = min($offset1, $offset2);
		if($strand=="-" && $offset1!=0 && $offset2!=0)	$offset1 = max($offset1, $offset2);
		$offset2 = $offset1;
		if($strand=="-" && empty($offset1) && empty($offset2))	$end = --$start;	//change of insertion site required for "-"-strand variants.
		
		//alleles
		$ref = "-";
		$obs = $matches["obs"];
	}
	else if(preg_match("/^c\.(?<start>\d+)(?<offset1>[\-\+]\d+)?_?(?<end>\d+)?(?<offset2>[\-\+]\d+)?del(?<ref>[CATG]+)?ins(?<obs>[CATG]+)$/i",$cdna,$matches)!=0)	//combined InDel
	{
		if(empty($matches["end"]))	$matches["end"] = $matches["start"];	//if no end position is given

		$result = convert_coding2genomic($transcript, $matches["start"], $matches["end"],$error);
		if(is_array($result))	list($chr,$start,$end,$strand) = $result;
		else	$e = $result;

		
		if(!empty($matches["offset1"]))	$offset1 = $matches["offset1"];
		if(!empty($matches["offset2"]))	$offset2 = $matches["offset2"];
		if(!empty($matches["ref"]))	$ref = $matches["ref"];
		if(empty($ref))	$ref = get_ref_seq($chr,$start,$end);
		if($strand=="-")	$ref = rev_comp ($ref);
		$obs = $matches["obs"];			
	}
	else if(preg_match("/^c\.(?<start>\d+)(?<offset1>[\-\+]\d+)?_?(?<end>\d+)?(?<offset2>[\-\+]\d+)?dup(?<obs>[CATG]+)?$/i",$cdna,$matches)!=0)	//Duplication
	{
		if(empty($matches["end"]))	$matches["end"] = $matches["start"];
		
		$result = convert_coding2genomic($transcript, $matches["start"], $matches["end"],$error);
		if(is_array($result))	list($chr,$start,$end,$strand) = $result;
		else	$e = $result;

		if($strand == "+")	$end = --$start;
		if($strand == "-")	$start = $end;
		//if on - strand move insertion to the right
		if(!empty($matches["offset1"]))	$offset1 = $matches["offset1"];
		if(!empty($matches["offset2"]))	$offset2 = $matches["offset2"];
		$ref = "-";
		$obs = get_ref_seq($chr,$start,$end);
		if(!empty($matches["obs"]))	$obs = $matches["obs"];			
		if(strlen($obs)==1)	$start=--$end;
	}
	else	//default (not identifiable)
	{
		if($error)	trigger_error("Could not identify HGVS for variant: $cdna.",E_USER_ERROR);
		return "Could not identify variant HGVS for variant: $cdna.";
	}

	if(!is_null($e))
	{
		return $e;
	}
	
	if($strand=="+")
	{
		$start += $offset1;
		$end += $offset2;
		if($obs=="-" && empty($ref))	$ref = get_ref_seq($chr,$start,$end);
		$ref = strtoupper($ref);
		$obs = strtoupper($obs);
	}
	if($strand == "-")	
	{
		$start -= $offset2;
		$end -= $offset1;
		
		//convert reference
		if($obs=="-" && empty($ref))	$ref = strtoupper(get_ref_seq($chr,$start,$end));
		else if($ref!="-")	$ref = strtoupper(rev_comp($ref));
		
		//convert obs
		if($obs!="-")	$obs = strtoupper(rev_comp($obs));
	}

	//check if reference is valid
	$r = get_ref_seq($chr,$start,$end);	//adopt for different builds
	if(!empty($chr) && !empty($ref) && $ref!="-" && strtoupper($r)!=strtoupper($ref))
	{
		if($error)	trigger_error("Wrong reference sequence for HGVS '$transcript:$cdna': is '$ref', should be '".$r."' ($chr:$start-$end).",E_USER_ERROR);
		return "Wrong reference sequence for HGVS '$transcript:$cdna': is '$ref', should be '".$r."' ($chr:$start-$end).";
	}

	//check
	$l = $end - $start + 1;
	$b = strlen($ref);
	if($l!=$b)
	{
		if($error)	trigger_error("HGVS ref does not match lenght of variant '$transcript:$cdna': $chr:$start-$end, ref is '$ref', obs is '$obs'.",E_USER_ERROR);
		return "HGVS ref does not match lenght of variant '$transcript:$cdna': $chr:$start-$end, ref is '$ref', obs is '$obs'.";
	}
	
	return array($chr,$start,$end,$ref,$obs);
}

function convert_coding2genomic($transcript,$cdna_start,$cdna_end, $error = true)
{
	//identify transcript ID, currently only refseq allowed
	if(!preg_match('/(?<id>NM_\d+)(\.\d+)?/',$transcript,$matches))
	{
		if($error)	trigger_error("Invalid coding transcript ID '$transcript' (currently only refseq supported).",E_USER_ERROR);
		return "Invalid coding transcript ID '$transcript' (currently only refseq supported).";
	}
	$refseq_id = $matches["id"];

	//get gene information, chromosome, strand
	$exons = array();
	$chr = null;
	$strand = null;
	$transcript_length = 0;
	$known_gene = get_path("data_folder")."/dbs/UCSC/refGene.txt";
	$handle = fopen($known_gene, "r");
	if($handle)
	{
		while(($buffer=fgets($handle)) !== FALSE)
		{
			$row = explode("\t", $buffer);
			if($row[1]==$refseq_id)
			{
				$chr = $row[2];
				$strand = $row[3];
				$cdsStart = $row[6];
				$cdsEnd = $row[7];
				$exons_start = explode(",", $row[9]);
				$exons_end = explode(",", $row[10]);
				for($i=0;$i<count($exons_start);++$i)
				{
					if(empty($exons_start[$i]))	continue;
					$exons[] = array('start'=>$exons_start[$i], 'end'=>$exons_end[$i]);
					$transcript_length += ($exons_end[$i]-$exons_start[$i]);
				}
				break;
			}
		}
		fclose($handle);
	}
	else
	{
		if($error)	trigger_error("Could not open file $known_gene.",E_USER_ERROR);
		return  "Could not open file $known_gene.";
	}
	if(count($exons)==0)
	{
		if($error)	trigger_error("Could not find exons for $ucsc_id.",E_USER_ERROR);
		return "Could not find exons for $ucsc_id.";
	}
		
	//get genomic positions for coding cDNA position, all coordinates are 0-based (knownGene.txt)
	$start = null;
	$end = null;
	$coding_start = false;
	$coding_end = false;
	$first_coding_exon = null;
	$coding_length = 0;
	$tmp_basepairs_start = $cdna_start;
	$tmp_basepairs_end = $cdna_end;
	if($strand == "+")
	{
		for($i=0;$i<count($exons);++$i)
		{
			//get coding length, adopt if translation start or end is within this exon
			$tmp_start = $exons[$i]["start"];
			$tmp_end = $exons[$i]["end"];
			if($cdsStart>=$tmp_start && $cdsStart<$tmp_end)	//find translation start
			{
				$tmp_start = $cdsStart;
				$coding_start = true;
			}
			if($cdsEnd>$tmp_start && $cdsEnd<=$tmp_end)	//find translation start
			{
				$tmp_end = $cdsEnd;
				$coding_end = true;
			}
			$length = $tmp_end - $tmp_start;
			
			//subtract coding basepairs of this exon from given cDNA positions
			if($coding_start)
			{
				$tmp_basepairs_start-=$length;
				$tmp_basepairs_end-=$length;
				$coding_length+=$length;
			}

			//identify cDNA start and end position (no more tmp_basepairs left)
			if($tmp_basepairs_start<=0 && $start==null)	$start = $tmp_end + $tmp_basepairs_start;
			if($tmp_basepairs_end<=0 && $end==null)	$end = $tmp_end + $tmp_basepairs_end;
			if($start!=null && $end!=null)	break;	//no need to check additional exons since coding start and end were already found
		}
	}
	else if($strand == "-")
	{
		$exons = array_reverse($exons);
		for($i=0;$i<count($exons);++$i)
		{
			//get coding length, adopt if translation start or end is within this exon
			$tmp_start = $exons[$i]["start"];
			$tmp_end = $exons[$i]["end"];
			if($cdsEnd>$exons[$i]["start"] && $cdsEnd<=$exons[$i]["end"])	//find translation start (= cdsEnd in reverse mode)
			{
				$tmp_end = $cdsEnd;
				$coding_start = true;
			}
			if($cdsStart>=$exons[$i]["start"] && $cdsStart<$exons[$i]["end"])	//find translation end (= cdsStart in reverse mode)
			{
				$tmp_start = $cdsStart;
				$coding_end = true;
			}
			$length = $tmp_end - $tmp_start;
		
			//subtract coding basepairs of this exon from given cDNA positions
			if($coding_start)
			{
				$tmp_basepairs_start-=$length;
				$tmp_basepairs_end-=$length;
				$coding_length+=$length;
			}
			
			//identify cDNA start and end position (no more tmp_basepairs left)
			if($tmp_basepairs_start<=0 && $end==null)	$end = $tmp_start - $tmp_basepairs_start + 1;	//convert 0-based start to 1-based start, convert strand
			if($tmp_basepairs_end<=0 && $start==null)	$start = $tmp_start - $tmp_basepairs_end + 1;	//convert 0-based start to 1-based start, convert strand
			if($start!=null && $end!=null)	break;
		}
	}
	if($start<$cdsStart || $end>=($cdsEnd+1))	//0-based the other way round...
	{
		if($error)	trigger_error("Given cDNA position is invalid (given cDNA start: $cdna_start, given cDNA end: $cdna_end; length of coding transcript: $coding_length bp).",E_USER_ERROR);
		return "Given cDNA position is invalid (given cDNA start: $cdna_start, given cDNA end: $cdna_end; length of coding transcript: $coding_length bp).";
	}

	return array($chr,$start,$end,$strand);
}

function indel_for_vcf($chr, $start, $ref, $obs)
{
	$ref = trim($ref,"-");
	$obs = trim($obs,"-");

	// handle indels
	$extra1 = "";
	
	// 1. simple insertion - add prefix
	if(strlen($ref)==0 && strlen($obs)>=1)
	{
		$extra1 = get_ref_seq($chr,$start,$start);
	}

	// 2. simple deletion - correct start and add prefix
	if(strlen($ref)>=1 && strlen($obs)==0)
	{
		$start -= 1;
		$extra1 = get_ref_seq($chr,$start,$start);
	}

	// 3. complex indel, nothing to do
	
	// combine all information
	$ref = strtoupper($extra1.$ref);
	$obs = strtoupper($extra1.$obs);
	
	return array($chr,$start,$ref,$obs);
}

function is_valid_ref_sample_for_cnv_analysis($file, $tumor_only = false)
{
	//check that sample is not NIST reference sample (it is a cell-line)
	if (contains($file, "NA12878")) return false;

	//no NGSD => error
	if (!db_is_enabled("NGSD"))
	{
		trigger_error("is_valid_ref_sample_for_cnv_analysis needs NGSD access!", E_USER_ERROR);
	}

	//check sample is in NGSD
	$db_conn = DB::getInstance("NGSD");
	$ps_id = get_processed_sample_id($db_conn, $file, false);
	if ($ps_id<0) return false;
	
	//check that sample is not tumor and not not FFPE
	$res = $db_conn->executeQuery("SELECT s.tumor, s.ffpe, ps.quality q1, r.quality q2, p.type FROM sequencing_run r, sample s, processed_sample ps, project p WHERE s.id=ps.sample_id AND ps.id='$ps_id' AND ps.sequencing_run_id=r.id AND ps.project_id = p.id");
	if ($tumor_only)
	{
		if ($res[0]['tumor']!="1") return false;
	}
	else
	{
		if ($res[0]['tumor']=="1") return false;
		if ($res[0]['ffpe']=="1") return false;
	}
	
	//check that run and processed sample do not have bad quality
	if ($res[0]['q1']=="bad") return false;
	if ($res[0]['q2']=="bad") return false;
	
	//check that project type is research/diagnostics
	if ($res[0]['type']!="research" && $res[0]['type']!="diagnostic") return false;
	
	return true;
}

function is_valid_ref_tumor_sample_for_cnv_analysis($file, $discard_ffpe = false)
{
	//check that sample is not NIST reference sample (it is a cell-line)
	if (contains($file, "NA12878")) return false;

	//no NGSD => error
	if (!db_is_enabled("NGSD"))
	{
		trigger_error("is_valid_ref_tumor_sample_for_cnv_analysis needs NGSD access!", E_USER_ERROR);
	}
	
	//check sample is in NGSD
	$db_conn = DB::getInstance("NGSD");
	$ps_id = get_processed_sample_id($db_conn, $file, false);
	
	//check that sample is not FFPE
	$res = $db_conn->executeQuery("SELECT s.tumor, s.ffpe, ps.quality q1, r.quality q2, p.type FROM sequencing_run r, sample s, processed_sample ps, project p WHERE s.id=ps.sample_id AND ps.id='$ps_id' AND ps.sequencing_run_id=r.id AND ps.project_id = p.id");
	if ($ps_id<0) return false;
	if ($discard_ffpe && $res[0]['ffpe']=="1") return false;
	
	//check that run and processed sample do not have bad quality
	if ($res[0]['q1']=="bad") return false;
	if ($res[0]['q2']=="bad") return false;
	
	//check that project type is research/diagnostics
	if ($res[0]['type']!="research" && $res[0]['type']!="diagnostic") return false;
	
	//check that sample is tumor sample
	if ($res[0]['tumor'] != "1") return false;
	
	return true;
}

///Loads a VCF file in a normalized manner, e.g. for comparing them
/// - sorts header lines
/// - sorts info column by name
/// - sorts and format/sample columns by format string
function load_vcf_normalized($filename)
{
	$comments = array();
	$header = "";
	$vars = array();
	
	//load and normalize data
	if(!is_file($filename))	trigger_error("Could not find file $filename.",E_USER_WARNING);
	$file = file($filename);
	foreach($file as $line)
	{
		$line = nl_trim($line);
		if ($line=="") continue;
		
		if (starts_with($line, "##"))
		{
			$comments[] = $line;
		}
		else if (starts_with($line, "#"))
		{
			$header = $line;
		}
		else
		{
			$parts = explode("\t", $line);
			if (count($parts)<10) trigger_error("VCF file $filename has line with less than 10 colums: $line", E_USER_ERROR);
			list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format) = $parts;
			
			
			//convert INFO data to associative array
			$info = explode(";", $info);
			$tmp = array();
			foreach($info as $entry)
			{
				if (!contains($entry, "="))
				{
					$tmp[$entry] = "";
				}
				else
				{
					list($key, $value) = explode("=", $entry, 2);
					$tmp[$key] = $value;
				}
			}
			$info = $tmp;
			ksort($info);
			$var = array($chr, $pos, $id, $ref, $alt, $qual, $filter, $info);
			
			//convert format/sample data (also additional sample columns of multi-sample VCF)
			for($i=9; $i<count($parts); ++$i)
			{
				$sample = array_combine(explode(":", $format), explode(":", $parts[$i]));
				ksort($sample);
				$var[] = $sample;
			}
			
			$vars[] = $var;
		}
	}

	//output: comments
	sort($comments);
	foreach($comments as $line)
	{
		if (!starts_with($line, "##INFO") && !starts_with($line, "##FORMAT")) $output[] = $line;
	}
	foreach($comments as $line)
	{
		if (starts_with($line, "##INFO")) $output[] = $line;
	}
	foreach($comments as $line)
	{
		if (starts_with($line, "##FORMAT")) $output[] = $line;
	}
	
	//output: header
	$output[] = $header;
	
	//output: variants
	foreach($vars as $var)
	{
		list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info) = $var;
		
		//info
		$tmp = array();
		foreach($info as $key => $value)
		{
			$tmp[] = "$key=$value";
		}
		$info = implode(";", $tmp);
		
		//format
		$format = implode(":", array_keys($var[8]));
				
		//sample columns (also additional sample columns of multi-sample VCF)
		$samples = array();
		for($i=8; $i<count($var); ++$i)
		{
			$samples[] = implode(":", array_values($var[$i]));
		}
		
		$output[] = "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\t$format\t".implode("\t", $samples);
	}
	
	return $output;
}


function vcf_strelka_snv($format_col, $sample_col, $obs)
{
	if (!preg_match("/^[acgtACGT]*$/", $obs))	trigger_error("Invalid observed allele (".$format_col." ".$sample_col." ".$obs.").", E_USER_ERROR);

	$f = explode(":", $format_col);
	$index_depth = array_search("DP",$f);
	$index_TU = array_search("TU",$f);
	$index_AU = array_search("AU",$f);
	$index_CU = array_search("CU",$f);
	$index_GU = array_search("GU",$f);
	if($index_depth===FALSE || $index_TU===FALSE || $index_AU===FALSE || $index_CU===FALSE || $index_GU===FALSE) trigger_error("Invalid strelka format: either field DP, TU, AU, CU or GU not available.", E_USER_ERROR);
	
	$entries = explode(":", $sample_col);
	$d = $entries[$index_depth];
	
	list($tu,) = explode(",", $entries[$index_TU]);
	list($au,) = explode(",", $entries[$index_AU]);
	list($cu,) = explode(",", $entries[$index_CU]);
	list($gu,) = explode(",", $entries[$index_GU]);
	$sum = $au + $tu + $cu + $gu;

	$o = 0;
	if($obs == "T") $o = $tu;
	else if($obs == "A") $o = $au;
	else if($obs == "C") $o = $cu;
	else if($obs == "G") $o = $gu;
	else trigger_error("Alternative allele '$obs' unknown (".$format_col." ".$sample_col." ".$obs.").", E_USER_WARNING);	// unknown alleles multiallelic or '.'

	$f = 0;
	if($sum!=0)	$f = number_format($o/$sum,4);

	return array($d,$f);
}

function vcf_strelka_indel($format_col, $sample_col)
{
	$f = explode(":", $format_col);
	$index_depth = array_search("DP",$f);
	$index_TIR = array_search("TIR",$f);
	$index_TAR = array_search("TAR",$f);
	if($index_depth===FALSE || $index_TIR===FALSE || $index_TAR===FALSE) trigger_error("Invalid strelka format: either field DP, TIR or TAR not available.", E_USER_ERROR);
	
	$entries = explode(":", $sample_col);
	$d = $entries[$index_depth];
	list($tir,) = explode(",", $entries[$index_TIR]);
	list($tar,) = explode(",", $entries[$index_TAR]);
	if(!is_numeric($tir) || !is_numeric($tar)) trigger_error("Could not identify numeric depth for strelka indel (".$format_col." ".$sample_col.")", E_USER_ERROR);

	//tir and tar contain strong supportin reads, tor (not considered here) contains weak supportin reads like breakpoints
	//only strong supporting reads are used for calculation of allele fraction
	$f = ($tir+$tar)==0 ? 0.0 : number_format($tir/($tir+$tar), 4);
	
	return array($d,$f);
}

function vcf_freebayes($format_col, $sample_col)
{
	$g = explode(":",$format_col);
	$index_DP = NULL;
	$index_AO = NULL;
	$index_GT = NULL;
	for($i=0;$i<count($g);++$i)
	{
		if($g[$i]=="DP")	$index_DP = $i;
		if($g[$i]=="AO")	$index_AO = $i;
		if($g[$i]=="GT")	$index_GT = $i;
	}

	if(is_null($index_DP) || is_null($index_AO) ||is_null($index_GT))	trigger_error("Invalid freebayes format; either field DP, GT or AO not available.",E_USER_ERROR);	
	
	$s = explode(":",$sample_col);

	// workaround for bug during splitting of multi-allelic variants - allele counts for multiple alleles are kept
	if(strpos($s[$index_AO],",")!==FALSE)
	{		
		$gt = $s[$index_GT];
		
		$sep = "/";
		if(strpos($s[$index_GT],"|")!==FALSE)	$sep = "|";
		
		$idx_al1 = min(explode($sep,$gt));
		$idx_al2 = max(explode($sep,$gt));
		if($idx_al1!=0 || $idx_al2!=1)	trigger_error("Unexpected error. Allele 1 is $idx_al1, Allele 2 is $idx_al2; expected 0 and 1  (".$format_col." ".$sample_col.").",E_USER_ERROR); 
				
		$tmp = 0;
		$tmp = explode(",",$s[$index_AO])[$idx_al2-1];
		$s[$index_AO] = $tmp;
	}
	
	if(!is_numeric($s[$index_AO]))	trigger_error("Invalid alternative allele count (".$format_col." ".$sample_col.").",E_USER_ERROR);	// currently no multiallelic variants supported
	if(!is_numeric($s[$index_DP]))	trigger_error("Could not identify numeric depth (".$format_col." ".$sample_col.").",E_USER_ERROR);

	$d1 = $s[$index_DP];
	$d2 = $s[$index_AO];
	
	$f = null;
	if($d1>0)	$f = number_format($d2/$d1, 4);
	return array($d1,$f);
}

function vcf_iontorrent($format_col, $sample_col, $idx_al)
{
	$g = explode(":",$format_col);
	$index_DP = NULL;
	$index_AF = NULL;
	for($i=0;$i<count($g);++$i)
	{
		if($g[$i]=="DP")	$index_DP = $i;
		if($g[$i]=="AF")	$index_AF = $i;
	}

	if(is_null($index_DP) || is_null($index_AF))	trigger_error("Invalid iontorrent format; either field DP or AF not available.",E_USER_ERROR);
	
	$d = explode(":",$sample_col)[$index_DP];
	$f = number_format(explode(",",explode(":",$sample_col)[$index_AF])[$idx_al], 4);

	return array($d,$f);
}

//Converts genotye 
function vcfgeno2human($gt, $upper_case=false)
{
	$gt = strtr($gt, "/.", "|0");
	if ($gt=="0|0")
	{
		$geno = "wt";
	}
	else if ($gt=="0|1" || $gt=="1|0")
	{
		$geno = "het";
	}
	else if ($gt=="1|1")
	{
		$geno = "hom";
	}
	else
	{
		trigger_error("Invalid VCF genotype '$gt'!", E_USER_ERROR);
	}
	
	return $upper_case ? strtoupper($geno) : $geno;
}

//Returns information from the NGSD about a processed sample (as key-value pairs), or null/error if the sample is not found.
function get_processed_sample_info(&$db_conn, $ps_name, $error_if_not_found=true)
{
	//get info from NGSD
	list($sample_name, $process_id) = explode("_", $ps_name."_");
	$res = $db_conn->executeQuery("SELECT p.name as project_name, p.type as project_type, p.analysis as project_analysis, p.id as project_id, ps.id as ps_id, r.name as run_name, d.type as device_type, r.id as run_id, ps.normal_id as normal_id, s.tumor as is_tumor, s.gender as gender, s.ffpe as is_ffpe, s.disease_group as disease_group, s.disease_status as disease_status, sys.type as sys_type, sys.target_file as sys_target, sys.name_manufacturer as sys_name, sys.name_short as sys_name_short, s.name_external as name_external ".
	                        "FROM project p, processed_sample ps, sample s, processing_system as sys, sequencing_run as r, device as d ".
				"WHERE ps.sequencing_run_id=r.id AND ps.project_id=p.id AND ps.sample_id=s.id AND s.name='$sample_name' AND ps.processing_system_id=sys.id AND ps.process_id='".(int)$process_id."' AND r.device_id = d.id");
	if (count($res)!=1)
	{
		if ($error_if_not_found)
		{
			trigger_error("Could not find information for processed sample with name '$ps_name' in NGSD!", E_USER_ERROR);
		}
		else
		{
			return null;
		}
	}
	$info = $res[0];
	
	if($info['normal_id']!="")
	{
		$info['normal_name'] = $db_conn->getValue("SELECT CONCAT(s.name,'_',LPAD(ps.process_id,2,'0')) FROM processed_sample as ps, sample as s WHERE ps.sample_id = s.id AND ps.id='".$info['normal_id']."'");
	}
	else
	{
		$info['normal_name'] = "";
	}
	
	//additional info
	$project_folder = get_path("project_folder");
	if (!ends_with($project_folder, "/")) $project_folder .= "/";
	$info['project_folder'] = $project_folder.$info['project_type']."/".$info['project_name']."/";
	$info['ps_name'] = $ps_name;
	$info['ps_folder'] = $info['project_folder']."Sample_{$ps_name}/";
	$info['ps_bam'] = $info['ps_folder']."{$ps_name}.bam";
	
	ksort($info);
	return $info;
}

//Converts NGSD sample meta data to a GSvar file header (using $override_map to allow replace NGSD info)
function gsvar_sample_header($ps_name, $override_map, $prefix = "##", $suffix = "\n")
{
	//get information from NGSD
	$parts = array();
	$parts['ID'] = $ps_name;
	
	if (db_is_enabled("NGSD"))
	{
		$db_conn = DB::getInstance("NGSD");
		$details = get_processed_sample_info($db_conn, $ps_name, false);
		$parts['Gender'] = $details['gender'];
		$parts['ExternalSampleName'] = strtr($details['name_external'], ",", ";");
		$parts['IsTumor'] = $details['is_tumor'] ? "yes" : "no";
		$parts['IsFFPE'] = $details['is_ffpe'] ? "yes" : "no";
		$parts['IsFFPE'] = $details['is_ffpe'] ? "yes" : "no";
		$parts['DiseaseGroup'] = $details['disease_group'];
		$parts['DiseaseStatus'] = $details['disease_status'];
	}
	
	//apply overwrite settings
	$parts = array_merge($parts, $override_map);
	
	//create and check output
	$output = array();
	$valid = array('ID','Gender','ExternalSampleName','SampleName','IsTumor','IsFFPE','DiseaseGroup','DiseaseStatus');
	foreach($parts as $key => $value)
	{
		if (!in_array($key, $valid))
		{
			trigger_error("Invalid GSvar sample header key '$key'! Valid are: ".implode(",", $valid), E_USER_ERROR); 
		}
		
		$output[] = "{$key}={$value}";
	}
	
	return "{$prefix}SAMPLE=<".implode(",", $output).">{$suffix}";
}

//Returns the index of the most similar column in a VCF header
function vcf_column_index($name, $header)
{
	//calculate distances (of prefix of same length)
	$dists = array();
	for ($i=9; $i<count($header); ++$i)
	{
		$min_len = min(strlen($name), strlen($header[$i]));
		$tmp_n = substr($name, 0, $min_len);
		$tmp_e = substr($header[$i], 0, $min_len);
		$dists[levenshtein($tmp_n, $tmp_e)][] = $i;
	}
	
	//determine minimum
	$min = min(array_keys($dists));
	$indices = $dists[$min];
	if (count($indices)>1)
	{
		$hits = array();
		foreach($indices as $index)
		{
			$hits[] = $header[$index];
		}
		trigger_error("Cannot determine sample column of '$name'. Samples '".implode("','", $hits)."' have the same edit distance ($min)!", E_USER_ERROR);
	}
	
	return $indices[0];
}

//checks whether gene names are up to date. Expects array with genes as input.
//Returns an array with approved symbols, obsolete gene names are replaced with up-to-data names
//If gene symbol is not found it is returned unaltered.
function approve_gene_names($input_genes)
{
	$genes_as_string = "";
	foreach($input_genes as $gene)
	{
		//set dummy if there are empty lines in input file
		if(trim($gene) == "")
		{
			$gene =  "NOT_AVAILABLE";
		}
		
		$genes_as_string = $genes_as_string.$gene."\n";
	}
	$non_approved_genes_file = tempnam(sys_get_temp_dir(),"temp_");
	file_put_contents($non_approved_genes_file,$genes_as_string);
	
	//write stdout to $approved_genes -> each checked gene is one array element
	$approved_genes = exec2(get_path("ngs-bits",true)."GenesToApproved -in $non_approved_genes_file")[0];
	
	$output = array();
	foreach($approved_genes as $gene)
	{
		//remove dummy before saving
		if($gene == "NOT_AVAILABLE") $gene = "";
		$output[] = (explode("\t",$gene))[0];
	}
	return $output;
}


//Returns information about an analysis job from the NGSD
function analysis_job_info(&$db_conn, $job_id, $error_if_not_found=true)
{
	$res = $db_conn->executeQuery("SELECT * FROM analysis_job WHERE id=:job_id", array("job_id"=>$job_id));
	if (count($res)==0)
	{
		if ($error_if_not_found)
		{
			trigger_error("Could not find information for analyis job with ID '$job_id' in NGSD!", E_USER_ERROR);
		}
		else
		{
			return null;
		}
	}
	$info = $res[0];
	
	//extract samples
	$info['samples'] = array();
	$res = $db_conn->executeQuery("SELECT processed_sample_id, info FROM analysis_job_sample WHERE analysis_job_id=:job_id ORDER BY id ASC", array("job_id"=>$job_id));
	foreach($res as $row)
	{
		$sample = $db_conn->getValue("SELECT CONCAT(s.name,'_',LPAD(ps.process_id,2,'0')) FROM processed_sample as ps, sample as s WHERE ps.sample_id = s.id AND ps.id='".$row['processed_sample_id']."'");
		$info['samples'][] = $sample."/".$row['info'];
	}
	
	//extract history
	$info['history'] = array();
	$res = $db_conn->executeQuery("SELECT status FROM analysis_job_history WHERE analysis_job_id=:job_id ORDER BY id ASC", array("job_id"=>$job_id));
	foreach($res as $row)
	{
		$info['history'][] = $row['status'];
	}
	
	return $info;
}

?>
