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

	if $cache_size is greater than 0, the requested region is extended by the given number and cached for further function calls

	@ingroup genomics
*/
function get_ref_seq($build, $chr, $start, $end, $cache_size=0, $use_local_data=true)
{
	// init vars for cache:
	static $cache_build = null;
	static $cache_chr = null;
	static $cache_start = null;
	static $cache_end = null;
	static $cache_sequence = null;

	//fix chromosome for GHCh37
	if ($chr=="chrM") $chr = "chrMT";

	if ($cache_size > 0)
	{
		if ($end < $start)
		{
			// report invalid chr position
			trigger_error("Error: Invalid chromosomal range {$chr}:{$start}-{$end} given for get_ref_seq()!", E_USER_ERROR);
		}
		// check if interval is in cache
		if ($build == $cache_build && $chr == $cache_chr && $start > $cache_start && $end < $cache_end)
		{
			// interval in cache -> return sequence from cache
			return substr($cache_sequence, $start - $cache_start, $end - $start + 1);
		}
		else
		{
			// interval not in cache -> create new cache
			$cache_build = $build;
			$cache_chr = $chr;
			$cache_start = max($start - $cache_size, 1);
			$cache_end = $end + $cache_size;

			// get sequence
			$output = array();
			list($output, $stderr, $ret) = execApptainer("samtools", "samtools faidx", genome_fasta($build, $use_local_data)." $chr:{$cache_start}-{$cache_end}", [genome_fasta($build, $use_local_data)], []);

			if ($ret!=0)
			{
				trigger_error("Error in get_ref_seq: ".implode("\n\t", $stderr), E_USER_ERROR);
			}
			elseif (!empty(array_filter($stderr)))
			{
				trigger_error("Warning in get_ref_seq: ".implode("\n\t", $stderr), E_USER_WARNING);
			}

			// check if chr range exceeds chr end:
			if (starts_with($output[0], "[faidx] Truncated sequence:"))
			{
				//skip warning
				$cache_sequence = trim(implode("", array_slice($output, 2)));
				// correct cached end position, if cache exceeds chr end
				$cache_end = $cache_start + strlen($cache_sequence);
			}
			else
			{
				$cache_sequence = trim(implode("", array_slice($output, 1)));
			}

			// return requested interval
			return substr($cache_sequence, $start - $cache_start, $end - $start + 1);
		}
	}
	else
	{
		// run without caching

		//get sequence
		$output = array();
		$genome = genome_fasta($build, $use_local_data);
		list($output, $stderr, $ret) = execApptainer("samtools", "samtools faidx", "{$genome} {$chr}:{$start}-{$end}", [$genome], []);

		if ($ret!=0)
		{
			trigger_error("Error in get_ref_seq: ".implode("\n\t", $stderr), E_USER_ERROR);
		}
		elseif (!empty(array_filter($stderr)))
		{
			trigger_error("Warning in get_ref_seq: ".implode("\n\t", $stderr), E_USER_WARNING);
		}

		return implode("", array_slice($output, 1));
	}
	
	
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
	@brief Returns the list of all chromosomes.
	@ingroup genomics
*/
function chr_list($include_chry = true)
{
	$chrs = range(1,22);
	$chrs[] = "X";
	if ($include_chry) $chrs[] = "Y";
	
	return preg_filter('/^/', 'chr', $chrs);
}

//Returns parsed INI file
function get_ini()
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
	
	return $parsed_ini;
}

/**
	@brief Function to get central organized paths to tools.
	
	Special handling for name 'megSAP-Tests' returns repository test/data/+analysis folder.
	
	@param name name of variable in ini-file.
	@return value of variable in ini-file (e.g. path).
		
	@ingroup helpers
*/
function get_path($name, $throw_on_error=true)
{
	$dir = repository_basedir();
	
	//get ini file
	$parsed_ini = get_ini();

	//check key name
	if (strpos($name, ".") !== false)
	{
		trigger_error("Key name '$name' contains '.' which is not supported!",E_USER_ERROR);
	}

	//get value (special handling for certain servers)
	$server = implode("", exec2("hostname")[0]);
	$name_server = $name.".".strtoupper(trim($server));
	if (isset($parsed_ini[$name_server]))
	{
		return $parsed_ini[$name_server];
	}
	
	if (!isset($parsed_ini[$name]) && $name == "megSAP-Tests")
	{
		return $dir."/test/data/+analysis/";
	}
	
	
	//get value
	if (!isset($parsed_ini[$name]) && $throw_on_error)
	{
		trigger_error("Could not find key '$name' in settings file!", E_USER_ERROR);
	}
	$value = isset($parsed_ini[$name]) ? $parsed_ini[$name] : "";

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
function get_db($db, $name, $default_value=null)
{
	$values = get_path($name);
	
	if (!isset($values[$db]))
	{
		if (is_null($default_value))
		{
			trigger_error("get_db could not find value '$name' for DB '$db'!", E_USER_ERROR);
		}
		else
		{
			return $default_value;
		}
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
		$db_conn = DB::getInstance("NGSD", false);
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
	$roi = trim($res[0]['target_file']);
	if ($roi!="") $roi = get_path("data_folder")."/enrichment/".$roi;
	$output[] = "target_file = \"".$roi."\"";
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
	static $cache = [];
	
	//not in cache > add it
	if (!isset($cache[$name]))
	{
		//split name
		$name = trim($name);
		$sep_idx = strpos($name, "_");
		if ($sep_idx==false)
		{
			$sname = $name;
			$id = "99";
		}
		else
		{
			$sname = substr($name, 0, $sep_idx);
			$id = substr($name, $sep_idx+1, 2);
			$id = ltrim($id, "0");
		}
		
		//query NGSD
		$res = $db_conn->executeQuery("SELECT ps.id FROM processed_sample ps, sample s WHERE ps.sample_id=s.id AND s.name=:name AND ps.process_id=:id", array('name' => $sname, "id"=>$id));
		
		if (count($res)<1)
		{
			if ($error_if_not_found) trigger_error("Could not find processed sample with name '$name' in NGSD!", E_USER_ERROR);
			return -1;
		}
		
		$cache[$name] = $res[0]['id'];
	}
	
	return $cache[$name];
}


///Returns variant ID or -1 if variant is not present in NGSD
function get_variant_id(&$db_conn, $chr, $start, $end, $ref, $obs)
{
	return $db_conn->getValue("SELECT id FROM variant WHERE chr='{$chr}' AND start='{$start}' AND end='{$end}' AND ref='{$ref}' AND obs='{$obs}'", -1);
}

///returns array with processed sample names in the form DX******_** from run name
function get_processed_samples_from_run(&$db, $run_id)
{
	$query = "SELECT ps.process_id, s.name FROM sequencing_run as run, processed_sample as ps, sample as s WHERE run.name = '{$run_id}' AND ps.sequencing_run_id = run.id AND s.id = ps.sample_id";
	
	$res = $db->executeQuery($query);
	
	$sample_names = array();
	
	foreach($res as $data)
	{
		$sample_names[] = $data['name'] . "_" . str_pad($data['process_id'],2,'0',STR_PAD_LEFT);
	}
	
	return $sample_names;
}

///Updates normal sample entry for given tumor sample.
function updateNormalSample(&$db_conn, $ps_tumor, $ps_normal, $overwrite = false)
{
	//tumor
	$ps_tid = get_processed_sample_id($db_conn, $ps_tumor);
	$n = trim($db_conn->getValue("SELECT normal_id FROM processed_sample WHERE id='$ps_tid'"));

	//normal
	$ps_nid = get_processed_sample_id($db_conn, $ps_normal);
	
	//check not set
	if($n!="" && $ps_nid!=$n) trigger_error("Normal sample for tumor {$ps_tumor} already set in NGSD (ID={$n}) - cannot set {$ps_normal} as normal sample!", E_USER_WARNING);
	
	//set
	if($n=="" || $overwrite)
	{
		$db_conn->executeStmt("UPDATE processed_sample SET normal_id='$ps_nid' WHERE id='$ps_tid'");
		return true;
	}

	return false;
}

//Convert variant in GSvar format to VCF format.
function gsvar_to_vcf($build, $chr, $start, $ref, $obs)
{
	$ref = trim($ref,"-");
	$obs = trim($obs,"-");

	//get extra base before indels from reference genome
	$extra1 = "";
	if(strlen($ref)==0 && strlen($obs)>=1) // simple insertion - add prefix
	{
		$extra1 = get_ref_seq($build, $chr,$start,$start);
	}
	else if(strlen($ref)>=1 && strlen($obs)==0) // simple deletion - correct start and add prefix
	{
		$start -= 1;
		$extra1 = get_ref_seq($build, $chr,$start,$start);
	}
	else if (strlen($ref)==strlen($obs)) //SNP or MNP - nothing to do
	{
		
	}
	else if(strlen($ref)>1 || strlen($obs)>1) // complex indel - correct start and add prefix
	{
		$start -= 1;
		$extra1 = get_ref_seq($build, $chr,$start, $start);
	}

	// combine all information
	$ref = strtoupper($extra1.$ref);
	$obs = strtoupper($extra1.$obs);
	
	return array($chr, $start, $ref, $obs);
}

function is_valid_ref_sample_for_cnv_analysis($file, $tumor_only = false, $include_test_projects = false, $include_ffpe = false)
{
	//check that it is not a NIST reference sample (cell-lines and sequenced very often)
	if (contains($file, "NA12878")) return false;
	if (contains($file, "NA12891")) return false;
	if (contains($file, "NA12892")) return false;
	if (contains($file, "NA24385")) return false;

	//no NGSD => error
	if (!db_is_enabled("NGSD"))
	{
		trigger_error("is_valid_ref_sample_for_cnv_analysis needs NGSD access!", E_USER_ERROR);
	}

	//check sample is in NGSD
	$db_conn = DB::getInstance("NGSD", false);
	$ps_id = get_processed_sample_id($db_conn, $file, false);
	if ($ps_id<0) return false;
	
	//check that sample is not tumor and not not FFPE
	$res = $db_conn->executeQuery("SELECT s.tumor, s.ffpe, ps.quality q1, r.quality q2, p.type FROM sequencing_run r, sample s, processed_sample ps, project p WHERE s.id=ps.sample_id AND ps.id='$ps_id' AND ps.sequencing_run_id=r.id AND ps.project_id = p.id");
	if (count($res)==0) return false;
	if ($tumor_only)
	{
		if ($res[0]['tumor']!="1") return false;
	}
	else
	{
		if ($res[0]['tumor']=="1") return false;
		if ($res[0]['ffpe']=="1" && !$include_ffpe) return false;
	}
	
	//check that run and processed sample do not have bad quality
	if ($res[0]['q1']=="bad") return false;
	if ($res[0]['q2']=="bad") return false;
	
	//include test projects only if requested
	if ($res[0]['type']=="test" && !$include_test_projects) return false;
	
	return true;
}

function is_valid_ref_tumor_sample_for_cnv_analysis($file, $discard_ffpe = false, $include_test_projects = false)
{
	//check that sample is not NIST reference sample (it is a cell-line)
	if (contains($file, "NA12878")) return false;

	//no NGSD => error
	if (!db_is_enabled("NGSD"))
	{
		trigger_error("is_valid_ref_tumor_sample_for_cnv_analysis needs NGSD access!", E_USER_ERROR);
	}
	
	//check sample is in NGSD
	$db_conn = DB::getInstance("NGSD", false);
	$ps_id = get_processed_sample_id($db_conn, $file, false);
	
	//check that sample is not FFPE
	$res = $db_conn->executeQuery("SELECT s.tumor, s.ffpe, ps.quality q1, r.quality q2, p.type FROM sequencing_run r, sample s, processed_sample ps, project p WHERE s.id=ps.sample_id AND ps.id='$ps_id' AND ps.sequencing_run_id=r.id AND ps.project_id = p.id");
	if ($ps_id<0) return false;
	if ($discard_ffpe && $res[0]['ffpe']=="1") return false;
	
	//check that run and processed sample do not have bad quality
	if ($res[0]['q1']=="bad") return false;
	if ($res[0]['q2']=="bad") return false;
	
	//check that project type is research/diagnostics
	if ($res[0]['type'] =="test" && !$include_test_projects) return false;
	
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
	if (!file_exists($filename)) return  [];

	$comments = array();
	$header = "";
	$vars = array();
	
	//load and normalize data
	foreach(file($filename) as $line)
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
			if (count($parts)<10) trigger_error("VCF file $filename has line with less than 10 columns: $line", E_USER_ERROR);
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
	if (!preg_match("/^[acgtACGT]$/", $obs))
	{
		trigger_error("Invalid observed allele ({$format_col} {$sample_col} {$obs}).", E_USER_ERROR);
	}

	$obs = strtoupper($obs);
	$f = explode(":", $format_col);
	$s = explode(":", $sample_col);
	$data = array_combine($f, $s);
	$required_fields = [ "DP", "TU", "AU", "CU", "GU" ];
	foreach ($required_fields as $req)
	{
		if (!array_key_exists($req, $data))
		{
			trigger_error("Invalid strelka format: field {$req} not available.", E_USER_ERROR);
		}
	}

	//extract tier 1 counts
	$observed_counts = [];
	foreach ([ "T", "A", "C", "G" ] as $base)
	{
		$observed_count = $data["{$base}U"];
		list($tier1, $tier2) = explode(",", $observed_count, 2);
		$observed_counts[$base] = $tier1;
	}

	$sum = array_sum($observed_counts);
	$af = 0;
	if ($sum != 0)
	{
		$af = $observed_counts[$obs] / $sum;
	}

	return [ $data["DP"], $af ];
}

/**
 * Identifies bases above frequency cutoff from strelka2 VCF records.
 *
 * @param format_col	VCF format column
 * @param sample_col	VCF sample column
 * @param ref			reference base
 * @param obs			observed base
 * @param min_taf		minimum frequency required to call alternativ observation
 *
 * @return array Alternative observation base and frequency.
 */
function vcf_strelka_snv_postcall($format_col, $sample_col, $ref, $obs, $min_taf)
{
	$postcalls = [];
	foreach (array_diff([ "T", "A", "C", "G" ], [ $ref, $obs ]) as $base)
	{
		list($depth, $freq) = vcf_strelka_snv($format_col, $sample_col, $base);

		if ($freq >= $min_taf)
		{
			$postcalls[] = [ $base, $freq ];
		}
	}
	return $postcalls;
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

	//tir and tar contain strong supporting reads, tor (not considered here) contains weak supportin reads like breakpoints
	//only strong supporting reads are used for calculation of allele fraction
	$f = ($tir+$tar)==0 ? 0.0 : $tir/($tir+$tar);
	
	return array($d,$f);
}

//Extracts depth and variant allele frequency for umivar2 output line.
function vcf_umivar2($filter_col, $info_col, $format_col, $sample_col)
{
	$filter_data = explode(";", $filter_col);
	$info_data = explode(";", $info_col);
	$format_data = explode(":", $format_col);
	$sample_data = explode(":", $sample_col);
	if(count($format_data) != count($sample_data))
	{
		trigger_error("FORMAT column and sample column have different count of entries.", E_USER_ERROR);
	}
	
	//parse filter column:
	$homopolymer = "false";
	foreach ($filter_data as $filter) 
	{
		if(trim($filter) == "Homopolymer") 
		{
			$homopolymer = "true";
		}
	}

	//parse info column
	$seq_context = NULL;
	foreach ($info_data as $info) 
	{
		if(starts_with($info, "Vicinity="))
		{
			$seq_context = substr($info, 9);
		}
	}
	if(is_null($seq_context))
	{
		trigger_error("Invalid umiVar2 format. Vicinity field is missing in INFO column!", E_USER_ERROR);
	}
	//parse format cols
	$i_alt_count = NULL; //index alt count
	$i_dp = NULL; //index depth
	$i_af = NULL; //index allele frequency
	$i_strand = NULL; //index strand
	$i_p_value = NULL; //index P-value
	$i_m_af = NULL; //index multi UMI allele frequency
	$i_m_ref_count = NULL; //index multi UMI ref count
	$i_m_alt_count = NULL; //index multi UMI alt count

	for($i=0;$i<count($format_data);++$i)
	{
		if($format_data[$i] == "AC") $i_alt_count = $i;
		if($format_data[$i] == "DP") $i_dp = $i;
		if($format_data[$i] == "AF") $i_af = $i;
		if($format_data[$i] == "Strand") $i_strand = $i;
		if($format_data[$i] == "Pval") $i_p_value = $i;

		if($format_data[$i] == "M_AF") $i_m_af = $i;
		if($format_data[$i] == "M_REF") $i_m_ref_count = $i;
		if($format_data[$i] == "M_AC") $i_m_alt_count = $i;
	}
	
	if(is_null($i_alt_count) || is_null($i_dp) || is_null($i_af) || is_null($i_strand) || is_null($i_p_value))
	{
		trigger_error("Invalid umiVar2 format. Missing one of the fields 'AC', 'DP', 'AF', 'Strand' or 'Pval'", E_USER_ERROR);
	}

	if(is_null($i_m_af) || is_null($i_m_ref_count) || is_null($i_m_alt_count))
	{
		trigger_error("Invalid umiVar2 format. Missing one of the fields 'M_AF', 'M_REF' or 'M_AC'", E_USER_ERROR);
	}

	return array($sample_data[$i_dp], number_format($sample_data[$i_af], 5), $sample_data[$i_p_value], $sample_data[$i_alt_count], $sample_data[$i_strand], $seq_context, 
				$homopolymer, number_format($sample_data[$i_m_af], 5), $sample_data[$i_m_ref_count], $sample_data[$i_m_alt_count]);
}

//Calculates depth and variant allele frequency for varscan2 output line.
function vcf_varscan2($format_col, $sample_col)
{
	//combine FORMAT and SAMPLE data
	$format_data = explode(":", $format_col);
	$sample_data = explode(":", $sample_col);
	if(count($format_data) != count($sample_data)) trigger_error("FORMAT column and sample column have different count of entries.", E_USER_ERROR);
	$data = array_combine($format_data, $sample_data);
	
	//check it's from VarScan2
	if (!isset($data['DP']) || !isset($data['RD']) || !isset($data['AD']) || !isset($data['FREQ']))
	{
		trigger_error("Invalid Varscan2 format. Missing one of the fields 'DP', 'RD', 'AD' or 'FREQ'", E_USER_ERROR);
	}
	
	//convert AF
	$af = trim(strtr($data['FREQ'],'%', ' '));
	$af  = number_format($af/100.0, 4);
	
	return array($data['DP'], $af);
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

function vcf_deepvariant($format_col, $sample_col)
{
	$g = explode(":",$format_col);
	$index_DP = NULL;
	$index_VAF = NULL;
	$index_GT = NULL;
	for($i=0;$i<count($g);++$i)
	{
		if($g[$i]=="DP")	$index_DP = $i;
		if($g[$i]=="VAF")	$index_VAF = $i;
		if($g[$i]=="GT")	$index_GT = $i;
	}
	if(is_null($index_DP) || is_null($index_VAF) ||is_null($index_GT))	trigger_error("Invalid DeepVariant format; either field DP, GT or VAF not available.",E_USER_ERROR);	
	
	$s = explode(":",$sample_col);

	if(!is_numeric($s[$index_DP]))	trigger_error("Could not identify numeric depth (".$format_col." ".$sample_col.").",E_USER_ERROR);
	if(!is_numeric($s[$index_VAF]))	trigger_error("Invalid variant allele frequency (".$format_col." ".$sample_col.").",E_USER_ERROR);	// currently no multiallelic variants supported

	$dp = $s[$index_DP];
	$fq = $s[$index_VAF];
	
	return array($dp,$fq);
}

function vcf_deepsomatic($info_col, $n_name)
{
	$g = explode(";", $info_col);
	$af = NULL;
	$dp = NULL;
	for($i=0;$i<count($g);++$i)
	{
		$parts = explode("=", $g[$i]);
		if ($parts[0] == "{$n_name}_AF") $af = $parts[1];
		if ($parts[0] == "{$n_name}_DP") $dp = $parts[1];
	}

	if(is_null($af) || is_null($dp)) trigger_error("Normal AF or DP not annotated in the info filed of Deepsomatic tumor-normal output. (info field: {$info_col})", E_USER_ERROR);

	if(!is_numeric($af)) trigger_error("Invalid allele frequenzy (AF: '".$af."' in info field: '".$info_col."').", E_USER_ERROR);	// currently no multiallelic variants supported
	if(!is_numeric($dp)) trigger_error("Could not identify numeric depth (DP: '".$dp."' in info field: '".$info_col."').", E_USER_ERROR);

	return array($dp,$af);
}

function vcf_dragen_var($format_col, $sample_col)
{
	$g = explode(":",$format_col);
	$index_DP = NULL;
	$index_AF = NULL;
	for($i=0;$i<count($g);++$i)
	{
		if($g[$i]=="DP")	$index_DP = $i;
		if($g[$i]=="AF")	$index_AF = $i;
	}

	if(is_null($index_DP) || is_null($index_AF)) trigger_error("Invalid dragon vcf format; either field DP or AF not available.", E_USER_ERROR);	
	
	$s = explode(":",$sample_col);	
	
	if(!is_numeric($s[$index_AF]))	trigger_error("Invalid allele frequenzy (".$format_col." ".$sample_col.").", E_USER_ERROR);	// currently no multiallelic variants supported
	if(!is_numeric($s[$index_DP]))	trigger_error("Could not identify numeric depth (".$format_col." ".$sample_col.").", E_USER_ERROR);

	$dp = $s[$index_DP];
	$af = $s[$index_AF];
	
	return array($dp,$af);
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
	
	if ($gt=="0" || $gt=="0|0")
	{
		$geno = "wt";
	}
	else if ($gt=="0|1" || $gt=="1|0")
	{
		$geno = "het";
	}
	else if ($gt=="1|1" || $gt=="1") //'1' is used by Dragen for hemizygous variants on chrX/chrY in males
	{
		$geno = "hom";
	}
	else
	{
		trigger_error("Invalid VCF genotype '$gt'!", E_USER_ERROR);
	}
	
	return $upper_case ? strtoupper($geno) : $geno;
}

function processed_sample_name(&$db_conn, $ps_id)
{
	return $db_conn->getValue("SELECT CONCAT(s.name,'_',LPAD(ps.process_id,2,'0')) FROM processed_sample as ps, sample as s WHERE ps.sample_id = s.id AND ps.id={$ps_id}");
}

//Returns information from the NGSD about a processed sample (as key-value pairs), or null/error if the sample is not found.
function get_processed_sample_info(&$db_conn, $ps_name, $error_if_not_found=true, $no_paths=false)
{
	//get info from NGSD
	$ps_name = trim($ps_name);
	list($sample_name, $process_id) = explode("_", $ps_name."_");
	$res = $db_conn->executeQuery("SELECT p.name as project_name, p.type as project_type, p.analysis as project_analysis, p.preserve_fastqs as preserve_fastqs, p.id as project_id, ps.id as ps_id, r.name as run_name, d.type as device_type, r.id as run_id, ps.normal_id as normal_id, s.sample_type as sample_type, s.tumor as is_tumor, s.gender as gender, s.ffpe as is_ffpe, s.disease_group as disease_group, s.disease_status as disease_status, s.tissue as tissue, s.received as s_received, s.comment as s_comments, s.patient_identifier, s.year_of_birth, s.order_date, sys.id as sys_id, sys.type as sys_type, sys.target_file as sys_target, sys.name_manufacturer as sys_name, sys.name_short as sys_name_short, sys.adapter1_p5 as sys_adapter1, sys.adapter2_p7 as sys_adapter2, sys.umi_type as sys_umi_type, g.build as sys_build, s.name_external as name_external, ps.comment as ps_comments, ps.lane as ps_lanes, r.recipe as run_recipe, r.fcid as run_fcid, ps.mid1_i7 as ps_mid1, ps.mid2_i5 as ps_mid2, sp.id as species_id, sp.name as species, s.id as s_id, ps.quality as ps_quality, ps.processing_input, d.name as device_name, psa.population as ancestry, ps.urgent ".
									"FROM project p, sample s, processing_system as sys, processed_sample ps LEFT JOIN sequencing_run as r ON ps.sequencing_run_id=r.id LEFT JOIN device as d ON r.device_id=d.id LEFT JOIN processed_sample_ancestry as psa ON ps.id=psa.processed_sample_id, species sp, genome g ".
									"WHERE ps.project_id=p.id AND ps.sample_id=s.id AND s.name='$sample_name' AND ps.processing_system_id=sys.id AND sys.genome_id=g.id AND s.species_id=sp.id AND ps.process_id='".(int)$process_id."'");
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
	$info['ps_name'] = $ps_name;
	
	//prefix target region
	$roi = trim($info['sys_target']);
	if ($roi!="") $roi = get_path("data_folder")."/enrichment/".$roi;
	$info['sys_target'] = $roi;
	
	
	//normal sample name (tumor reference sample)
	if($info['normal_id']!="")
	{
		$info['normal_name'] = processed_sample_name($db_conn, $info['normal_id']);
	}
	else
	{
		$info['normal_name'] = "";
	}
	
	//replace fields to make them human_readable
	if($info['ps_mid1']!="")
	{
		$info['ps_mid1'] = $db_conn->getValue("SELECT name FROM mid WHERE id='".$info['ps_mid1']."'");
	}
	if($info['ps_mid2']!="")
	{
		$info['ps_mid2'] = $db_conn->getValue("SELECT name FROM mid WHERE id='".$info['ps_mid2']."'");
	}
	
	//split fields if required
	$info['ps_comments'] = array_map("trim", explode("\n", $info['ps_comments']));
	$info['ps_lanes'] = array_map("trim", explode(",", $info['ps_lanes']));
	$info['s_comments'] = array_map("trim", explode("\n", $info['s_comments']));
	
	//add paths if requested (quite slow)
	if (!$no_paths)
	{
		if ($info['project_type'] == "megSAP-Tests") //special handling for tests
		{
			$info['project_folder'] = get_path($info['project_type'])."/".$info['project_name']."/";
			$info['ps_folder'] = $info['project_folder']."Sample_{$ps_name}/";
			$info['ps_bam'] = $info['ps_folder']."{$ps_name}.cram";
			if (!file_exists($info['ps_bam'])) $info['ps_bam'] = $info['ps_folder']."{$ps_name}.bam"; //fallback to BAM
		}
		else
		{
			//get BAM path via ngs-bits SamplePath tool (handles override paths for project/processed sample and BAM/CRAM)
			$args = [];
			$args[] = "-ps {$ps_name}";
			$args[] = "-type SAMPLE_FOLDER";
			if ($db_conn->name()=="NGSD_TEST") $args[] = "-test";
			list ($stdout, $stderr) = execApptainer("ngs-bits", "SamplePath", implode(" ", $args));
			$ps_folder = trim(implode("", $stdout));
			if(!ends_with($ps_folder, '/')) $ps_folder .= '/';
			$info['ps_folder'] = $ps_folder;
			$ps_bam_or_cram = $ps_folder.$ps_name.".cram";
			if (!file_exists($ps_bam_or_cram)) $ps_bam_or_cram = $ps_folder.$ps_name.".bam"; //fallback to BAM
			$info['ps_bam'] = $ps_bam_or_cram;
			$project_folder = dirname($ps_folder);
			if(!ends_with($project_folder, '/')) $project_folder .= '/';
			$info['project_folder'] = $project_folder;
		}
	}
	
	ksort($info);
	return $info;
}

//Returns the value of the given QC term. Throws an error if the QC term is unknown or the sample is not found. Returns a empty string if the QC entry is not in NGSD.
function get_processed_sample_qc(&$db, $ps, $qc_name_or_accession)
{
	//determine QC term ID
	$term_id = $db->getValue("SELECT id FROM qc_terms WHERE name LIKE '{$qc_name_or_accession}' OR qcml_id LIKE '{$qc_name_or_accession}'");
	
	$ps_id = get_processed_sample_id($db, $ps);
	
	return $db->getValue("SELECT value FROM processed_sample_qc WHERE processed_sample_id='{$ps_id}' AND qc_terms_id={$term_id}", "");
}

//Returns the value of the given QC term from the file. If default is not set, throws an error if the file does not exist, or does not contain the term.
function get_qcml_value($qcml_file, $accession, $default=null)
{
	if (!file_exists($qcml_file)) 
	{
		if (is_null($default)) trigger_error("qcML file does not exists: {$qcml_file}", E_USER_ERROR);
		return $default;
	}
	
	
	$value = "";
	$xml = simplexml_load_file($qcml_file);
	if ($xml===false) trigger_error("Could not load XML file: $qcml_file", E_USER_ERROR);
	foreach ($xml->runQuality->qualityParameter as $qp)
	{
		if ($qp["accession"]==$accession)
		{
			return $qp["value"];
		}
	}
	
	if (is_null($default)) trigger_error("QC term with ID '{$accession}' not found in qcML file: {$qcml_file}", E_USER_ERROR);
	return $default;
}

//Creates a samples header using meta data from NGSD (using $override_map to allow replace NGSD info)
function gsvar_sample_header($ps, $override_map, $prefix = "##", $suffix = "\n")
{
	//init
	$parts = [];
	$parts['ID'] = $ps;
	$parts['Gender'] = "";
	$parts['ExternalSampleName'] = "";
	$parts['IsTumor'] = "";
	$parts['IsFFPE'] = "";
	$parts['DiseaseGroup'] = "";
	$parts['DiseaseStatus'] = "";
	
	//get data from NGSD
	if (db_is_enabled("NGSD"))
	{
		$db = DB::getInstance("NGSD", false);
		$info = get_processed_sample_info($db, $ps, false);
		if ($info!=NULL)
		{
			$parts['Gender'] = $info['gender'];
			if ($parts['Gender'] == 'n/a') trigger_error("Gender of sample {$ps} not set in NGSD!", E_USER_WARNING);
			$parts['ExternalSampleName'] = strtr($info['name_external'], ",", ";");
			$parts['IsTumor'] = $info['is_tumor'] ? "yes" : "no";
			$parts['IsFFPE'] = $info['is_ffpe'] ? "yes" : "no";
			$parts['DiseaseGroup'] = strtr($info['disease_group'], ",", ";");
			$parts['DiseaseStatus'] = $info['disease_status'];
		}
	}
	
	//apply overwrite settings
	foreach($override_map as $key => $value)
	{
		if (!isset($parts[$key])) trigger_error("Invalid GSvar sample header key '$key'! Valid are: ".implode(",", array_keys($parts)), E_USER_ERROR); 
		
		$parts[$key] = $value;
	}
	
	//create output
	$output = [];
	foreach($parts as $key => $value)
	{
		//skip empty entries
		if ($value=="") continue;
		
		$output[] = "{$key}=".trim($value);
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
	//no NGSD > return input
	if (!db_is_enabled("NGSD")) return $input_genes;
	
	//write file with one line per gene
	$genes_as_string = "";
	foreach($input_genes as $gene)
	{
		$gene = trim($gene);
		
		//set dummy if there are empty lines in input file
		if($gene=="") $gene = "NOT_AVAILABLE";
		
		$genes_as_string .= $gene."\n";
	}
	$non_approved_genes_file = temp_file(".txt", "approve_gene_names");
	file_put_contents($non_approved_genes_file, $genes_as_string);
	
	//convert to approved gene names
	list($approved_genes) = execApptainer("ngs-bits", "GenesToApproved", "-in $non_approved_genes_file");
	
	//get new names (first element of tab-separated line)
	$output = [];
	foreach($approved_genes as $line)
	{
		//remove dummy before saving
		if(starts_with($line, "NOT_AVAILABLE")) $line = "";
		
		$output[] = explode("\t", $line)[0];
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
		$info['samples'][] = processed_sample_name($db_conn, $row['processed_sample_id'])."/".$row['info'];
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

//Returns if a BED, VCF or VCF.GZ file contains mito lines
function contains_mito($filename)
{
	$h = gzopen($filename, 'r');
	while(!gzeof($h))
	{
		$line = gzgets($h);
		if (starts_with($line, "chrMT\t"))
		{
			return true;
		}
	}
	gzclose($h);
	
	return false;
}

//Returns if special variant calling for mitochondrial DNA should be performed
function enable_special_mito_vc($sys)
{
	//only for GRCh38
	if ($sys['build']!="GRCh38") return false;
	
	//only exome/genome
	if ($sys['type']!="WES" && $sys['type']!="WGS") return false;
	
	//no ROI > chrMT is called anyway > no special calling
	if ($sys['target_file']=="") return false;
	
	//check if chrMT is already contained in the ROI
	if (contains_mito($sys['target_file'])) return false;
	
	return true;
}

//Returns the FASTA file of the genome
function genome_fasta($build, $use_local_data=true, $use_local_ramdrive=true)
{
	//allow overriding genome fasta
	$genome_fasta_override = trim(get_path("genome_fasta_override", false));
	if  ($genome_fasta_override!="")
	{
		return $genome_fasta_override;
	}
	
	//use genome FASTA of the local maching
	if ($use_local_data)
	{
		//prefer ram drive
		if($use_local_ramdrive)
		{
			$ramdrive_fasta = "/mnt/genome_ramdrive/".$build.".fa";
			if (file_exists($ramdrive_fasta)) return $ramdrive_fasta;
		}
		
		//use local copy in tmp
		$local_fasta = get_path("local_data")."/".$build.".fa";
		if (file_exists($local_fasta)) return $local_fasta;
		else trigger_error("Use of local genome file requested, but '$local_fasta' does not exist.", E_USER_NOTICE);
	}
	
	//use the genome FASTA from the megSAP installation
	return get_path("data_folder")."/genomes/".$build.".fa";
}

//Create Bed File that contains off target regions of a target region
function create_off_target_bed_file($out,$target_file,$ref_genome_fasta)
{
	//generate bed file that contains edges of whole reference genome
	$handle_in = fopen2("{$ref_genome_fasta}.fai","r");
	$ref_bed = temp_file(".bed");
	$handle_out = fopen2($ref_bed,"w");
	while(!feof($handle_in))
	{
		$line = trim(fgets($handle_in));
		if(empty($line)) continue;
		list($chr,$chr_length) = explode("\t",$line);
		if(chr_check($chr,22,false) === false) continue;
		$chr_length--; //0-based coordinates
		fputs($handle_out,"{$chr}\t0\t{$chr_length}\n");
	}
	fclose($handle_in);
	fclose($handle_out);
	
	//Create off target bed file
	$tmp_bed = temp_file(".bed");

	$command_bed_extend = execApptainer("ngs-bits", "BedExtend", "-in ".$target_file." -n 1000 -fai {$ref_genome_fasta}.fai", [$target_file, $ref_genome_fasta], [], true);
	$command_bed_merge = execApptainer("ngs-bits", "BedMerge", "-out {$tmp_bed}", [], [], true);

	exec2("{$command_bed_extend} | {$command_bed_merge}");

	$command_bed_subtract = execApptainer("ngs-bits", "BedSubtract", "-in ".$ref_bed." -in2 {$tmp_bed}", [], [], true);
	$command_bed_chunk = execApptainer("ngs-bits", "BedChunk", "-n 100000", [], [], true);
	$command_bed_shrink = execApptainer("ngs-bits", "BedShrink", "-n 25000", [], [], true);
	$command_bed_extend = execApptainer("ngs-bits", "BedExtend", "-n 25000 -fai {$ref_genome_fasta}.fai", [$ref_genome_fasta], [], true);
	$command_bed_annotate_gc = execApptainer("ngs-bits", "BedAnnotateGC", "-ref {$ref_genome_fasta}", [$ref_genome_fasta], [], true);
	$command_bed_annotate_genes = execApptainer("ngs-bits", "BedAnnotateGenes ", "-out {$out}", [], [dirname($out)], true);

	exec2 ("{$command_bed_subtract} | {$command_bed_chunk} | {$command_bed_shrink} |{$command_bed_extend} | {$command_bed_annotate_gc} | {$command_bed_annotate_genes}");
	execApptainer("ngs-bits", "BedSort", "-uniq -in $out -out $out", [], [dirname($out)]);
}

//returns the allele counts for a sample at a certain position as an associative array, reference skips and start/ends of read segments are ignored
function allele_count($bam, $chr, $pos)
{
	if (ends_with($bam, ".cram"))
	{
		trigger_error("This function uses a samtools call without a reference. If you're using CRAM files check first if samtools downloads the reference chromosomes implicitly.", E_USER_ERROR);
		//output of "samtools mpileup" also changes with reference if you add a reference you need to make the function compatible with the new output.
	}
	
	//get pileup
	list($output) = execApptainer("samtools", "samtools mpileup", "-aa -r $chr:$pos-$pos $bam", [$bam]);
	list($chr2, $pos2, $ref2,, $bases) = explode("\t", $output[0]);
	
	//count bases
	$bases = strtoupper($bases);
	$counts = array("A"=>0, "C"=>0, "G"=>0, "T"=>0, "*"=>0); //4 bases and "*" denoting an insertion or deletion
	
	for($i=0; $i<strlen($bases); ++$i)
	{
		$char = $bases[$i];
		
		if($char == "^") break; //skip qual bases/segment start at the end of pileup string
		if($char == "$") continue;
		
		if (isset($counts[$char]))
		{
			++$counts[$char];
		}
		elseif($char == "-" || $char == "+") //start of an deletion or insertion
		{
			$indel_size_as_string = "";
			++$i;
			while(is_numeric($bases[$i])) //Determine size of the indel, is decimal with multiple places in orgiinal pileup format, e.g. "+99AGTC...."
			{
				$indel_size_as_string .= $bases[$i];
				++$i;
			}

			$indel_size = intval($indel_size_as_string);
			
			$counts["*"] += 1;
			$i += $indel_size - 1; //skip bases that belong to the descroption of the indel, -1 is neccessary due to for loop
		}
	}
	arsort($counts);
	
	return $counts;
}

/**
	@brief Returns an array with (1) report configuration id of processed sample or -1 (2) if small variants report config exists and (3) CNV report config exists
 */
function report_config(&$db_conn, $name, $error_if_not_found=false)
{
	$ps_id = get_processed_sample_id($db_conn, $name, $error_if_not_found);
		
	$rc_id = $db_conn->getValue("SELECT id FROM report_configuration WHERE processed_sample_id=".$ps_id, -1);
	
	$var_ids = $db_conn->getValues("SELECT id FROM report_configuration_variant WHERE report_configuration_id=".$rc_id);
	$cnv_ids = $db_conn->getValues("SELECT id FROM report_configuration_cnv WHERE report_configuration_id=".$rc_id);
	$sv_ids  = $db_conn->getValues("SELECT id FROM report_configuration_sv WHERE report_configuration_id=".$rc_id);
	$re_ids  = $db_conn->getValues("SELECT id FROM report_configuration_re WHERE report_configuration_id=".$rc_id);
	
	
	return array($rc_id, count($var_ids)>0, count($cnv_ids)>0, count($sv_ids)>0, count($re_ids)>0);
}

function somatic_report_config(&$db_conn, $t_ps, $n_ps, $error_if_not_found=false)
{
	$t_ps_id = get_processed_sample_id($db_conn, $t_ps, $error_if_not_found);
	$n_ps_id = get_processed_sample_id($db_conn, $n_ps, $error_if_not_found);
	
	$config_id = $db_conn->getValue("SELECT id FROM somatic_report_configuration WHERE ps_tumor_id='{$t_ps_id}' AND ps_normal_id='{$n_ps_id}'", -1);
	
	if($error_if_not_found && $config_id == -1)
	{
		trigger_error("Could not find somatic report id for $t_ps_id {$n_ps_id}.", E_USER_ERROR);
	}
	
	$var_ids = $db_conn->getValues("SELECT id FROM somatic_report_configuration_variant WHERE somatic_report_configuration_id=".$config_id);
	$cnv_ids = $db_conn->getValues("SELECT id FROM somatic_report_configuration_cnv WHERE somatic_report_configuration_id=".$config_id);
	$sv_ids  = $db_conn->getValues("SELECT id FROM somatic_report_configuration_sv WHERE somatic_report_configuration_id=".$config_id);
	
	return array($config_id, count($var_ids)>0, count($cnv_ids) > 0, count($sv_ids) > 0);	
}

//Returns cytobands of given genomic range as array
function cytobands($chr, $start, $end)
{
	$handle = fopen2(repository_basedir()."/data/misc/cytoBand.txt", "r");
	
	$out = array();
	
	while(!feof($handle))
	{
		$line = trim(fgets($handle));
		if(empty($line)) continue;
		list($cchr, $cstart, $cend, $cname) = explode( "\t", $line );
		if($chr != $cchr) continue;
		
		if(range_overlap($start, $end, $cstart, $cend))
		{
			$out[] =  $cname;
		}
	}
	fclose($handle);

	asort($out);
	
	return $out;
}

//checks if any contig line is given, if not adds all contig lines from reference genome (using fasta index file)
function add_missing_contigs_to_vcf($build, $vcf)
{
	$file = fopen2($vcf, 'c+');
	$new_file_lines = array();
	if($file)
	{
		$new_file_lines = explode("\n", fread($file, filesize($vcf)));
		fseek($file, 0);
	}
	else
	{
		trigger_error("Could not open file ".$vcf.": no new contig lines written.", E_USER_WARNING);
	}

	$contains_contig = false;
	$line_below_reference_info = 0;
	$count = 0;

	$new_contigs = array();

	while(!feof($file))
	{
		$count += 1;
		$line = trim(fgets($file));

		if(starts_with($line, "##reference"))
		{
			$line_below_reference_info = $count;
		}
		else if (starts_with($line, "##"))
		{
			if(starts_with($line, "##contig"))
			{
				$contains_contig = true;
				break;
			}
		}
		else
		{
			//write contigs in second line if no reference genome line is given
			if($line_below_reference_info == 0)
			{
				$line_below_reference_info = 1;
			}
			break;
		}
	}
	if(!$contains_contig)
	{
		$fai_file_path = genome_fasta($build).".fai";
		if (!file_exists($fai_file_path))
		{
			trigger_error("Fasta index file \"{$fai_file_path}\" is missing!", E_USER_ERROR);
		}
		$fai_file_content = file($fai_file_path, FILE_IGNORE_NEW_LINES);
		foreach ($fai_file_content as $line) 
		{
			// skip empty lines
			if (trim($line) == "") continue;

			// split table
			$parts = explode("\t", $line);
			if (count($parts) != 5)
			{
				trigger_error("Error parsing Fasta index file!", E_USER_ERROR);
			}
			$chr = trim($parts[0]);
			$len = intval($parts[1]);
			$new_contigs[] = "##contig=<ID={$chr},length={$len}>";
		}

		if(empty($new_contigs))
		{
			trigger_error("No new contig lines were written for ".$vcf, E_USER_WARNING);
		}
		else
		{
			array_splice( $new_file_lines, $line_below_reference_info, 0, $new_contigs);  
			file_put_contents($vcf, implode("\n", $new_file_lines));
		}
	}

}

$aa1_to_aa3 = array( 'A'=>"Ala", 'R'=>"Arg", 'N'=>"Asn", 'D'=>"Asp", 'C'=>"Cys", 'E'=>"Glu", 'Q'=>"Gln", 'G'=>"Gly", 'H'=>"His", 'I'=>"Ile", 'L'=>"Leu", 'K'=>"Lys", 'M'=>"Met", 'F'=>"Phe", 'P'=>"Pro", 'S'=>"Ser", 'T'=>"Thr", 'W'=>"Trp", 'Y'=>"Tyr", 'V'=>"Val", '*'=>"Ter");
$aa3_to_aa1 = array_flip($aa1_to_aa3);

//converts amino acid from three letter notation to 1 letter notation
function aa3_to_aa1($three_letter_notation)
{
	return strtr($three_letter_notation, $GLOBALS["aa3_to_aa1"]);
}

function aa1_to_aa3($one_letter_notation)
{
	return strtr($one_letter_notation, $GLOBALS["aa1_to_aa3"]);
}

/**
 * get_related_processed_samples
 * 
 * Find related processed samples.
 *
 * @param  mixed $db
 * @param  string $ps_name
 * @param  string $relation
 * @param  string $systype
 * @param  bool $exclude_bad
 * @return array
 */
function get_related_processed_samples(&$db, $ps_name, $relation, $systype="", $exclude_bad=true)
{
	$ps_name = trim($ps_name);
	list($sample_name, $process_id) = explode("_", $ps_name."_");

	$res = $db->executeQuery("SELECT id, name FROM sample WHERE name=:name", ["name" => $sample_name]);
	if (count($res)!=1)
	{
		return [];
	}
	$sample_id = $res[0]['id'];

	//related samples
	$res = $db->executeQuery("SELECT sample1_id, sample2_id FROM sample_relations WHERE relation=:rel AND (sample1_id=:sid OR sample2_id=:sid)", ["rel" => $relation, "sid" => $sample_id]);
	$related_samples_ids = array_diff(array_merge(array_column($res, "sample1_id"), array_column($res, "sample2_id")), [$sample_id]);
	
	
	//collect processed samples
	$query = <<<SQL
SELECT
	CONCAT(s.name, '_', LPAD(ps.process_id, 2, '0')) as ps_name,
	ps.sample_id, ps.process_id, ps.processing_system_id, ps.quality,
	sys.id, sys.type
FROM
	processed_sample as ps
LEFT JOIN merged_processed_samples mps ON mps.processed_sample_id = ps.id
LEFT JOIN sample s ON s.id = ps.sample_id
LEFT JOIN processing_system sys ON sys.id = ps.processing_system_id
WHERE
	mps.merged_into IS NULL AND
	ps.sample_id=:sid
SQL;

	if ($systype !== "")	$query .= " AND sys.type='{$systype}'";
	if ($exclude_bad)		$query .= " AND ps.quality!='bad'";

	$psamples = [];
	foreach ($related_samples_ids as $rel_sample_id)
	{

		$res = $db->executeQuery($query, ["sid" => $rel_sample_id]);
		$psamples = array_merge($psamples, array_column($res, "ps_name"));
	}

	return $psamples;
}

/**
 * annotate_gsvar_by_gene
 * 
 * Annotate GSvar file based on 'gene' column.
 *
 * @param  object $gsvar reference to GSvar Matrix object
 * @param  string $outfile_f output GSvar file
 * @param  string $annotation_f annotation file
 * @param  string $key column in annotation file to use as key
 * @param  string $column column in annotation file to use as value
 * @param  string $column_name output column name
 * @param  string $column_description output column description
 * @return void
 */
function annotate_gsvar_by_gene(&$gsvar, $annotation_f, $key, $column, $column_name, $column_description, $numeric=true)
{
	$genes = $gsvar->getCol($gsvar->getColumnIndex("gene"));

	$annotation = Matrix::fromTSV($annotation_f);
	$values = array_combine($annotation->getCol($annotation->getColumnIndex($key)), $annotation->getCol($annotation->getColumnIndex($column)));

	foreach($genes as $k => $v)
	{
		$annotated_genes = explode(',', $v);
		$vals = [];
		foreach ($annotated_genes as $g)
		{
			if ($numeric)
			{
				$vals[] = (isset($values[$g]) && $values[$g] != "n/a") ? number_format($values[$g], 4) : "n/a";
			}
			else
			{
				$vals[] = isset($values[$g]) ? $values[$g] : "na";
			}
		}

		$genes[$k] = implode(",", $vals);
	}
	$gsvar->removeColByName($column_name);
	$gsvar->addCol($genes, $column_name, $column_description);
}

//checks that the genome build of a BAM, VCF (small variants or SVs) or TSV (CNVs) matches the expected build.
//If genome cannot be determined, a E_USER_NOTICE is triggered and 0 is returned.
//If the genome matches 1 is returned.
//If the genome does not match a E_USER_ERROR is triggered.
//If @throw_error is false, no error is triggered and -1 is returned. 
function check_genome_build($filename, $build_expected, $throw_error = true)
{
	$builds = array();
	
	//check file exists
	if (!file_exists($filename))
	{
		trigger_error("Cannot check genome build of file '{$filename}'. The file does not exist!",  E_USER_ERROR);
	}
	
	//BAM file
	if (ends_with($filename, ".bam") || ends_with($filename, ".cram"))
	{
		$samtools_command = execApptainer("samtools", "samtools view", "-H $filename", [$filename], [], true);
		list($stdout, $stderr, $exit_code) = exec2("{$samtools_command} | egrep '^@PG' ");
		if  ($exit_code==0)
		{
			foreach($stdout as $line)
			{
				$split_line = explode("\t", trim($line));
				if ($split_line[0] == "@PG")
				{
					if (($split_line[1] == "ID:bwa") || ($split_line[1] == "ID:bwa-mem2"))
					{
						$build = "";
						// parse genome build from bwa command line
						foreach($split_line as $column)
						{
							if (starts_with($column, "CL:"))
							{
								$ref_file_path = explode(" ", $column)[2];
								$build = basename($ref_file_path, ".fa");
								break;
							}
						}
						if ($build!="") 
						{
							$builds[] = $build;
						}
					}
					else if ($split_line[1] == "ID: Hash Table Build")
					{
						$build = "";
						// parse genome build from ABRA2 command line
						foreach($split_line as $column)
						{
							if (starts_with($column, "CL:"))
							{
								$cl = explode(" ", $column);
								$ref_file_path = "";
								for ($i=0; $i < count($cl); ++$i) 
								{ 
									if($cl[$i] == "--ht-reference")
									{
										$ref_file_path = $cl[$i + 1];
										break;
									}
								}
								if (($ref_file_path != "") && basename($ref_file_path, ".fa") != "genome") //special case NovaSeq X: always uses genome.fa as genome file
								{
									$build = basename($ref_file_path, ".fa");
									break;
								}
							}
						}
						if ($build!="") 
						{
							$builds[] = $build;
						}
					}
					else if ($split_line[1] == "ID:STAR")
					{
						$build = "";
						// parse genome build from STAR command line
						foreach($split_line as $column)
						{
							if (starts_with($column, "CL:"))
							{
								while(contains($column, "  ")) $column = strtr($column, ["  "=>" "]);
								$cl = explode(" ", $column);
								for ($i=0; $i < count($cl); ++$i) 
								{ 
									if($cl[$i] == "--genomeDir")
									{
										$path = trim($cl[$i + 1]);
										if (!ends_with($path, "/")) $path .= "/";
										$path_parts = explode("/", $path);
										$build = $path_parts[count($path_parts)-2];
										break;
									}
								}
							}
						}
						if ($build!="") 
						{
							$builds[] = $build;
						}
					}
					else if ($split_line[1] == "ID:minimap2")
					{
						$build = "";
						// parse genome build from minimap2 command line
						foreach($split_line as $column)
						{
							if (starts_with($column, "CL:"))
							{
								while (contains($column, "  ")) $column = strtr($column, ["  "=>" "]);
								$cl = explode(" ", $column);
								//use the first entry that ends with '.fa' when iteration through the parameter list in reverse order (normally second-to-last)
								$idx = count($cl);
								while ($idx)
								{
									$parameter = $cl[--$idx];
									if (ends_with($parameter, ".fa"))
									{
										$build = basename($parameter, ".fa");
										break;
									} 
								}
								if ($build!="") break;								
							}
						}
						if ($build!="") 
						{
							$builds[] = $build;
						}
					}
					else if	($split_line[1] == "ID: DRAGEN SW build")
					{
						//NovaSeq X: Genome has to be stripped from the command line
						$build = "";
						// parse genome build from bwa command line
						foreach($split_line as $column)
						{
							if (starts_with($column, "CL:"))
							{
								$cl = explode(" ", $column);
								for ($i=0; $i < count($cl); ++$i) 
								{
									if($cl[$i] == "--ref-dir")
									{
										$ref_file_path = $cl[$i + 1];
										break;
									}
								}
								if (starts_with($ref_file_path, "/usr/local/illumina/install/genomes/"))
								{
									$build = trim(explode("/", $ref_file_path)[6]);
									break;
								}
							}
						}
						if ($build!="") 
						{
							$builds[] = $build;
						}
					}

				}
			}
		}
	}
	
	//small variants and unannotated structural variants (unannotated)
	if (ends_with($filename, ".vcf.gz") || ends_with($filename, ".vcf"))
	{
		if (ends_with($filename, ".vcf.gz"))
		{
			list($stdout, $stderr, $exit_code) = exec2("zcat $filename | egrep '^##reference='");
		}
		else
		{
			list($stdout, $stderr, $exit_code) = exec2("egrep '^##reference=' {$filename}");
		}
		
		if ($exit_code==0)
		{
			foreach($stdout as $line)
			{
				list(, $fasta) = explode("=", $line);
				if (ends_with($fasta, ".fa"))
				{
					$builds[] = basename($fasta, ".fa");
				}
				else if (contains($fasta, "/dragen/")) //special handling for Dragen (e.g. file://staging/genomes/GRCh38/dragen/reference.bin)
				{
					$fasta = strtr($fasta, ["//"=>"/"]);
					$parts = explode("/", $fasta);
					$builds[] = $parts[count($parts)-3];
				}
			}
		}
	}
		
	//GSvar
	if (ends_with($filename, ".GSvar"))
	{
		list($stdout, $stderr, $exit_code) = exec2("egrep '^##GENOME_BUILD=' $filename", false);
		if ($exit_code==0)
		{
			foreach($stdout as $line)
			{
				list(, $name) = explode("=", $line);
				$builds[] = $name;
			}
		}
	}
	
	//CNVs file
	if (ends_with($filename, ".tsv"))
	{
		list($stdout, $stderr, $exit_code) = exec2("egrep '^##GENOME_BUILD=' $filename", false);
		if ($exit_code==0)
		{
			foreach($stdout as $line)
			{
				list(, $name) = explode("=", $line);
				$builds[] = $name;
			}
		}
	}
	
	//structural variants (annotated)
	if (ends_with($filename, ".bedpe"))
	{
		list($stdout, $stderr, $exit_code) = exec2("cat $filename | egrep '^##reference='", false);
		if ($exit_code==0)
		{
			foreach($stdout as $line)
			{
				list(, $fasta) = explode("=", $line);
				
				if (ends_with($fasta, ".fa"))
				{
					$builds[] = basename($fasta, ".fa");
				}
				else if (contains($fasta, "/dragen/") || contains($fasta, "/DRAGEN/")) //special handling for Dragen (e.g. file:///staging/human/reference/GRCh38/dragen/ or file:///usr/local/illumina/install/genomes/GRCh38/DRAGEN/10)
				{
					$fasta = strtr($fasta, ["//"=>"/"]);
					$parts = explode("/", $fasta);
					$builds[] = $parts[count($parts)-3];
				}
			}
		}
	}
	
	//check that there is not more/less than one genome match
	$builds = array_map('strtolower', $builds);
	$builds = array_unique($builds);
	if (count($builds) < 1) 
	{
		trigger_error("File '$filename' does not contain genome build information!", E_USER_NOTICE);
		return 0;
	}
	if (count($builds) > 1)
	{
		trigger_error("File '$filename' contains inconsistent genome build information!", E_USER_NOTICE);
		return 0;
	}
	
	//compare found and expected build
	$build_found = $builds[0];
	$build_expected = strtolower($build_expected);
	if (!starts_with($build_found, $build_expected))
	{
		if ($throw_error)
		{
			trigger_error("Genome build of file '{$filename}' is '{$build_found}'. It does not match expected genome build '{$build_expected}'!",  E_USER_ERROR);
		}
		return -1;
	}
	
	return 1;
}

//converts a proessing system genome build string to a string compatible with ngs-bits
function ngsbits_build($system_build)
{
	$system_build = strtolower($system_build);	
	if ($system_build=="hg38" || $system_build=="grch38")
	{
		return "hg38";
	}
	else if (starts_with($system_build, "hg38_") || starts_with($system_build, "grch38"))
	{
		return "hg38";
	}
	else if ($system_build=="hg19" || $system_build=="grch37")
	{
		return "hg19";
	}
	else if (starts_with($system_build, "hg19_") || starts_with($system_build, "grch37_"))
	{
		return "hg19";
	}
	else 
	{
		return "non_human";
	}
}

//returns an array of processed sample names that are currently being analyzed via the SGE queue
function ps_running_in_sge()
{
	//extract job ids
	$job_ids = [];
	list($tmp) = exec2("qstat -u '*'");
	foreach($tmp as $line)
	{
		$line = trim($line);
		if ($line=="") continue;
		$id = explode(" ", $line)[0];
		if (is_numeric($id))
		{
			$job_ids[] = $id;
		}
	}
	
	$output = [];
	foreach($job_ids as $job_id)
	{
		list($args) = exec2("qstat -j {$job_id} | tr ',' '\n'");
		foreach($args as $arg)
		{
			$arg = trim($arg);
			
			if (contains($arg, "Sample_")) //sample folder
			{
				$parts = explode("/", $arg);
				foreach($parts as $part)
				{
					if (!contains($part, "Sample_")) continue;
					
					$ps = strtr($part, ["Sample_"=>""]);
					$output[$ps] = true;
				}
			}
			else if (preg_match("/_[0-9][0-9]*$/", $arg)) //processed sample name
			{
				$output[$arg] = true;
			}
		}
	}
	
	return array_keys($output);
}

function bed_size($filename)
{
	$container_command = execApptainer("ngs-bits", "BedInfo", "-in $filename", [$filename], [], true);
	list($stdout) = exec2("$container_command | grep -i bases");
	return trim(explode(":", $stdout[0])[1]);
}

function phenotype_roi(&$db_conn, $ps_name)
{		
	$info = get_processed_sample_info($db_conn, $ps_name);
	$s_id = $info['s_id'];
	$hpos = $db_conn->getValues("SELECT disease_info FROM sample_disease_info WHERE sample_id='$s_id' AND `type`='HPO term id'");
	
	//store HPO term ids
	$tmp_hpos = temp_file(".txt");
	file_put_contents($tmp_hpos, $hpos);
	
	//convert HPO terms to gene symbols
	$tmp_genes = temp_file(".tsv");
	exec2("PhenotypesToGenes -in $tmp_hpos -ignore_invalid -ignore_non_phenotype | cut -f1 > $tmp_genes");
	
	//convert genes to roi
	$output = [];
	list($stdout) = exec2("GenesToBed -in $tmp_genes -source ensembl -mode gene -fallback | BedMerge");
	foreach($stdout as $line)
	{
		$line = trim($line);
		if ($line=="" || $line[0]=='#') continue;
		
		list($chr, $start, $end) = explode("\t", $line);
		
		$output[] = [$chr, $start, $end];
	}
	
	return $output;
}


function phenotype_roi_overlaps(&$roi, $chr, $start, $end)
{
	foreach($roi as list($chr2, $start2, $end2))
	{
		if ($chr==$chr2 && range_overlap($start, $end, $start2, $end2))
		{
			return true;
		}
	}
	
	return false;
}

//Returns the base file name without path and extension. Used e.g. to extract the base name independent of BAM/CRAM format.
function basename2($filename)
{
	$filename = basename($filename);
	
	//not file extension => return full name
	$sep_idx = strrpos($filename, '.');
	if ($sep_idx===FALSE) return $filename;
	
	return substr($filename, 0, $sep_idx);
}

//If the given file is a CRAM, a temporary BAM is created and the path to it is returned. Otherwise returns the given filename.
function convert_to_bam_if_cram($filename, $parser, $build, $threads, $tmp_folder="")
{
	if (ends_with(strtolower($filename), ".cram"))
	{
		if ($tmp_folder=="") $tmp_folder = $parser->tempFolder();
		$bam = "{$tmp_folder}/".substr(basename($filename), 0, -5).".bam";
		$parser->execTool("Tools/cram_to_bam.php", "-cram $filename -bam $bam -build $build -threads $threads");
		return $bam;
	}
	
	return $filename;
}

//check if run is created by NovaSeq X (Plus)
function is_novaseq_x_run($run_parameters_xml)
{
	$xml = simplexml_load_file($run_parameters_xml);
	if ($xml===false) trigger_error("Could not load XML file: $run_parameters_xml", E_USER_ERROR);
	if(empty($xml->InstrumentType)) return false; //No entry => e.g. NovaSeq 6000
	$instrument_type = $xml->InstrumentType;
	if(($instrument_type == "NovaSeqXPlus") || ($instrument_type == "NovaSeqX")) return true;
	//else
	return false; //unknown Sequencer
}

//checks if the genome used to map the BAM/CRAM file had masked false duplications
function genome_masked($bam, $build)
{
	$genome = genome_fasta($build);
	$samtools_command = execApptainer("samtools", "samtools view", "-T {$genome} {$bam} chr21:6110084-6124379", [$genome, $bam], [], true);
	list($stdout) = exec2("{$samtools_command} | wc -l", false);
	
	$read_count = trim(implode("", $stdout));
	
	return $read_count==0;
}

//check if BAM file contains methylation data (only check the first $n_rows)
function contains_methylation($bam_file, $n_rows=100, $build="GRCh38")
{
	if (!file_exists($bam_file)) trigger_error("BAM file '{$bam_file}'", E_USER_ERROR);
	$genome = genome_fasta($build);
	// ignore errors occuring of unknown reason (broken pipe)
	$samtools_command = execApptainer("samtools", "samtools view", "-T {$genome} {$bam_file}", [$genome, $bam_file], [], true);
	list($stdout) = exec2("{$samtools_command} | head -n {$n_rows}", false);
	//additional testing since we cannot rely on samtools error reporting
	if (count($stdout) != $n_rows) trigger_error("Couldn't extract the first {$n_rows} rows of the BAM file!", E_USER_ERROR);

	$n_mm = 0;
	$n_ml = 0;
	foreach ($stdout as $row) 
	{
		$columns = explode("\t", $row);
		foreach ($columns as $cell) 
		{
			if(starts_with($cell, "MM:")) $n_mm++;
			elseif(starts_with($cell, "ML:")) $n_ml++;
		}
	}

	//sanity checks
	if($n_mm != $n_ml) trigger_error("Number of ML tags and MM tags is not equal!\nMM:\t{$n_mm}\nML:\t{$n_ml}", E_USER_ERROR);

	// no methylation
	if($n_mm == 0) return false;
	
	// methylation
	if($n_mm == $n_rows) return true;

	//else: something is wrong
	trigger_error("Ambiguous tag counts. Please check BAM file!\nMM:\t{$n_mm}/{$n_rows}\nML:\t{$n_ml}/{$n_rows}", E_USER_ERROR);
}

//extracts base-calling description from BAM header
function get_read_group_description($bam_file)
{
	$rg_description = array();
	$samtools_command = execApptainer("samtools", "samtools view", "-H $bam_file", [$bam_file], [], true);
	list($stdout, $stderr, $exit_code) = exec2("{$samtools_command} | egrep '^@RG' ", false);
	if  ($exit_code==0 || $exit_code==1)	{
		foreach($stdout as $line)
		{
			$split_line = explode("\t", trim($line));
			if ($split_line[0] == "@RG")
			{
				foreach($split_line as $column)
				{
					if (starts_with($column, "DS")) $rg_description[] = $column;
				}
			}
		}
	}
	else
	{
		trigger_error("ERROR: Extraction of ReadGroup description failed! \n\tExit code: {$exit_code} \n\t".implode("\n\t", $stderr), E_USER_ERROR);
	}

	if (count($rg_description) < 1) trigger_error("WARNING: No ReadGroup description found!", E_USER_WARNING);

	return $rg_description;
}

//returns the model which was used for base-calling
function get_basecall_model($bam_file)
{
	$basecall_model = array();
	$samtools_command = execApptainer("samtools", "samtools view", "-H $bam_file", [$bam_file], [], true);
	list($stdout, $stderr, $exit_code) = exec2("{$samtools_command} | egrep '^@RG' ", false);
	if  ($exit_code==0 || $exit_code==1)	{
		foreach($stdout as $line)
		{
			$split_line = explode("\t", trim($line));
			if ($split_line[0] == "@RG")
			{
				foreach($split_line as $column)
				{
					if (starts_with($column, "DS:"))
					{
						$entries = explode(" ", substr($column, 3));
						foreach($entries as $entry)
						{
							if (starts_with($entry, "basecall_model=")) $basecall_model[] = trim(explode("=", $entry)[1]);
						}
					}
				}
			}
		}
	}
	else
	{
		trigger_error("ERROR: Extraction of basecall model failed! \n\tExit code: {$exit_code} \n\t".implode("\n\t", $stderr), E_USER_ERROR);
	}

	//remove duplicates
	$basecall_model = array_unique($basecall_model);

	if (count($basecall_model) < 1)
	{
		return "";
	}
	elseif  (count($basecall_model) > 1)
	{
		trigger_error("Multiple different basecall model entries found! Using the first one ('".$basecall_model[0]."').", E_USER_WARNING);
	}
	return $basecall_model[0];

}

//returns the project path, or NULL if the project does not exist in NGSD
function get_project_folder(&$db, $project_name)
{
	//get project id
	$project_name = trim($project_name);
	$project_id = $db->getValue("SELECT id FROM project WHERE name LIKE '$project_name'", -1);
	if ($project_id==-1) return NULL;
	
	//folder override
	$folder_override = $db->getValue("SELECT folder_override FROM project WHERE id='{$project_id}'");
	$folder_override = trim($folder_override);
	if ($folder_override!="")
	{
		$output = $folder_override;
	}
	else //fallback to default path
	{
		$type = $db->getValue("SELECT type FROM project WHERE id='{$project_id}'");
		$output = get_path('project_folder')[$type];
		if (!ends_with($output, "/")) $output .= "/";
		$output .= $project_name;
	}
	
	//output
	$output = strtr($output, ["//"=>"/"]);
	if (!ends_with($output, "/")) $output .= "/";
	return $output;
}

//updates the GSvar SAMPLE entry with NGSD info (inplace)
function update_gsvar_sample_header($file_name, $disease_status_map)
{
	$file_content = Matrix::fromTSV($file_name);
	$new_comments = [];
	foreach ($file_content->getComments() as $line) 
	{
		//only update SAMPLE lines
		if(!starts_with($line, "#SAMPLE="))
		{
			$new_comments[] = $line;
			continue;
		}
		
		//determine overrides
		$override_map = [];
		list($ps) = explode(",", substr($line, 12, -1), 2);			
		if (isset($disease_status_map[$ps]))
		{
			$override_map["DiseaseStatus"] = $disease_status_map[$ps];
			
			//keep gender (for trios)
			if (contains($line, "Gender=male")) $override_map["Gender"] = "male";
			if (contains($line, "Gender=female")) $override_map["Gender"] = "female";
		}

		//generate new header line
		$new_comments[] = gsvar_sample_header($ps, $override_map, "#", "\n"); 
	}
	$file_content->setComments($new_comments);
	$file_content->toTSV($file_name);
}

//get bam read count
function get_read_count($bam_file, $threads = 4, $samtools_params = array(), $build = "GRCh38", $region = "")
{
	$read_count = -1;
	list($stdout, $stderr, $exit_code) = execApptainer("samtools", "samtools view", "-@ {$threads} -c -T ".genome_fasta($build)." ".implode($samtools_params)." {$bam_file} {$region}", [$bam_file, genome_fasta($build)]);
	if ($exit_code == 0)
	{
		$read_count = (int) $stdout[0];
	}
	else
	{
		trigger_error("Error calculating read counts: \n".implode("\n", $stderr), E_USER_ERROR);
	}

	return $read_count;
}

//compare bam read count
function compare_bam_read_count($bam_file1, $bam_file2, $threads = 4, $throw_error = true, $chr_wise = false, $relative_tolerance = 0.0, $samtools_params = array(), $build = "GRCh38")
{
	if ($chr_wise)
	{
		$differences = array();
		foreach (chr_list() as $chr) 
		{
			$read_count1 = get_read_count($bam_file1, $threads, $samtools_params, $build, $chr);
			$read_count2 = get_read_count($bam_file2, $threads, $samtools_params, $build, $chr);
			$diff = abs($read_count1 - $read_count2);
			//count are identical
			if ($diff == 0) continue;

			$relative_error = $diff / (($read_count1 + $read_count2)/2);
			if ($relative_error <= $relative_tolerance) continue;

			$differences[] = "\t{$chr}:\t {$read_count1} vs. {$read_count2} reads";
		}

		if (count($differences) == 0) return true;

		if ($throw_error) trigger_error("Bam file read count do not match!\n".implode("\n", $differences), E_USER_ERROR);
		return false;
	}
	else
	{
		$read_count1 = get_read_count($bam_file1, $threads, $samtools_params, $build);
		$read_count2 = get_read_count($bam_file2, $threads, $samtools_params, $build);
		$diff = abs($read_count1 - $read_count2);
		//count are identical
		if ($diff == 0) return true;

 		$relative_error = $diff / (($read_count1 + $read_count2)/2);
		if ($relative_error <= $relative_tolerance) return true;

		if ($throw_error) trigger_error("Bam file read count do not match! ({$read_count1} vs. {$read_count2} reads)", E_USER_ERROR);
		return false;
	}
}

/**
    @brief Determines the sequencing platform type for a processing system.
    @param string $name_short Short name of the processing system.
    @return string 'PB', 'ONT'.
*/
function get_longread_sequencing_platform($name_short)
{
	//TODO: use better way than short name comparison
    if (starts_with($name_short, 'LR-PB-')) {
        return 'PB';
    } else {
        return 'ONT';
    }
}

/**
	@brief	Determines the number of reads of a set of ORA files
	@param	array containing ORA file names
	@return	int number of reads
*/
function get_ora_read_count($ora_files)
{
	$read_count = 0;
	foreach ($ora_files as $ora_file) 
	{
		$rc_current_file = -1;
		list($stdout, $stderr, $exit_code) = execApptainer("orad", "orad", "-i {$ora_file}", [$ora_file]);
		foreach ($stdout as $line) 
		{
			if (starts_with($line, "Total number of sequences"))
			{
				$rc_current_file = intval(explode(":", $line)[1]);
				break;
			}
		}
		if ($rc_current_file < 0) trigger_error("Sequences count not found in ORA stats of file '{$ora_file}'!", E_USER_ERROR);
		$read_count += $rc_current_file;
	}
	return $read_count;
}

/**
	@brief	Determines the number of reads of a set of FASTQ files
	@param	array containing FASTQ file names
	@return	int number of reads
*/
function get_fastq_read_count($fastq_files)
{
	$read_count = 0;
	foreach ($fastq_files as $fastq_file) 
	{
		$rc_current_file = -1;
		list($stdout, $stderr, $exit_code) = exec2("zcat -f {$fastq_file} | wc -l");
		$rc_current_file = intval($stdout[0]) / 4;
		$read_count += $rc_current_file;
	}
	return $read_count;
}

//Checks if a BED file is sorted
function bed_is_sorted($filename)
{
	$is_sorted = "";
	list($stdout) = execApptainer("ngs-bits", "BedInfo", "-in {$filename} -nomerge", [$filename]);
	foreach($stdout as $line)
	{
		$line = trim($line);
		if (starts_with($line, "Is sorted"))
		{
			$is_sorted = trim(explode(":", $line, 2)[1]);
		}
	}
	if ($is_sorted=="yes") return true;
	if ($is_sorted=="no") return false;
	trigger_error("Could not determine if file is sorted: $filename", E_USER_ERROR);
}

?>
