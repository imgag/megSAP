<?php
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

$parser = new ToolBase("cgi_send_data", "Sends annotation file to CGI.");
$parser->addInfile("mutations", "Input .vcf-file with mutation annotations, hg19 coordinates.", true);
$parser->addInfile("cnas","File containing the CNV data. ",true);
$parser->addInfile("translocations","File containing the translocation data",true);
$parser->addString("cancertype", "cancer type, see cancergenomeinterpreter.org for nomenclature",true,"CANCER");
$parser->addString("out", "Output file for zipped CGI result file",false);
$parser->addInfile("t_region", ".txt-File which contains genes that shall be included in CNV analysis. If unset, all genes in the CNV file will be sent to CGI", true);
$parser->addFlag("no_del", "Do not delete Job on Cancer Genome Interpreter after submission");
$parser->addFlag("no_snv_limit","Allow uploading more than 1000 SNPs");
extract($parser->parse($argv));

//discard if no input files given
if(!isset($mutations) && !isset($cnas) && !isset($translocations))
{
	$parser->printUsage();
	trigger_error("At least one variant file must be uploaded to CGI",E_USER_ERROR);
}

//set user credentials
$user =  get_path("cgi_user");
$token = get_path("cgi_token");
$url =  get_path("cgi_url");
$header = "Authorization: ".$user." ".$token;


//Creates header for cURL requests
function generateHeader($user,$token)
{
	return "Authorization: ".$user." ".$token;
}

//sends annotation files to CGI and returns job identifier, mutation_file is in vcf format, cnv_file in CGI specific format
function sendData($title,$cancer_type = "",$mutation_file = "",$cnv_file = "",$translocation_file = "")
{
	global $token,$user,$url,$header,$parser;
	if($cancer_type == "")
	{
		$cancer_type = "CANCER";
	}
	if($mutation_file == "" && $cnv_file == "" && $translocation_file == "")
	{
		trigger_error("At least one sample file containing data should be provided. Cannot send data to CGI.", E_USER_ERROR);
		exit(1);
	}

	$parameters = "-k --silent --show-error --request POST --url $url --header \"$header\" ";
	if($mutation_file !== "")
	{
		$parameters = $parameters. "-F \"mutations=@$mutation_file\" ";
	}
	if($cnv_file !== "")
	{
		$parameters = $parameters."-F \"cnas=@$cnv_file\" ";
	}
	if($translocation_file !== "")
	{
		$parameters = $parameters."-F \"translocations=@$translocation_file\" ";
	}
	$parameters = $parameters."-F \"cancer_type=$cancer_type\" ";
	$parameters = $parameters."-F \"title=$title\" ";
	
	$result = $parser->exec("curl",$parameters,true,false);

	//Exception handling for connection errors
	$error_code = $result[2];
	if($error_code != 0)
	{
		print("Error while sending data to CGI. CURL return error code ".$error_code."\n");
		
		//retry three times
		for($i=0;$i<3;$i++)
		{
			print("Retrying\n");
			sleep(20);
			$result = $parser->exec("curl",$parameters,true,false);
			$error_code = $result[2];
			
			//finish loop if retry was successful, end program 
			if($error_code == 0)
			{
				print "Success. Data transfered successfully to CancerGenomeInterpreter.org\n";
				break;
			}
		}
		if($error_code != 0)
		{
			trigger_error("Error while sending data to CGI. CURL return error code ".$error_code,E_USER_ERROR);
		}
	}
	
	//array with the answer received by CGI
	$cgi_answer = (array)json_decode($result[0][0],true);
	
	//Exception handling for CGI errors
	$error_in_job = false;
	if(array_key_exists("error_code",$cgi_answer)) 
	{
		$error_in_job = true;
		switch($cgi_answer["error_code"])
		{
			//Error code 403 is returned if all slots are full -> check for finished jobs which can be removed / wait 
			case "403":			
				$identifiers = getIdentifiers($user,$token);
				//Delete finished jobs
				foreach($identifiers as $jobID)
				{
					$jobStatus = getJobStatus($jobID,$user,$token,$url);
					//Delete jobs which are marked as finsihed ("Done") or failed ("Error")
					if($jobStatus == "Done" || $jobStatus == "Error")
					{
						print("CancerGenomeInterpreter.org queue is full. Deleting job ".$jobID.".\n");
						//wait 5 seconds to make sure job data can be downloaded by a potential other job
						sleep(5);
						//Check job status a second time to make sure job has not yet deleted by another jobID
						$jobStatus = getJobStatus($jobID,$user,$token,$url);
						if($jobStatus == "Done" || $jobStatus == "Error")
						{
							delJob($jobID,$user,$token,$url);
						}
					}
				}
				break;
				
			default:
				trigger_error("CancerGenomeInterpreter.org returned an error. Error code: ".$cgi_answer["error_code"],E_USER_WARNING);
		}
		
		//Try again in case of previous error
		if($error_in_job)
		{
			print("\n Retrying to send data\n");
			$result = $parser->exec("curl",$parameters,true);
			$cgi_answer = (array)json_decode($result[0][0],true);
			if(array_key_exists("error_code",$cgi_answer))
			{
				trigger_error("CancerGenomeInterpreter.org returned an error. Error code: ".$cgi_answer["error_code"],E_USER_ERROR);
			}
		}
	}
	//returns job identifier as array
	return json_decode($result[0][0],true);
}

//Request identifiers of all jobs, return identifiers as array
function getIdentifiers($user,$token,$url='https://www.cancergenomeinterpreter.org/api/v1')
{
	global $parser;
	$header = generateHeader($user,$token);
	$parameters = "-k --silent --show-error --request GET --url $url --header \"$header\" ";
	
	$result = $parser->exec("curl",$parameters,true);
	
	//Transform tokens to array
	$tokens = json_decode($result[0][0],true);
	
	return $tokens;
}

//Requests basic information about a job in CGI account
function getJobInformation($job_id,$user,$token,$url='https://www.cancergenomeinterpreter.org/api/v1')
{
	global $parser;
	$header = generateHeader($user,$token);
	$parameters = "-k --silent --show-error --request GET --url $url/".$job_id." --header \"$header\" ";
	$result = $parser->exec("curl",$parameters,true);
	
	//return information as array
	return json_decode($result[0][0],true);
}

//get job status
function getJobStatus($job_id,$user,$token,$url='https://www.cancergenomeinterpreter.org/api/v1')
{
	global $parser;
	$header = generateHeader($user,$token);	
	$parameters = "-k --silent --show-error --request GET --url $url/".$job_id." --header \"$header\"  -G --data \"action=logs\"";
	
	$result = $parser->exec("curl",$parameters,false);

	//$results[0][0] contains array with processing status
	return json_decode($result[0][0],true)["status"];
}

//Download results and store them as zip-files
function downloadJobResults($destination,$job_id,$user,$token,$url='https://www.cancergenomeinterpreter.org/api/v1')
{
	global $parser;
	$header = generateHeader($user,$token);	
	$parameters = "-k --silent --show-error --request GET --url $url/".$job_id." --header \"$header\"  -G --data \"action=download\" ";
	$parameters = $parameters."-o $destination";
	
	$result = $parser->exec("curl",$parameters,true);
}

//delete a job
function delJob($job_id,$user,$token,$url='https://www.cancergenomeinterpreter.org/api/v1')
{
	global $parser;
	$header = generateHeader($user,$token);
	$parameters = "-k --silent --show-error --request DELETE --url $url/".$job_id." --header \"$header\" ";	
	$result = $parser->exec("curl",$parameters,true);
}

//returns associative array (key=>value) with genes => alteration type
function get_genes_from_file($file)
{
	$cnvs = Matrix::fromTSV($file);
	$i_genes = $cnvs->getColumnIndex("genes");
	
	$i_copy_numbers = $cnvs->getColumnIndex("region_copy_numbers",false,false); // CNVHunter file
	$i_cn_change = $cnvs->getColumnIndex("CN_change",false,false); //ClinCnv file
	
	if($i_copy_numbers === false && $i_cn_change === false) trigger_error("Unknown CNV format. Aborting.",E_USER_ERROR);
	
	$genes_type = array(); // ass. array with $genes => CNV_type
	
	if($i_copy_numbers !== false) //CNVHunter file
	{
		for($row=0;$row<$cnvs->rows();++$row)
		{
			$tmp_region_copy_numbers = median(explode(',',$cnvs->get($row,$i_copy_numbers)));
			
			$genes = $cnvs->get($row,$i_genes);
			
			if($tmp_region_copy_numbers>2.) $genes_type[$genes] = "AMP";
			else $genes_type[$genes] = "DEL";
		}
	} else //ClinCNV
	{
		for($row=0;$row<$cnvs->rows();++$row)
		{
			$genes = $cnvs->get($row,$i_genes);
			$cn_change = $cnvs->get($row,$i_cn_change);
			if($cn_change > 2.) $genes_type[$genes] = "AMP";
			elseif($cn_change < 2.) $genes_type[$genes] = "DEL";
			else $genes_type[$genes] = "NA";
		}
	}
	
	//explode genes
	$output = array();
	foreach($genes_type as $genes_per_cnv => $type)
	{
		$genes = explode(',',$genes_per_cnv);
		foreach($genes as $gene)
		{
			//skip genes that are annotated in CNVs with inconsistent CN type
			if(array_key_exists($gene,$output) && $type != $output[$gene]) 
			{
				$output[$gene] = "NA";
				continue;
			}
			$output[$gene] = $type;
		}
	}
	
	return $output;
}

//read CNV .tsv-file and transform it into CGI specific format. If target region is set: remove genes that do not lie inside the target region
function transform_cnv_annotations($tsv_in_file,$tsv_out_filename)
{
	global $t_region;
	$genes = get_genes_from_file($tsv_in_file);
	$cnv_to = new Matrix();
	$cnv_to->addRow(array("gene","cna"));
	foreach($genes as $gene => $cn_type)
	{
		$cnv_to->addRow(array($gene,$cn_type));
	}

	//transform all outdated gene names
	$approved_gene_names = approve_gene_names($cnv_to->getCol(0));
	if($approved_gene_names[0] == "GENE") $approved_gene_names[0] = "gene";
	$cnv_to->setCol(0,$approved_gene_names);
	//there can be duplicates in CNV list if a gene lies in two reported regions -> remove
	$cnv_to->unique();
	//discard all genes which do not lie in target region
	if(isset($t_region))
	{
		$genes_in_region = explode("\n",file_get_contents($t_region));
		if(count($genes_in_region) == 0)
		{
			trigger_error("Target region file $t_region does not contain any genes.",E_USER_ERROR);
		}
		
		//approve gene names in target region
		$genes_in_region = approve_gene_names($genes_in_region);
		
		//Only use genes which are listed in target region
		$index_discarded_genes = array();
		for($i=0;$i<$cnv_to->rows();$i++)
		{
			if(!in_array($cnv_to->get($i,0),$genes_in_region))
			{
				$index_discarded_genes[] = $i;
			}
		}

		//Remove genes that lie outside target region, begin at the end of cnv_to
		for($i=count($index_discarded_genes)-1;$i>0;$i--)
		{
			$cnv_to->removeRow($index_discarded_genes[$i]);
		}
		
		if($cnv_to->rows() <= 1)
		{
			trigger_error("Warning: CNV file does not contain any genes after filtering for target region.\n",E_USER_WARNING);
		}
	
		global $cnas;
		print(count($index_discarded_genes) . " genes from ".$cnas." were discarded because they do not lie in target region $t_region\n");
	}
	//write to temporary file
	$cnv_to->toTSV($tsv_out_filename);
}

/********
 * MAIN *
 ********/

//Parse SNV file
$temp_mut_file = "";
if(isset($mutations))
{
	if(strpos($mutations,".gz") !== false)
	{
		exec2("gzip -d -k -f  $mutations");
		$mutations = str_replace(".gz","",$mutations);
	}

	$vcf_in = file($mutations, FILE_IGNORE_NEW_LINES);
	$temp_mut_cont = array("sample\tchr\tpos\tref\talt");
	foreach($vcf_in as $line)
	{
		if(starts_with($line,"#")) continue;
		list($chr,$pos,,$ref,$alt) = explode("\t",$line);
		
		//unique identifier for each variant (needed to allocate drugs precisely to a SNPs)
		$temp_id = "{$chr}_{$pos}_{$ref}_{$alt}";
		
		$temp_mut_cont[] = "{$temp_id}\t{$chr}\t{$pos}\t{$ref}\t{$alt}";
	}
	
	if(!$no_snv_limit && count($temp_mut_cont) > 2000)
	{
		trigger_error("Too many variants in {$mutations}. Please filter SNPs before passing them to this tool.",E_USER_ERROR);
	}
	$temp_mut_file = $parser->tempFile(".tsv");
	file_put_contents($temp_mut_file,implode("\n",$temp_mut_cont));
}


//convert CNVs to CGI format
$temp_cnv_file = $parser->tempFile(".tsv");
if(isset($cnas))
{
	transform_cnv_annotations($cnas,$temp_cnv_file);
} 
else 
{
	$temp_cnv_file = "";
}

//check whether temp_cnv_file contains data
if($temp_cnv_file != "" && count(file($temp_cnv_file)) <= 1)
{
	$temp_cnv_file = "";
}

/****************************
 * PARSE TRANSLOCATION FILE *
 ****************************/
$temp_translocations_file = "";
if(isset($translocations))
{
	$in = Matrix::fromTSV($translocations);
	
	$i_genes = $in->getColumnIndex("genes");
	$i_mate_genes = $in->getColumnIndex("mate_genes");
	
	$i_type = $in->getColumnIndex("type");
	
	$mate_pairs = array("fus");
	for($i = 0;$i<$in->rows();++$i)
	{
		$line = $in->getRow($i);
		
		if(trim($line[$i_type]) != "BND") continue;
		$genes = explode(",",trim($line[$i_genes]));
		$mate_genes = explode(",",trim($line[$i_mate_genes]));
		
		if(!empty($genes) && !empty($mate_genes) )
		{
			foreach($genes as $gene)
			{
				if($gene == "") continue;
				foreach($mate_genes as $mate_gene)
				{
					if($mate_gene == "") continue;
					$mate_pairs[] = "{$gene}__{$mate_gene}";
				}
			}
		}
		
	}
	
	if(count($mate_pairs) > 1)
	{
		$temp_translocations_file = $parser->tempFile(".tsv");	
		file_put_contents($temp_translocations_file,implode("\n",$mate_pairs));
	}
}

//Send data to CGI, save job_id 
$jobId = sendData("IMGAG",$cancertype,$temp_mut_file,$temp_cnv_file,$temp_translocations_file);
//display job status
do
{
	$status = getJobStatus($jobId,$user,$token);
	if($status == "Error")
	{
		trigger_error("CGI returned an error. See logfile on CancerGenomeInterpreter.org for Details.", E_USER_ERROR);
	}
	print("status of CGI job ". $jobId . ": " . $status . "\n");

	if($status != "Done")
	{
		sleep(10);
	}
	
}while($status != "Done");

//download CGI results
downloadJobResults($out,$jobId,$user,$token,$url);

//Delete job from remote server
if(!$no_del)
{
	delJob($jobId,$user,$token,$url);
}
?>
