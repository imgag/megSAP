<?php
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

$parser = new ToolBase("cgi_send_data", "Sends annotation file to CGI.");
$parser->addString("cancertype", "cancer type, see cancergenomeinterpreter.org for nomenclature",true);
$parser->addInfile("mutations", "Input .vcf-file with mutation annotations, hg19 coordinates.", true);
$parser->addInfile("cnas","File containing the CNV data. ",true);
$parser->addInfile("translocations","File containing the translocation data",true);
$parser->addString("o_folder", "Output folder for CGI files. If not set, files are stored in the same folder as the input files.",true);
$parser->addInfile("t_region", ".txt-File which contains the gene names of the target region, only these genes will be sent to CGI. If not set all genes will be sent to CGI", true);
$parser->addFlag("no_del", "Do not delete Job in Cancer Genome Interpreter after submission");
extract($parser->parse($argv));

//discard if no input file given
if(!isset($mutations) && !isset($cnas) && !isset($translocations))
{
	trigger_error("At least one file must be uploaded to CGI",E_USER_ERROR);
	exit(1);
}

//determine destination folder 
$destination_folder = "";
if(isset($o_folder))
{

	$destination_folder = $o_folder;
}
else
{
	if(isset($mutations))
	{
		$destination_folder = realpath(dirname($mutations));
	}elseif(isset($cnas)) {
		$destination_folder = realpath(dirname($cnas));
	}elseif(isset($translocations)){
		$destination_folder = realpath(dirname($translocations));
	}
}

//if cancertype not set, use generic CGI type
if(!isset($cancertype))
{
	$cancertype = "CANCER";
}

//set user credentials
$user =  get_path("cgi_user");
$token = get_path("cgi_token");
$url =  get_path("cgi_url");
$header = "Authorization: ".$user." ".$token;

//get Sample IDs

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

//results are stored as zip-files: ToDo: Error Check, if authorization fails, the error message is written to the destination file
function downloadJobResults($destination,$job_id,$user,$token,$url='https://www.cancergenomeinterpreter.org/api/v1')
{
	global $parser;
	$header = generateHeader($user,$token);	
	$parameters = "-k --silent --show-error --request GET --url $url/".$job_id." --header \"$header\"  -G --data \"action=download\" ";
	$parameters = $parameters."-o $destination";
	
	$result = $parser->exec("curl",$parameters,true);
}


//delete job
function delJob($job_id,$user,$token,$url='https://www.cancergenomeinterpreter.org/api/v1')
{
	global $parser;
	$header = generateHeader($user,$token);
	$parameters = "-k --silent --show-error --request DELETE --url $url/".$job_id." --header \"$header\" ";	
	$result = $parser->exec("curl",$parameters,true);
}

//read CNV .tsv-file and transform it into CGI specific format. if target region is set -> remove genes that do not lie inside the target region
function transform_cnv_annotations($tsv_in_file,$tsv_out_filename)
{
	global $t_region;
	
	$cnvs = Matrix::fromTSV($tsv_in_file);
	$i_genes = $cnvs->getColumnIndex("genes");
	$i_copy_numbers = $cnvs->getColumnIndex("region_copy_numbers");
	
	$genes = array();
	$cnv_type = array();

	for($row=0;$row<$cnvs->rows();$row++)
	{
		$tmp_region_copy_numbers = median(explode(',',$cnvs->get($row,$i_copy_numbers)));
		if($tmp_region_copy_numbers>2)
		{
			$cnv_type[] = "AMP";
		}
		else
		{
			$cnv_type[] = "DEL";
		}
		
		$genes[] = $cnvs->get($row,$i_genes);
	}
	
	$cnv_to = new Matrix();
	
	$cnv_to->addRow(array("gene","cna"));	
	
	for($row=0;$row<$cnvs->rows();$row++)
	{
		$tmp_region_genes = explode(',',$genes[$row]);
		for($i=0;$i<count($tmp_region_genes);$i++)
		{
			$cnv_to->addRow(array($tmp_region_genes[$i],$cnv_type[$row]));
		}
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

//return Matrix which includes only filtered variants (-> "PASS" or "freq-tum" in filter column)
function filter_vcf_file($vcf_file)
{
	$vcf_filtered = new Matrix();
	$vcf_filtered->setHeaders($vcf_file->getHeaders());
	$vcf_filtered->setComments($vcf_file->getComments());
	$i_filter = $vcf_file->getColumnIndex("FILTER");
	
	for($i=0;$i<$vcf_file->rows();$i++)
	{
		if($vcf_file->get($i,$i_filter) == "PASS" || $vcf_file->get($i,$i_filter) == "freq-tum")
		{
			$filtered_row = $vcf_file->getRow($i);
			$vcf_filtered->addRow($filtered_row);
		}
	}
	return $vcf_filtered;
}

//adds a # in the first line of a file
function addCommentCharInHeader($filename)
{
	$file = fopen($filename,"r+");
	$old_contents = file_get_contents($filename);
	fwrite($file,"#");
	fwrite($file,$old_contents);
}

//Main
//get sample ID
$sample_id = "";
if(isset($mutations))
{
	if(strpos($mutations,".gz") !== false) 
	{	
		$sample_id = basename($mutations,"_var_annotated.vcf.gz");
	}else{
		$sample_id = basename($mutations,"_var_annotated.vcf");
	}
}elseif(isset($cnas)){
	$sample_id = basename($cnas,"_cnvs.tsv");
}
//TODO: sample_id if only translocations are given
if(strpos($mutations,".gz") !== false)
{
	exec2("gzip -d -k -f  $mutations");
	$mutations = str_replace(".gz","",$mutations);
}

//only send filtered SNVs to CGI (->those with filter "PASS" entry)
$temp_mutation_file = tempnam(sys_get_temp_dir(),"temp_");
if(isset($mutations))
{
	$vcf_file = Matrix::fromTSV($mutations);
	$vcf_file_filtered = filter_vcf_file($vcf_file);
	
	
	$vcf_file_filtered->toTSV($temp_mutation_file);
	chmod($temp_mutation_file,0666);
	
	//check whether filtered vcf file contains data
	if($vcf_file_filtered->rows() < 1)
	{
		$temp_mutation_file = "";
	}
	
} else {
	//if no mutations given, make $temp_mutation_file empty
	$temp_mutation_file =  "";
}

//convert CNVs to CGI format
$temp_cnv_file = tempnam(sys_get_temp_dir(),"temp_");
chmod($temp_cnv_file,0666);

if(isset($cnas))
{
	transform_cnv_annotations($cnas,$temp_cnv_file);
} else {
	$temp_cnv_file = "";
}
//check whether temp_cnv_file contains data
if($temp_cnv_file != "" && count(file($temp_cnv_file)) <= 1)
{
	$temp_cnv_file = "";
}

//Send data to CGI, save job_id 
$jobId = sendData("Sample: $sample_id",$cancertype,$temp_mutation_file,$temp_cnv_file);
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

//save results
$tmp_file = $parser->tempFile(".zip","temp");
downloadJobResults($tmp_file,$jobId,$user,$token,$url);

exec2("unzip -n $tmp_file -d $destination_folder");


//change file names
if(file_exists($destination_folder."/"."mutation_analysis.tsv"))
{	
	addCommentCharInHeader($destination_folder."/"."mutation_analysis.tsv");
	$parameters = $destination_folder."/"."mutation_analysis.tsv $destination_folder"."/".$sample_id."_cgi_mutation_analysis.tsv";
	$command = $parser->exec("mv",$parameters,true); 
}
if(file_exists($destination_folder."/"."cna_analysis.tsv"))
{
	addCommentCharInHeader($destination_folder."/"."cna_analysis.tsv");
	$parameters = $destination_folder."/"."cna_analysis.tsv $destination_folder/".$sample_id."_cgi_cnv_analysis.tsv";
	$parser->exec("mv",$parameters,true);
}
if(file_exists($destination_folder."/"."drug_prescription.tsv"))
{
	addCommentCharInHeader($destination_folder."/"."drug_prescription.tsv");
	$parameters = $destination_folder."/"."drug_prescription.tsv ".$destination_folder."/".$sample_id."_cgi_drug_prescription.tsv";
	$parser->exec("mv",$parameters,true);
}
if(file_exists($destination_folder."/"."drug_prescription_bioactivities.tsv"))
{
	addCommentCharInHeader($destination_folder."/"."drug_prescription_bioactivities.tsv");
	$parameters = $destination_folder."/drug_prescription_bioactivities.tsv ".$destination_folder."/".$sample_id."_cgi_drug_prescription_bioactivities.tsv";
	$parser->exec("mv",$parameters,true);
}

$parser->exec("rm",$destination_folder."/"."input0*.tsv",true);

if(!$no_del)
{
	delJob($jobId,$user,$token,$url);
}
?>
