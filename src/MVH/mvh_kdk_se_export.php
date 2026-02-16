<?php
/** 
	@page mvh_grz_export
*/

require_once("mvh_functions.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

function json_metadata($cm_data, $tan_k, $rc_data_json, $se_data, $se_data_rep)
{
	global $db_mvh;
	global $sub_id;
	global $parser;
	global $submission_type;
	
	//Teilnahmeerklärung
	$te_present_date = xml_str($cm_data->mvconsentpresenteddate);
	$te_date = xml_str($cm_data->datum_teilnahme);
	$output = [
			"type" => $submission_type,
			"transferTAN" => $tan_k,
			"modelProjectConsent" => [
				"version" => "",
				"date" => $te_present_date,
				"version" => xml_str($cm_data->version_teilnahme),
				"provisions" => [
					[
						"date" => $te_date,
						"purpose" => "sequencing",
						"type" => xml_str($cm_data->particip_4)=="Ja" ? "permit" : "deny"
					],
					[
						"date" => $te_date,
						"purpose" => "case-identification",
						"type" => xml_str($cm_data->particip_4_1)=="Ja" ? "permit" : "deny"
					],
					[
						"date" => $te_date,
						"purpose" => "reidentification",
						"type" => xml_str($cm_data->particip_4_2)=="Ja" ? "permit" : "deny"
					]
				],
			],
		];
	
	//Forschungseinwilligung
	if (xml_str($cm_data->bc_signed)=="Ja")
	{
		//get BC meta data from SE RedCap
		list($bc_type, $bc_item) = get_bc_data_se($se_data_rep);
		
		$active_rcs = [];
		if (starts_with($bc_type, "Kinder")) //generate consent JSON from SE RedCap
		{
			$active_rcs[] = convert_se_kids_bc_to_fhir($bc_item, $se_data, $parser);
		}
		else //adult > search for consent data from SAP
		{
			if (isset($rc_data_json["entry"]))
			{
				foreach($rc_data_json["entry"] as $entry)
				{
					$entry = $entry['resource'];

					//skip not active
					if ($entry['status']!='active') continue;
					
					//skip if V9
					if ($entry['identifier'][0]['system']=="source.ish.document.v09") continue;
					
					$active_rcs[] = $entry;
				}
			}
		}
		if (count($active_rcs)>1) trigger_error("Several active research consents found!", E_USER_ERROR);
		if (count($active_rcs)==0) trigger_error("No active consent resources found, but CM RedCap says that the consent was signed!", E_USER_ERROR);
		
		if (count($active_rcs)>0)
		{
			$output["researchConsents"] = [ $active_rcs[0] ];
		}
	}
	else
	{
		$output["reasonResearchConsentMissing"] = convert_bc_missing(xml_str($cm_data->bc_reason_missing));
	}
			
	return $output;
}

function json_patient($cm_data, $se_data)
{
	global $patient_id;
	
	$output = [
			"id" => $patient_id,
			"gender" => [ "code" => convert_gender(xml_str($se_data->gender))],
			"birthDate" => xml_str($se_data->birthdate),
			"healthInsurance" => [
				"type" => [
					"code"=> convert_coverage($cm_data->coveragetype)
					]
				]
			];
	
	//optional
	$ags = get_raw_value($se_data->psn, "ags");
	if ($ags!="")
	{
		$output["address"] = [ "municipalityCode" => $ags ];
	}
			
	return $output;
}

function json_patient_ref()
{
	global $patient_id;
	return [
			"id" => $patient_id,
			"type" => "Patient"
			];
}

function json_episode_of_care($se_data)
{
	$output = [
			"id"=>"ID_EOC_1",
			"patient" => json_patient_ref(),
			"period" => [
				"start" => xml_str($se_data->datum_kontakt_zse)
				]
			];
	
	return $output;
}

function json_diagnoses($se_data, $se_data_rep)
{
	global $no_seq;
	
	//prepare codes
	$codes = [];
	$icd10 = xml_str($se_data->diag_icd10);
	$icd10_ver = xml_str($se_data->diag_icd10_ver);
	if ($icd10!="")
	{
		$code = get_raw_value($se_data->psn, "diag_icd10");
		if (ends_with($code, '+')) $code = substr($code, 0, -1);
		
		$codes[] = [
			"code" => $code,
			"display" => $icd10,
			"system" => "http://fhir.de/CodeSystem/bfarm/icd-10-gm",
			"version" => $icd10_ver
			
			
		];
	}
	$orpha = xml_str($se_data->diag_orphacode);
	$orpha_ver = xml_str($se_data->diag_orphacode_ver);
	if ($orpha!="")
	{
		$code = get_raw_value($se_data->psn, "diag_orphacode");
		$code = strtr($code, ["ORDO:Orphanet_"=>"ORPHA:"]);
		$codes[] = [
			"code" => $code,
			"display" => "ORPHA:".$orpha,
			"system" => "https://www.orpha.net",
			"version" => $orpha_ver
		];
	}
	if (xml_bool($se_data->orphacode_undiagnostiziert___1, false, "orphacode_undiagnostiziert"))
	{
		$codes[] = [
			"code" => "ORPHA:616874",
			"display" => "Fully investigated rare disorder without a determined diagnosis",
			"system" => "https://www.orpha.net",
		];
	}
	
	$alpha_id_se = xml_str($se_data->diag_se_code);
	if ($alpha_id_se!="")
	{
		$codes[] = [
			"code" => get_raw_value($se_data->psn, "diag_se_code"),
			"display" => $alpha_id_se,
			"system" => "https://www.bfarm.de/DE/Kodiersysteme/Terminologien/Alpha-ID-SE",
		];
	}
	
	//determine onset date from HPO terms
	$onset_date = [];
	foreach($se_data_rep->item as $item)
	{		
		$onset = substr(xml_str($item->beginn_symptome), 0, 7);
		if ($onset=="" || $onset=="0000-01") continue;

		$onset_date[] = $onset;
	}
	if (count($onset_date)==0)
	{
		$onset_date = "";
	}
	else
	{
		asort($onset_date);
		$onset_date = $onset_date[0];
	}
	
	$output = [
			"id"=>"ID_DIAG_1",
			"patient" => json_patient_ref(),
			"recordedOn" => xml_str($se_data->datum_fallkonferenz),
			"verificationStatus" => [
				"code" => ($no_seq ? "unconfirmed" : convert_diag_status(xml_str($se_data->bewertung_gen_diagnostik))),
				],
			"codes" => $codes,
			//missing fields: notes
			];
	if ($onset_date!="") $output["onsetDate"] = $onset_date;
	if (!$no_seq)
	{
		$output["familyControlLevel"] = [
				"code" => convert_diag_recommendation(xml_str($se_data->diagnostik_empfehlung)),
				];
	}
	
	if (count($codes)!=3)
	{
		$output["missingCodeReason"] = [
			"code" => "no-matching-code",
			"display" => "Kein geeigneter Code (ICD-10-GM, ORDO, Alpha-ID-SE) verfügbar",
			"system" => "dnpm-dip/rd/diagnosis/missing-code-reason"
			];
	}
	
	return $output;
}

function json_hpo_history($instance, $se_data_rep)
{
	$output = [];
	
	//search for updates in repeat instances
	foreach($se_data_rep->item as $item)
	{
		if (xml_str($item->hpo_wiedervorst)!=$instance) continue;
		
		$change = xml_str($item->hpo_wiedervorst_aend);

		$element = [
			"date" => xml_str($item->datum_wiedervorst),
			"status" => [
					"code" => convert_hpo_change($change),
				]
			];
		
		$output[] = $element;
	}
	
	return $output;
}

function json_hpos($se_data, $se_data_rep)
{
	global $db_ngsd;
	$output = [];
	
	$num = 1;
	foreach($se_data_rep->item as $item)
	{
		$hpo = xml_str($item->hpo);
		if ($hpo=="") continue;
		
		//convert HPO name to id
		$hpo_id = $db_ngsd->getValue("SELECT hpo_id FROM hpo_term WHERE name LIKE '{$hpo}'", "");
		if ($hpo_id=="") //fallback to ID (RedCap ontology handling bug that causes terms to have ID instead of name)
		{
			$hpo_id = $db_ngsd->getValue("SELECT hpo_id FROM hpo_term WHERE hpo_id LIKE '{$hpo}'", "");
		}
		if ($hpo_id=="")
		{
			trigger_error("Could not convert HPO name '{$hpo}' to ID using the NGSD hpo_term table!", E_USER_WARNING);
			continue;
		}
		
		//create entry
		$entry = [
		"id" => "ID_HPO_{$num}",
		"patient" => json_patient_ref(),
		"recordedOn" => xml_str($se_data->datum_fallkonferenz),
		"value" => [
			"code" => $hpo_id,
			"display" => $hpo,
			"system" => "https://hpo.jax.org",
			]
		];
		
		//add optional stuff to entry
		$onset = substr(xml_str($item->beginn_symptome), 0, 7);
		if ($onset!="") $entry["onsetDate"] = $onset;
		
		$hpo_ver = xml_str($item->version_hpo);
		if ($hpo_ver!="") $entry["value"]["version"] = $hpo_ver;
		
		$history = json_hpo_history(xml_str($item->redcap_repeat_instance), $se_data_rep);
		if (count($history)>0) $entry["status"] = [ "history" => $history ];
		
		$output[] = $entry;
		++$num;
	}
	
	return $output;
}

function json_gmfcs($se_data)
{
	$output = [];
	$num = 1;
	
	$gmf = xml_str($se_data->diag_gmfcs);
	if ($gmf!="")
	{
		$output[] = [
			"id" => "ID_GMFCS_{$num}",
			"patient" => json_patient_ref(),
			"effectiveDate" => xml_str($se_data->datum_fallkonferenz),
			"value" => [
				"code" => $gmf,
				"display" => "Level ".$gmf,
				"system" => "Gross-Motor-Function-Classification-System",			
				]
			];
	}
	
	//TODO no example with more than one GMFCS entry, i.e. changed GMFCS during follow-up
	
	return $output;
}

function json_followups($se_data_rep)
{
	$output = [];
	
	foreach($se_data_rep->item as $item)
	{
		$date = xml_str($item->datum_wiedervorst);
		if ($date=="") continue;
		
		$output[] = [
			"date" => $date,
			"patient" => json_patient_ref()
			];
	}
	
	return $output;
}


function json_hospitalization($se_data)
{
	$output = [
		 "numberOfStays" => [
			"code" => convert_hospitalization_stays($se_data->anzahl_stat_behandlungen_en)
			],
		 "numberOfDays" => [
			"code" => convert_hospitalization_days($se_data->dauer_stat_vortherapie_en)
			],
		];
	
	return $output;
}

function json_care_plan_1($se_data) //This object models the decisions of the first "Fallkonferenz". i.e. if the patient is included into the Modellvorhaben and sequencing is performed
{
	global $no_seq;
	
	$output = [
			"id"=>"ID_CARE_PLAN_1",
			"patient" => json_patient_ref(),
			"issuedOn" => xml_str($se_data->datum_fallkonferenz),
			];
			
	if ($no_seq)
	{
		$output["noSequencingPerformedReason"] = [ "code" => convert_noseq_reason($se_data->fallkonferenz_grund)];
	}
	
	return $output;	
}

function json_supporting_variants($variant_repeat_instances, $se_data_rep)
{
	global $var2id;
	$output = [];
	
	foreach($variant_repeat_instances as $instance)
	{
		if ($instance=="") continue;
		
		//search for variant in repeat instances
		foreach($se_data_rep->item as $item)
		{
			if (xml_str($item->redcap_repeat_instrument)!="Varianten") continue;
			if (xml_str($item->redcap_repeat_instance)!=$instance) continue;
			
			$var = xml_str($item->variante);
			$id = $var2id[$var];
			$output[] = [
				"display" => $var,
				"variant" =>
					[
						"id"=> $id,
						"display" => $var
					]
				];
			
			//missing: gene (already listed in variant section)
		}
	}
	
	return $output;
}

function json_care_plan_2($se_data, $se_data_rep) //'carePlan' is misleading. This object captures the decisions of the second "Fallkonferenz".
{	
	//prepare therapy recommendations
	$therapy_recoms = [];
	$num = 1;
	foreach($se_data_rep->item as $item)
	{
		$therapy_type = xml_str($item->therapie_art);
		if ($therapy_type=="") continue;
		
		//create entry
		$entry = [
		"id" => "ID_THERAPY_{$num}",
		"patient" => json_patient_ref(),
		"issuedOn" => xml_str($se_data->klin_datum_fallkonferenz),
		"category" => [
			"code" => convert_therapy_category($therapy_type),
			],
		"type" => [
			"code" => convert_therapy_type($item->therapie_beschreibung),
			]
		];
		
		//supporting variants (optional)
		$variant_repeat_instances = [
			xml_str($item->variante1),
			xml_str($item->variante2),
			xml_str($item->variante3),
			xml_str($item->variante4),
			xml_str($item->variante5),
			];
		$support = json_supporting_variants($variant_repeat_instances, $se_data_rep);
		if (count($support)>0) $entry['supportingVariants'] = $support;
		
		//missing: medication (not in SE data)

		$therapy_recoms[] = $entry;
		++$num;
	}
	
	//prepare study enrollment recommendations
	$study_recoms = [];
	$num = 1;
	foreach($se_data_rep->item as $item)
	{
		$study_register = xml_str($item->studien_register);
		if ($study_register=="") continue;
		
		//create entry
		$entry = [
		"id" => "ID_STUDY_{$num}",
		"patient" => json_patient_ref(),
		"issuedOn" => xml_str($se_data->klin_datum_fallkonferenz),
		"study" => [
				[
					"id" => xml_str($item->studien_id),
					"system" => convert_study_register($study_register),
					"type" => "Study",
					"display" => xml_str($item->studienname)
				]
			]
		];
		
		//supporting variants (optional)
		$variant_repeat_instances = [
			xml_str($item->studie_variante1),
			xml_str($item->studie_variante2),
			xml_str($item->studie_variante3),
			xml_str($item->studie_variante4),
			xml_str($item->studie_variante5),
			];
		$support = json_supporting_variants($variant_repeat_instances, $se_data_rep);
		if (count($support)>0) $entry['supportingVariants'] = $support;

		$study_recoms[] = $entry;
		++$num;
	}
	
	//prepare clinical management recommendations
	$clin_recoms = [];
	$num = 1;
	foreach($se_data_rep->item as $item)
	{
		$clinical_mgmt = xml_str($item->klinisches_management_dr);
		if ($clinical_mgmt=="") continue;
		
		//create entry
		$entry = [
		"id" => "ID_THERAPY_{$num}",
		"patient" => json_patient_ref(),
		"issuedOn" => xml_str($se_data->klin_datum_fallkonferenz),
		"type" => [
			"code" => convert_clinincal_management($clinical_mgmt),
			]
		];

		$clin_recoms[] = $entry;
		++$num;
	}
	if (count($clin_recoms)>1) trigger_error("More than one clinical management recommendation found. Only one is allowed in KDK-SE. Using first clinical management recommendation!", E_USER_WARNING);
	
	$output = [
			"id"=>"ID_CARE_PLAN_2",
			"patient" => json_patient_ref(),
			"issuedOn" => xml_str($se_data->klin_datum_fallkonferenz),
			"therapyRecommendations" => $therapy_recoms,
			"studyEnrollmentRecommendations" => $study_recoms,
			];
	
	//optional stuff
	$recom_counceling = xml_bool($se_data->empf_hg_beratung, true, "empf_hg_beratung");
	if (!is_null($recom_counceling))
	{
		$output["geneticCounselingRecommended"] = $recom_counceling;
	}
	
	$recom_reeval = xml_bool($se_data->empf_reeval, true, "empf_reeval");
	if (!is_null($recom_reeval))
	{
		$output["reevaluationRecommended"] = $recom_reeval;
	}
	
	if (count($clin_recoms)>=1)
	{
		$output["clinicalManagementRecommendation"] = $clin_recoms[0];
	}
	
	//missing: notes
	
	return $output;	
}

function json_ngs_report($cm_data, $se_data, $info, $results)
{
	global $db_ngsd;
	
	$output = [
			"id" => "ID_NGS_REPORT_1",
			"patient" => json_patient_ref(),
			"issuedOn" => xml_str($cm_data->gen_finding_date),
			"type" => [
				"code" => ($info['sys_type']=="lrGS" ? "genome-long-read" : "genome-short-read"),
				],
			"sequencingInfo" =>[
				"platform" => [
					"code" => convert_sequencing_platform($info["device_type"])
					],
				"kit" => $info['sys_name']
				]
		];
		
	//outcome (only for unsolved - if a case is solved, is documented in "diagnoses>verificationStatus")
	$outcome = convert_outcome($db_ngsd->getValue("SELECT outcome FROM diag_status WHERE processed_sample_id='".$info['ps_id']."'"));
	if ($outcome!="")
	{
		$output["conclusion"] = ["code"=>$outcome];
	}
	
	//results, i.e. causal variants
	if (count($results)>0)
	{
		$output['results'] = $results;
	}
	
	return $output;
}

//returns a map from gene symbol to HGNC id. Genes need to be calculated because the variant text shows "(many)" if more than 10 genes are affected by the variant.
function genes_overlapping($chr, $start, $end)
{
	global $db_ngsd;
	global $parser;
	
	//determine overlapping genes
	$tmp_file = $parser->tempFile(".bed");
	file_put_contents($tmp_file, "{$chr}\t".($start-1)."\t{$end}");
	list($stdout) = $parser->execApptainer("ngs-bits", "BedAnnotateGenes", "-in {$tmp_file}");
	$tmp = explode("\t", nl_trim($stdout[0]))[3];
	
	//no gene > extend by 5000 bases
	if ($tmp=="")
	{
		file_put_contents($tmp_file, "{$chr}\t".($start-5000)."\t".($end+5000));
		list($stdout) = $parser->execApptainer("ngs-bits", "BedAnnotateGenes", "-in {$tmp_file}");
		$tmp = explode("\t", nl_trim($stdout[0]))[3];
	}
	
	$genes = [];
	foreach(explode(",", $tmp) as $gene)
	{
		$gene = trim($gene);
		if ($gene=="") continue;
		
		//determine HGNC ID
		$genes[$gene] = "HGNC:".$db_ngsd->getValue("SELECT hgnc_id FROM gene WHERE symbol='{$gene}'");;
	}
	
	return $genes;
}

function json_results($se_data, $se_data_rep, $info)
{
	global $var2id; //we need variant IDs to link to variants, so we store the association in a global variable
	global $db_ngsd;
	global $parser;
	
	$results = [];
	
	foreach($se_data_rep->item as $item)
	{
		$var = xml_str($item->variante);
		if ($var=="") continue;
		
		//init
		$ps_id = $info['ps_id'];
		$rc_id = $db_ngsd->getValue("SELECT id FROM report_configuration WHERE processed_sample_id='{$ps_id}'");
		
		//determine ID
		$id = "ID_VAR_".(count($var2id)+1);
		$var2id[$var] = $id;
		
		//initialize variant datastructure
		$var_data = [];
		$var_data['id'] = $id;
		$var_data['patient'] = json_patient_ref();
		
		//add significance
		$var_data['significance'] = ['code' => convert_var_type(xml_str($item->typ_variante))];
		
		//small variants - example: chr8:10623069-10623069 G>A (RP1L1)
		if (starts_with($var, "chr"))
		{
			//convert GSvar format to VCF format
			for($i=0; $i<strpos($var, ' '); ++$i)
			{
				if ($var[$i]=='-') $var[$i] = ' ';
			}
			$var = strtr($var, ":>", "  ");
			list($chr, $start, $end, $ref, $alt) = explode(" ", $var);
			list($chr2, $start2, $ref2, $alt2) = gsvar_to_vcf("GRCh38", $chr, $start, $ref, $alt);
			$var_data['chromosome'] = $chr2;
			$var_data['startPosition'] = (int)$start2;
			$var_data['endPosition'] = $start2 + strlen($ref)-1; //not really defined how to calcualte it...
			$var_data['ref'] = $ref2;
			$var_data['alt'] = $alt2;
			
			//zygosity
			$var_id = get_variant_id($db_ngsd, $chr, $start, $end, $ref, $alt);
			$genotype = $db_ngsd->getValue("SELECT genotype FROM detected_variant WHERE processed_sample_id='{$ps_id}' AND variant_id='{$var_id}'");
			$var_data['zygosity'] = ['code' => convert_genotype($genotype)];
			
			//genes
			$genes = genes_overlapping($chr, $start, $end);
			foreach($genes as $gene => $hgnc_id)
			{
				if (!isset($var_data['genes'])) $var_data['genes'] = [];
				$var_data['genes'][] = ['code'=>$hgnc_id, 'disply'=>$gene];
			}
			
			//acmgClass
			$class = $db_ngsd->getValue("SELECT class FROM variant_classification WHERE variant_id='{$var_id}'", "");
			if (is_numeric($class))
			{
				$var_data['acmgClass'] = ['code' => $class];
			}
			
			//modeOfInheritance
			$inheritance = $db_ngsd->getValue("SELECT inheritance FROM report_configuration_variant WHERE report_configuration_id='{$rc_id}' AND variant_id='{$var_id}'");
			$var_data['inheritance'] = ['code' => convert_inheritance($inheritance)];
			
			//cDNAChange, proteinChange
			$tmp_file = $parser->tempFile(".bed");
			$transcripts = []; //best transcripts of overlapping 
			foreach($genes as $gene => $hgnc_id)
			{
				file_put_contents($tmp_file, $gene);
				list($stdout) = $parser->execApptainer("ngs-bits", "GenesToTranscripts", "-in {$tmp_file} -mode best");
				foreach($stdout as $line)
				{
					$line = nl_trim($line);
					if (!starts_with($line, $gene)) continue;
					
					$transcript = trim(explode("\t", $line)[1]);
					if ($transcript=="") continue;

					$transcripts[] = $transcript;
				}
			}
			$consequences = [];
			$parts = explode(",", $db_ngsd->getValue("SELECT coding FROM variant WHERE id='{$var_id}'"));
			foreach($parts as $part)
			{
				$hit = false;
				foreach($transcripts as $transcript)
				{
					if (contains($part, $transcript.'.')) $hit = true;
				}
				if ($hit)
				{
					$consequences[] = explode(":", $part);
				}
			}
			if (count($consequences)>0)
			{
				$parts = $consequences[0];
				if (count($parts)==7)
				{
					list($t_gene, $t_trans, $t_type, $t_impact, $t_exon, $t_cdna, $t_prot) = $consequences[0];
					if ($t_cdna!="")
					{
						$var_data['cDNAChange'] = $t_trans.":".$t_cdna;
					}
					if($t_prot!="")
					{
						$var_data['proteinChange'] = $t_trans.":".$t_prot;
					}
				}
				else trigger_error("Could not determine cDNA and protein change of variant (not in NGSD)", E_USER_NOTICE);
				
				if (count($consequences)>1) trigger_error("More than one variant consequence found. Only one is allowed in KDK-SE. Using first variant consequence!", E_USER_WARNING);
			}
			
			//add
			if (isset($results['smallVariants'])) $results['smallVariants'] = [];
			$results['smallVariants'][] = $var_data;
			
			//missing: publications (not in NGSD), externalIds (which?), segregationAnalysis (not in NGSD), acmgCriteria (not in NGSD), localization (depends on gene/transcript), gDNAChange (duplicate of chr, start, ref, alt)
		}
		//CNV - example: CNV: chr6:64839979-64993979 CN=0 size=154.000kb (EYS)
		else if (starts_with($var, "CNV:"))
		{
			//parse variant text
			$var = strtr(substr($var, 5), ":-", "  ");
			list($chr, $start, $end, $cn) = explode(" ", $var);
			$var_data['chromosome'] = $chr;
			$var_data['startPosition'] = (int)$start;
			$var_data['endPosition'] = (int)$end;
			$var_data['type'] = ['code' => convert_cn_to_type($cn, $chr, convert_gender(xml_str($se_data->gender)))];
			
			//genes
			$genes = genes_overlapping($chr, $start, $end);
			foreach($genes as $gene => $hgnc_id)
			{
				if (!isset($var_data['genes'])) $var_data['genes'] = [];
				$var_data['genes'][] = ['code'=>$hgnc_id, 'disply'=>$gene];
			}
			
			//acmgClass
			$cs_id = $db_ngsd->getValue("SELECT id FROM cnv_callset WHERE processed_sample_id='{$ps_id}'");
			$var_id = $db_ngsd->getValue("SELECT id FROM cnv WHERE cnv_callset_id='{$cs_id}' AND chr='{$chr}' AND start='{$start}' AND end='{$end}'");
			$class = $db_ngsd->getValue("SELECT class FROM report_configuration_cnv WHERE cnv_id='{$var_id}'", "");
			if (is_numeric($class))
			{
				$var_data['acmgClass'] = ['code' => $class];
			}
			
			//modeOfInheritance
			$inheritance = $db_ngsd->getValue("SELECT inheritance FROM report_configuration_cnv WHERE cnv_id='{$var_id}'", "");
			$var_data['inheritance'] = ['code' => convert_inheritance($inheritance)];
			
			//add
			if (isset($results['copyNumberVariants'])) $results['copyNumberVariants'] = [];
			$results['copyNumberVariants'][] = $var_data;
			
			//missing: publications (not in NGSD), externalIds (which?), segregationAnalysis (not in NGSD), acmgCriteria (not in NGSD), localization (depends on gene/transcript), gDNAChange (duplicate of chr, start, end) , cDNAChange/proteinChange (how?)
		}
		//SV - example: SV-DUP at chr9:17051726-138024890 genotype=het (many)
		else if (starts_with($var, "SV-"))
		{
			//parse variant text
			
			list(, $type, , $chr, $start, $end) = explode(" ", strtr($var, ":-", "  "));
			if ($type!="DEL" && $type!="DUP" && $type!="INV" && $type!="INS") trigger_error("Unsupported structural variant type '{$type}'!", E_USER_ERROR);
			
			//gDNAChange
			$var_data['gDNAChange'] = chr2NC($chr).":g.{$start}_{$end}".strtolower($type);
			
			//genes
			$genes = genes_overlapping($chr, $start, $end);
			foreach($genes as $gene => $hgnc_id)
			{
				if (!isset($var_data['genes'])) $var_data['genes'] = [];
				//$var_data['genes'][] = ['code'=>$hgnc_id, 'disply'=>$gene];
			}
			
			//acmgClass
			$cs_id = $db_ngsd->getValue("SELECT id FROM sv_callset WHERE processed_sample_id='{$ps_id}'");
			$table = sv_type_to_table($type);
			if ($type=="INS")
			{
				$var_id = $db_ngsd->getValue("SELECT id FROM {$table} WHERE sv_callset_id='{$cs_id}' AND chr='{$chr}' AND pos=>'{$start}' AND pos<='{$end}'");
			}
			else
			{
				$var_id = $db_ngsd->getValue("SELECT id FROM {$table} WHERE sv_callset_id='{$cs_id}' AND chr='{$chr}' AND start_min='{$start}' AND end_max='{$end}'");
			}
			
			$class = $db_ngsd->getValue("SELECT class FROM report_configuration_sv WHERE {$table}_id='{$var_id}'", "");
			if (is_numeric($class))
			{
				$var_data['acmgClass'] = ['code' => $class];
			}
			
			//modeOfInheritance
			$inheritance = $db_ngsd->getValue("SELECT inheritance FROM report_configuration_sv WHERE {$table}_id='{$var_id}'", "");
			$var_data['inheritance'] = ['code' => convert_inheritance($inheritance)];
			
			//zygosity
			$genotype = $db_ngsd->getValue("SELECT genotype FROM {$table} WHERE id='{$var_id}'");
			$var_data['zygosity'] = ['code' => convert_genotype($genotype)];
			
			//add
			if (isset($results['structuralVariants'])) $results['structuralVariants'] = [];
			$results['structuralVariants'][] = $var_data;
			
			//missing: publications (not in NGSD), externalIds (which?), segregationAnalysis (not in NGSD), acmgCriteria (not in NGSD), localization (depends on gene/transcript), cDNAChange/proteinChange (how?), iscnDescription (can be calculated from HGVS.g)
		}
		else if (starts_with($var, "RE:")) //RE
		{
			//missing: REs (not modelled in SE-DIP format)
		}
		else trigger_error("Unhandled variant type for variant '{$var}'!", E_USER_ERROR);	

	}
	
	//missing: autozygosity (there is no definition how to calculate it), REs (not modelled in SE-DIP format)
	
	return $results;
}

//parse command line arguments
$parser = new ToolBase("mvh_kdk_se_export", "KDK-SE export for Modellvorhaben.");
$parser->addInt("cm_id", "ID in case management RedCap database.", false);
$parser->addFlag("clear", "Clear export and QC folder before running this script.");
$parser->addFlag("test", "Test mode.");
extract($parser->parse($argv));

//init
$db_mvh = DB::getInstance("MVH");
$db_ngsd = DB::getInstance("NGSD");
$mvh_folder = get_path("mvh_folder");
$time_start = microtime(true);

//check that case ID is valid
$id = $db_mvh->getValue("SELECT id FROM case_data WHERE cm_id='{$cm_id}'");
if ($id=="") trigger_error("No case with id '{$cm_id}' in MVH database!", E_USER_ERROR);

//get patient identifer (Fallnummer from CM)
$cm_data = get_cm_data($db_mvh, $id);

//start export
print "CM ID: {$cm_id} (MVH DB id: {$id} / CM Fallnummer: ".xml_str($cm_data->case_id).")\n";
print "Export start: ".date('Y-m-d H:i:s')."\n";
$folder = realpath($mvh_folder)."/kdk_se_export/{$cm_id}/";
if ($clear) exec2("rm -rf {$folder}");
exec2("mkdir -p {$folder}/metadata/");
print "export folder: {$folder}\n";

//determine tanG==VNg
$sub_ids = $db_mvh->getValues("SELECT id FROM `submission_kdk_se` WHERE ".($test ? "" : "status='pending' AND")." case_id='{$id}'");
if (count($sub_ids)!=1) trigger_error(count($sub_ids)." pending KDK-SE submissions for case {$cm_id}. Must be one!", E_USER_ERROR);
$sub_id = $sub_ids[0];
print "ID in submission_kdk_se table: {$sub_id}\n";
$tan_k = $db_mvh->getValue("SELECT tank FROM submission_kdk_se WHERE id='{$sub_id}'");
print "TAN: {$tan_k}\n";
$patient_id = $db_mvh->getValue("SELECT pseudok FROM submission_kdk_se WHERE id='{$sub_id}'");
print "patient pseudonym: {$patient_id}\n";


//get data from MVH database
$se_data = get_se_data($db_mvh, $id);
$se_data_rep = get_se_data($db_mvh, $id, true);
$rc_data_json = get_rc_data_json($db_mvh, $id);

//check if that patient was included into the Modellvorhaben, i.e. sequencing is/was performed
$no_seq = !xml_bool($se_data->aufnahme_mvh, false, "aufnahme_mvh");

//get data from NGSD
if ($no_seq)
{
	print "index sample: none - no sequending was performed in MVH for this patient!\n";
}
else
{
	$ps = $db_mvh->getValue("SELECT ps FROM case_data WHERE id='{$id}'");
	print "index sample: {$ps}\n";
	$info = get_processed_sample_info($db_ngsd, $ps);
}

//create base JSON
print "\n";
print "### creating JSON ###\n";

//create results section (contained variants are referenced in therapy and study recommendations, so we need to create them first to have to the IDs later on)
$submission_type = $db_mvh->getValue("SELECT type FROM submission_kdk_se WHERE id='{$sub_id}'");
$var2id = [];
$results = $no_seq ? [] : json_results($se_data, $se_data_rep, $info);
$json = [
	"metadata" => json_metadata($cm_data, $tan_k, $rc_data_json, $se_data, $se_data_rep),
	"patient" => json_patient($cm_data, $se_data),
	"episodesOfCare" => [ json_episode_of_care($se_data) ],
	"diagnoses" => [ json_diagnoses($se_data, $se_data_rep) ],
	"hpoTerms" => json_hpos($se_data, $se_data_rep),
	"hospitalization" => json_hospitalization($se_data),
	"carePlans" => [ json_care_plan_1($se_data) ],
	//TODO: therapies - wait for answer on email "Docu eigentliche Therapie" from 31.08.2025
	];

//add optional parts to JSON
if (!$no_seq)
{
	$json["carePlans"][] = json_care_plan_2($se_data, $se_data_rep);
	$json["ngsReports"] = [ json_ngs_report($cm_data, $se_data, $info, $results)];
}
$gmfcs = json_gmfcs($se_data, $se_data_rep);
if (count($gmfcs)>0)
{
	$json["gmfcsStatus"] = $gmfcs;
}

$followups = json_followups($se_data_rep);
if (count($followups)>0)
{
	$json["followUps"] = $followups;
}

//write JSON
$json_file = "{$folder}/metadata/metadata.json";
file_put_contents($json_file, json_encode($json, JSON_PRETTY_PRINT));

print "creating JSON took ".time_readable(microtime(true)-$time_start)."\n";
$time_start = microtime(true);

//validate JSON
print "\n";
print "### validating JSON ###\n";

$url = "https://".($test ? "preview.dnpm-dip.net" : "dnpm-dip.med.uni-tuebingen.de")."/api/rd/etl/patient-record:validate";
print "URL: {$url}\n";
$ch = curl_init();
curl_setopt($ch, CURLOPT_POST, true);
curl_setopt($ch, CURLOPT_URL, $url);
curl_setopt($ch, CURLOPT_POSTFIELDS, file_get_contents($json_file));
curl_setopt($ch, CURLOPT_RETURNTRANSFER, true);
if (!$test) //production server is in UKT network
{
	curl_setopt($ch, CURLOPT_PROXY, '');
}
$result = curl_exec($ch);
if ($result===false) trigger_error('CURL ERROR: '.curl_error($ch), E_USER_ERROR);
curl_close($ch);

//parse/show validation result
print "exit code: ".curl_getinfo($ch, CURLINFO_HTTP_CODE)."\n";
print "\n";
$validation_error = false;
foreach(json_decode($result)->issues as $issue)
{
	if (isset($issue->message) && contains($issue->message, "Fehlende Angabe 'Krankenkassen-IK'")) continue;
	
	$type = strtoupper($issue->severity);
	if ($type=="ERROR") $validation_error = true;
	print "{$type}: ".(isset($issue->details) ? $issue->details : $issue->path.":".$issue->message)."\n";
	print "\n";
}
if ($validation_error) trigger_error("Validation failed!", E_USER_ERROR);

print "validating JSON took ".time_readable(microtime(true)-$time_start)."\n";
$time_start = microtime(true);

//upload JSON
print "\n";
print "### uploading JSON ###\n";

$url = "https://".($test ? "preview.dnpm-dip.net" : "dnpm-dip.med.uni-tuebingen.de")."/api/rd/etl/patient-record";
print "URL: {$url}\n";
$ch = curl_init();
curl_setopt($ch, CURLOPT_POST, true);
curl_setopt($ch, CURLOPT_URL, $url);
curl_setopt($ch, CURLOPT_POSTFIELDS, file_get_contents($json_file));
curl_setopt($ch, CURLOPT_RETURNTRANSFER, true);
if (!$test) //production server is in UKT network
{
	curl_setopt($ch, CURLOPT_PROXY, '');
}
$result = curl_exec($ch);
if ($result===false) trigger_error('CURL ERROR: '.curl_error($ch), E_USER_ERROR);
curl_close($ch);

//parse/show upload result
$exit_code = curl_getinfo($ch, CURLINFO_HTTP_CODE);
print "exit code: {$exit_code}\n";
print "\n";
foreach(json_decode($result)->issues as $issue)
{
	if (isset($issue->message) && contains($issue->message, "Fehlende Angabe 'Krankenkassen-IK'")) continue;
	
	$type = strtoupper($issue->severity);
	if ($type=="ERROR") $validation_error = true;
	print "{$type}: ".(isset($issue->details) ? $issue->details : $issue->path.":".$issue->message)."\n";
	print "\n";
}
if ($exit_code!="200" && $exit_code!="201") trigger_error("Upload failed!", E_USER_ERROR);

print "uploading JSON took ".time_readable(microtime(true)-$time_start)."\n";
$time_start = microtime(true);

//if upload successfull, add 'Pruefbericht' to CM RedCap
if ($submission_type=='initial' && !$test)
{
	print "Adding Pruefbericht to CM RedCap...\n";
	add_submission_to_redcap($cm_id, "K", $tan_k);
}

//archive metadata JSON
copy($json_file, $mvh_folder."/metadata_archive/KDK_SE/{$cm_id}.json");

//clean up export folder if successfull
if (!$test)
{
	exec2("rm -rf {$folder}");
}

print "cleanup took ".time_readable(microtime(true)-$time_start)."\n";

?>