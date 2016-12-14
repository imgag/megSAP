<?php

require_once("functions.php");

/**
	@brief Access to the sequence ontology.
	
	@ingroup base
*/
class Miso
{
	private $terms = array();
	private $typedefs = array();
	
	/// Constructor.
	function __construct()
	{
	}
	
	function addTerm($term)
	{
		$this->terms[] = $term;
	}

	function getTermIdByName($term_name, $error = TRUE)
	{
		$terms = $this->terms;
		for($i=0;$i<count($terms);++$i)
		{
			$term = $terms[$i];
			if($term_name==$term->getName())
			{
				return $term->getID();
			}
		}
		
		if($error)	trigger_error("Term '$term_name' was not found.", E_USER_ERROR);
		return null;
	}
	
	
	function getTermNameByID($term_id, $error = TRUE)
	{
		$terms = $this->terms;
		for($i=0;$i<count($terms);++$i)
		{
			$term = $terms[$i];
			if($term_id==$term->getID())
			{
				return $term->getName();
			}
		}
		
		if($error)	trigger_error("Term '$term_id' was not found.", E_USER_ERROR);
		return null;
	}
	
	function getTermChildrenNames($term_id, $recursive = false)
	{
		$children_names = array();
		
		$terms = $this->terms;
		for($i=0;$i<count($terms);++$i)
		{
			$term = $terms[$i];
			if($term->childOf($term_id))
			{
				$children_names[] = $term->getName();
				if($recursive)	$children_names = array_merge($children_names, $this->getTermChildrenNames($term->getID(),true));
			}
		}
		
		return $children_names;		
	}
	
	function getTermChildrenIDs($term_id, $recursive = false)
	{
		$children_ids = array();
		
		$terms = $this->terms;
		for($i=0;$i<count($terms);++$i)
		{
			$term = $terms[$i];
			if($term->childOf($term_id))
			{
				$children_ids[] = $term->getID();
				if($recursive)	$children_ids = array_merge($children_ids, $this->getTermChildrenIDs($term->getID(),true));
			}
		}
		
		return $children_ids;		
	}

	public static function getMisoOBO($thow_if_fails = true)
	{
		$file = repository_basedir()."/data/dbs/Ontologies/so-xp_2_5_3_v2.obo";
		if (file_exists($file))
		{
			return $file;
		}
		else if ($thow_if_fails)
		{
			trigger_error("Miso OBO file '$file' does not exist!", E_USER_ERROR);
		}
		else
		{
			return "";
		}
	}

	public static function isValidTermByID($term_id)
	{
		$miso = self::fromOBO();
		
		$terms = $miso->terms;
		for($i=0;$i<count($terms);++$i)
		{
			$term = $terms[$i];
			if($term->getID()==$term_id && !$term->getObsolete())	return true;
			else if($term->getID()==$term_id && $term->getObsolete())	return false;
		}
		
		return false;
	}
	
	public static function isValidTermByName($term_name)
	{
		$miso = self::fromOBO();
		
		$term_id = $miso->getTermIdByName($term_name,FALSE);
		
		$terms = $miso->terms;
		for($i=0;$i<count($terms);++$i)
		{
			$term = $terms[$i];
			if($term->getID()==$term_id && !$term->getObsolete())	return true;
			else if($term->getID()==$term_id && $term->getObsolete())	return false;
		}
		
		return false;
	}
	
	public static function fromOBO()
	{
		$miso = new Miso();

		$handle = gzopen(Miso::getMisoOBO() , "r");
		$count = 0;
		while (!feof($handle)) 
		{
			$line = nl_trim(fgets($handle));
			
			if($line=="[Term]")
			{
				$term = new Term();

				while(!empty($line))
				{
					$line = nl_trim(fgets($handle));			
					$field_id = trim(substr($line,0,strpos($line,":")));
					$field_value = trim(substr($line,strpos($line,":")+1));
					
					switch($field_id)
					{
						case "id":
							$term->setID($field_value);
							break;
						case "name":
							$term->setName($field_value);
							break;
						case "def":
							$term->setDef($field_value);
							break;
						case "is_a":
							$term->addIsA($field_value);
							break;
						case "is_obsolete":
							$term->setObsolete($field_value);
							break;
					}
				}
				$miso->addTerm($term);
			}
		}
		gzclose($handle);
		
		return $miso;
	}	
}	

class Term
{
	private $id = "";
	private $name = "";
	private $def = "";
	private $is_as = array();
	private $is_obsolete = false;
	
	public function __construct()
	{
	}
	
	function setID($id)
	{
		$this->id = $id;
	}
	
	function getID()
	{
		return $this->id;
	}
	
	function setName($name)
	{
		$this->name = $name;
	}
	
	function getName()
	{
		return $this->name;
	}
	
	function setObsolete($obsolete)
	{
		$this->is_obsolete = $obsolete;
	}
	
	function getObsolete()
	{
		return $this->is_obsolete;
	}
	
	function setDef($def)
	{
		$this->def = $def;
	}
	
	function getIsAs()
	{
		return $this->is_as;
	}
	
	function addIsA($is_a)
	{
		$this->is_as[] = $is_a;
	}
	
	function childOf($term_id)
	{
		for($i=0;$i<count($this->is_as);++$i)
		{
			$is_a = $this->is_as[$i];
			if(strpos($is_a,$term_id)!==false)	return true;
		}
		
		return false;
	}
}	

?>
