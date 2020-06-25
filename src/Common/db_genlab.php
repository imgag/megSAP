<?php

require_once("genomics.php");

//Connection to GenLab database (MSSQL)
class GenLabDB
{
	protected static $instance = null;
	protected $connection = null;
	
	public static function isEnabled()
	{
		$db_properties = array('genlab_host', 'genlab_user', 'genlab_pass');
		foreach($db_properties as $db_prop)
		{
			if (trim(get_path($db_prop))=="")
			{
				return false;
			}
		}
		
		return true;
	}
	
	//Returns current instance of the DB connection.
	public static function getInstance()
	{		
		//create a new connection if no connection exists OR the connection is no longer alive
		if(is_null(self::$instance) || !self::$instance->isAlive())
		{
			$host = get_path('genlab_host');
			$user = get_path('genlab_user');
			$pass = get_path('genlab_pass');
			
			self::$instance = new self($host, $user, $pass);
		}
		return self::$instance;
	}
	
	protected function __construct($host, $user, $pass)
	{		
		try 
		{
			$this->connection = new PDO("sqlsrv:Driver=ODBC Driver 17 for SQL Server;Server={$host}", $user, $pass);
			$this->connection->setAttribute(PDO::ATTR_ERRMODE, PDO::ERRMODE_EXCEPTION);	
		}
		catch(PDOException $e)
		{
			trigger_error("GenLabDB constructor encountered error: ".$e->getMessage(), E_USER_ERROR);
		}
	}
	
	public function __destruct()
	{
	}
	
	///quotes special characters for use in a SQL query
	function quote($string)
	{
		return $this->connection->quote($string);
	}
	
	///checks if the connection is still alive (it can be closed by the server if idle too long)
	function isAlive() 
	{
        try 
		{
            $this->connection->query('SELECT 1');
        } 
		catch (PDOException $e) 
		{
			return false;
        }
 
        return true;
    }
	
	public function executeQuery($query)
	{
		try
		{		
			$prepared_query = $this->connection->prepare($query);
			$prepared_query->execute();
			$result = $prepared_query->fetchAll(PDO::FETCH_ASSOC);	
		}
		catch(PDOException $e)
		{
			trigger_error("GenLabDB::executeQuery encountered error: ".$e->getMessage(), E_USER_ERROR);
		}
		
		return $result;
	}
	
	//Returns the 
	public function getValue($query, $default = null)
	{
		$result = $this->executeQuery($query);
			
		if (is_null($default) && count($result)==0) trigger_error("GenLabDB::getValue query returned zero rows: $query", E_USER_ERROR);
		if (is_null($default) && count($result)>1) trigger_error("GenLabDB::getValue query returned more than one row: $query", E_USER_ERROR);
		@$result = $result[0];
		if (is_null($default) && count($result)==0) trigger_error("GenLabDB::getValue query returned zero columns: $query", E_USER_ERROR);
		if (is_null($default) && count($result)>1) trigger_error("GenLabDB::getValue query returned more than one column: $query", E_USER_ERROR);
		@$result = array_values($result);
		@$result = $result[0];
		
		if ($result=="" || is_null($result)) return $default;
		
		return $result;
	}

	///Returns the query result as an array (works only for queries that return one value per row)
	public function getValues($query)
	{
		$output = array();
		
		$result = $this->executeQuery($query);
		foreach($result as $row)
		{
			if (count($row)!=1) trigger_error("GenLabDB::getValues query returned row with count!=1: ".implode("\t", $row), E_USER_ERROR);
			
			foreach($row as $key => $value)
			{
				$output[] = $value;
			}
		}
		
		return $output;
	}
}
