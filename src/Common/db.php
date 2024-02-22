<?php

require_once("genomics.php");

//Connection to NGSD database (MySQL)
class DB
{
	protected static $instances = array();
	protected $name = "";
	protected $connection = null;
	protected $prep_stmts = array();
	protected $stmt_counter = 0;
	protected $log_messages = array();
	const PERSISTENT_STORAGE = TRUE;
	
	///Returns a new instance
	protected static function get_instance($db)
	{
		$host = get_db($db, 'db_host');
		$name = get_db($db, 'db_name');
		$user = get_db($db, 'db_user');
		$pass = get_db($db, 'db_pass');
		
		return new self($host, $name, $user, $pass, $db);
	}
	
	/**
	 * Returns current instance of the DB connection and - if not yet available - creates a new one.
 	 * 
	 * @param $db Name of database to access
	 * @param $persistent If true, a persistent singleton instance is returned. If false, an instance is returned that is deleted when it leaves the scope.
	 */
	public static function getInstance($db, $persistent = true)
	{
		if (!$persistent) return self::get_instance($db);
		
		//create a new connection if no connection exists OR the connection is no longer alive
		if(!isset(self::$instances[$db]) || !self::$instances[$db]->isAlive())
		{
			self::$instances[$db] = self::get_instance($db);
		}
		return self::$instances[$db];
	}
	
	protected function __construct($host, $name, $user, $pass, $ini_name)
	{		
		try 
		{
			$this->connection = new PDO("mysql:host={$host};dbname={$name}", $user, $pass, array(PDO::MYSQL_ATTR_INIT_COMMAND => "SET NAMES utf8"));
			$this->connection->setAttribute(PDO::ATTR_ERRMODE, PDO::ERRMODE_EXCEPTION);	
			$hash = $this->prepare("SET SESSION sql_mode = 'STRICT_TRANS_TABLES';");
			$this->execute($hash);
			
			$this->name = $ini_name;
		}
		catch(PDOException $e)
		{
			$this->error($hash, $e->getMessage());
		}
	}
	
	public function __destruct()
	{
	}
	
	function name()
	{
		return $this->name;
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
	
	/**
	 * Prepares a statement for execution.
	 * 
	 * @param string $stmt SQL query statement
	 * @return string Hash for identification of the current prepared query statement
	 */
	public function prepare($stmt)
	{
		$hash = md5($stmt);
		$this->log_messages[$hash] = array('query' => "", 'bind' => array());
		
		try
		{
			$this->log_messages[$hash]['query'] = $stmt;
			$this->prep_stmts[$hash] = $this->connection->prepare($stmt);			
		}
		catch(PDOException $e)
		{
			$this->error($hash, $e->getMessage());			
		}
		
		return $hash;
	}
	
	/**
	 * Executes a prepared SQL query.
	 * 
	 * @param string $hash Hash to identify prepared statement
	 * @param boolean $persist Whether query should be kept or deleted after execution
	 * @return array Results of the query are returned 
	 */
	public function execute($hash, $persist = FALSE)
	{
		try
		{
			$result = $this->prep_stmts[$hash]->execute();
			
			if(!$persist)
			{
				$this->unsetStmt($hash);
			}
		}
		catch(PDOException $e)
		{
			$this->error($hash, $e->getMessage());
		}
		
		return $result;
	}

	/**
	 *Executes a query without binding of parameters and returns nothing
	 *  
	 * @param string $sql SQL Query to execute
	 * @param array $par associative array with parameter bindings
	 * @param boolean $log Should logging be activated for this query
	 * @return int The number of affected rows
	 */
	public function executeStmt($sql, $par = NULL)
	{
		//prepare
		$hash = $this->prepare($sql);
		
		//bind
		if(!empty($par))
		{
			foreach($par as $key => $value)
			{
				$this->bind($hash, $key, $value);
			}			
		}
		
		//execute
		try
		{
			$this->prep_stmts[$hash]->execute();
			$affected_rows = $this->prep_stmts[$hash]->rowCount();
		}
		catch(PDOException $e)
		{
			$this->error($hash, $e->getMessage());
		}
		
		//delete
		$this->unsetStmt($hash);
		
		return $affected_rows;
	}
	
	/**
	 * Executes a query without binding of parameters and returns the result
	 * 
	 * @param string $query SQL query to execute
	 * @param boolean $log Sets logging
	 * @return array Array of query results 
	 */
	public function executeQuery($query, $par = NULL)
	{
		try
		{
			$hash = $this->prepare($query);
			if(!empty($par))
			{
				foreach($par as $key => $value)
				{
					$this->bind($hash, $key, $value);
				}
			}
			$this->execute($hash, true);
			$result = $this->fetch($hash);
			$this->unsetStmt($hash);
		}
		catch(PDOException $e)
		{
			$extra = "";
			if (starts_with($query, "INSERT") || starts_with($query, "UPDATE"))
			{
				$extra = " - DB::executeQuery expects results";
			}
			$this->prep_stmts[$hash]->debugDumpParams();
			$this->error($hash, $e->getMessage().$extra);
		}
		
		return $result;
	}

	/**
	 * Executes a query that returns a single value.
	 * 
	 * @param string $query SQL query to execute
	 */
	public function getValue($query, $default = null)
	{
		try
		{
			$hash = $this->prepare($query);
			$this->execute($hash, true);
			$result = $this->fetch($hash);
			$this->unsetStmt($hash);
		}
		catch(PDOException $e)
		{
			$extra = "";
			if (starts_with($query, "INSERT") || starts_with($query, "UPDATE"))
			{
				$extra = " - DB::getValue expects results";
			}
			$this->prep_stmts[$hash]->debugDumpParams();
			$this->error($hash, $e->getMessage().$extra);
		}
			
		if (is_null($default) && count($result)==0) trigger_error("getValue query returned zero rows: $query", E_USER_ERROR);
		if (is_null($default) && count($result)>1) trigger_error("getValue query returned more than one row: $query", E_USER_ERROR);
		@$result = $result[0];
		if (is_null($default) && count($result)==0) trigger_error("getValue query returned zero columns: $query", E_USER_ERROR);
		if (is_null($default) && count($result)>1) trigger_error("getValue query returned more than one column: $query", E_USER_ERROR);
		@$result = array_values($result);
		@$result = $result[0];
		if ($result=="" || is_null($result)) return $default;
		return $result;
	}

	///Returns the query result as an array (works only for queries that return one value per row)
	public function getValues($query)
	{
		$output = array();
		$res = $this->executeQuery($query);
		foreach($res as $row)
		{
			if (count($row)!=1)
			{
				trigger_error("getValues query returned row with count!=1: ".implode("\t", $row), E_USER_ERROR);
			}
			foreach($row as $key => $value)
			$output[] = $value;
		}
		
		return $output;
	}
	
	///Returns all possible enum values of an enum column.
	public function getEnum($table, $column)
	{
		$res = $this->executeQuery("DESCRIBE $table $column");
		$type = $res[0]['Type'];
		$type = strtr($type, array("enum('"=>"", "')"=>""));
		return explode("','", $type);
	}
	
	///Returns the ID of a table entry (returns -1 or throws an error if the entry is not found)
	public function getId($table, $column, $value, $error_if_not_found = true)
	{
		$id = $this->getValue("SELECT id FROM {$table} WHERE {$column}='{$value}'", -1);
		if ($id==-1 && $error_if_not_found)
		{
			trigger_error("No entry in table {$table} where {$column} is '{$value}'", E_USER_ERROR);
		}
		return $id;
	}
	
	/**
	 * Begin of a transaction.
	 */
	public function beginTransaction()
	{
		$this->connection->beginTransaction();	
	}
	
	/**
	 * Rollback within a transaction.
	 */
	public function rollBack()
	{
		$this->connection->rollBack();	
	}
    
	/**
	 * End of a transaction.
	 */
	public function endTransaction()
	{
		try
		{
			$this->connection->commit();
		}
		catch(PDOException $e)
		{
			$this->error($hash, $e->getMessage());
		}
	}
	
	/**
	 * Removes a specific statement from the storage.
	 * 
	 * @param string $hash Identifier for the prepared statement that should be removed.
	 */
	public function unsetStmt($hash)
	{
		unset($this->prep_stmts[$hash]);
		unset($this->log_messages[$hash]);
	}
	
	/**
	 * Fetches results of a specific query
	 * 
	 * @param string $hash
	 * @return array
	 */
	public function fetch($hash)
	{
		return $this->prep_stmts[$hash]->fetchAll(PDO::FETCH_ASSOC);
	}
	
	/**
	 * Binds specific variables to a prepared statement.
	 * 
	 * @param string $hash Identifies the prepared statement
	 * @param string $id Identifier of the column
	 * @param mixed $value Value for the column
	 * @param array $null This array contains strings that are synonyms for the NULL value (array('null1',...)
	 */
	public function bind($hash, $id, $value, $null = array())
	{
		if(in_array($value, $null))
		{
			$value = null;
		}
		
		try
		{
			$this->log_messages[$hash]['bind'][] = "$id = $value;";
			$this->prep_stmts[$hash]->bindValue(":".$id, $value);
		}
		catch(PDOException $e)
		{
			$this->error($hash, $e->getMessage());
		}
	}

	/**
	 * Inserts a dataset into DB and returns the ID of the last dataset inserted.
	 * 
	 * @param string $hash identifies the prepared statement
	 * @return int ID of the last dataset Inserts
	 */
	public function insertGetID($hash)
	{
		$this->execute($hash, DB::PERSISTENT_STORAGE);
		$last_id = $this->connection->lastInsertId();

		return $last_id;
	}
	
	///Returns the ID of the last insert
	public function lastInsertId()
	{
		return $this->connection->lastInsertId();
	}
	
	/**
	 * Logs and writes error messages on std-out.
	 * 
	 * @param string $hash Identifies sql-statement
	 * @param string $message Error message
	 */
	protected function error($hash, $message)
	{
		//generate message
		$messages = [];
		$messages[] = trim($message);
		if (isset($this->log_messages[$hash]))
		{
			if (isset($this->log_messages[$hash]['query']))
			{
				$messages[] = "QUERY:";
				$messages[] = trim($this->log_messages[$hash]['query']);
			}
			if (isset($this->log_messages[$hash]['bind']))
			{
				$messages[] = "BOUND PARAMETERS:";
				foreach($this->log_messages[$hash]['bind'] as $binding)
				{
					$messages[] = trim($binding);
				}
			}
		}
		
		//throw message
		trigger_error(implode("\n", $messages)."\n", E_USER_ERROR);	
	}
}

?>