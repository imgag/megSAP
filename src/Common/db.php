<?php

require_once("genomics.php");

class DB
{
	protected static $instance = null;
	protected $connection = null;
	protected $prep_stmts = array();
	protected $stmt_counter = 0;
	protected $log = false;
	protected $path_to_log = '';
	protected $log_messages = array();
	protected $db_info = "";
	const PERSISTENT_STORAGE = TRUE;
	
	/**
	 * Returns current instance of the DB connection and - if
	 * not yet available - creates a new one.
 	 * 
	 * @param $db string Name of database to access
	 * @return PDO_Connection
	 */
	public static function getInstance($db)
	{
		$host = "mysql:host=".get_db($db, 'db_host');
		$name = "dbname=".get_db($db, 'db_name');
		
		//create a new connection if no connection exists OR if the DB host/name have changed
		if(is_null(self::$instance) || self::$instance->getDBInfo()!=$host.";".$name || !self::$instance->ping())
		{
			$user = get_db($db, 'db_user');
			$pass = get_db($db, 'db_pass');
			$log = get_db($db, 'db_log');
			$path_to_log = get_db($db, 'db_log_path');
			self::$instance = new self($host, $name, $user, $pass, $log, $path_to_log);
		}
		
		return self::$instance;
	}
	
	protected function __construct($host, $name, $user, $pass, $log, $path_to_log)
	{
		// setup connection
		$this->connection = new PDO($host.";".$name, $user, $pass, array(PDO::MYSQL_ATTR_INIT_COMMAND => "SET NAMES utf8"));
		$hash = $this->prepare("SET SESSION sql_mode = 'STRICT_TRANS_TABLES';");
		
		try 
		{
			$this->execute($hash);
			
			// error reporting
			$this->log = $log;
			$this->path_to_log = $path_to_log;
			$this->connection->setAttribute(PDO::ATTR_ERRMODE, PDO::ERRMODE_EXCEPTION);		
			
			$this->db_info = $host.";".$name;
		}
		catch(PDOException $e)
		{
			$this->error($hash, $e->getMessage());
		}
	}
	
	///custom ping functoin
	function ping() 
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
	
	
	///Returns the DB info string ($host.";".$name)
	function getDBInfo()
	{
		return $this->db_info;
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
		$this->log_messages[$hash] = array('sql' => "", 'par' => array());
		
		try
		{
			$this->log_messages[$hash]['sql'] = "SQL: ".$stmt;
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
			
			//normal log message (!= error)
			if($this->log)
			{
				if(stripos($this->log_messages[$hash]['sql'],'insert into')!==false || 
				   stripos($this->log_messages[$hash]['sql'],'delete from')!==false || 
				   stripos($this->log_messages[$hash]['sql'],'update')!==false)
				{
					$this->log($hash);
				}
			}
			
			if(!$persist)
			{
				$this->unsetStmt($hash);
			}
		}
		catch(PDOException $e)
		{
			print $e->getMessage();
			//error log message
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
	 */
	public function executeStmt($sql, $par = NULL, $log = FALSE)
	{
		$this->log = $log;
		$hash = $this->prepare($sql);
		if(!empty($par))
		{
			foreach($par as $key => $value)
			{
				$this->bind($hash, $key, $value);
			}			
		}
		$this->execute($hash, false);
	}
	
	/**
	 * Executes a query without binding of parameters and returns the result
	 * 
	 * @param string $query SQL query to execute
	 * @param boolean $log Sets logging
	 * @return array Array of query results 
	 */
	public function executeQuery($query, $par = NULL, $log = FALSE)
	{
		$this->log = $log;
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
				$extra = " - EXECUTEQUERY expects results";
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
				$extra = " - EXECUTEQUERY expects results";
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
			$this->log_messages[$hash]['par'][] = "$id = $value;";
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
		$this->log_messages[$hash]['sql'] = "ERROR: ".$message."; ".$this->log_messages[$hash]['sql'];
		$this->log($hash);
	}
	
	/**
	 * Logs all kind of messages to a log-file
	 * 
	 * @param string $hash Identifies sql-statement
	 */
	public function log($hash)
	{
		//generate message
		$user = get_current_user()." (SYSTEM)";
		if(isset($_SESSION['user']))
		{
			$tmp = $_SESSION['user'];
			$user = $tmp->user_id." (WEBGUI)";
		}

		$message = "";
		if (isset($this->log_messages) && isset($this->log_messages[$hash]))
		{
			$message = get_timestamp()."\t".$user."\t".str_replace(array("\r\n", "\n", "\r"), " ", $this->log_messages[$hash]['sql'])."; PAR: ".implode(" ", $this->log_messages[$hash]['par'])."\n";
		}
		
		//if log set to true write to file
		if($this->log)
		{
			$pathtolog = $this->path_to_log;
			file_put_contents($pathtolog, $message, FILE_APPEND);
		}
		
		//handle real errors
		if(isset($_SESSION['user']) && strpos($message, "ERROR") !== FALSE)
		{
			lib::setError($message);
			lib::sendTo();
		}
		elseif(strpos($message, "ERROR") !== FALSE)
		{
			trigger_error($message, E_USER_ERROR);			
		}
	}
}


