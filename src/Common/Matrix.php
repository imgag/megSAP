<?php

require_once("functions.php");

/**
	@brief A convenient two-dimensional array implementation, i.e. a matrix.
	
	@note The content is stored row-wise. Thus, row-rise operations are much faster than column-wise operations. It might be good to transpose the matrix if you perfome a lot of column-wise operations.
	
	@ingroup base
*/
class Matrix
{
	private $rows = 0;
	private $cols = 0;
	private $data = array();
	private $headers = array();
	private $comments = array();
	
	/// Constructor.
	function __construct()
	{
	}

	/// Returns the column count
	function cols()
	{
		return $this->cols;
	}

	/// Returns the row count
	function rows()
	{
		return $this->rows;
	}
	
	/// Returns the comment count
	function comments()
	{
		return count($this->comments);
	}

	/// Returns a data value correspinding to the given row and column index.
	function get($row, $col)
	{
		// Check valid indices
		if ($col >= $this->cols)
		{
			trigger_error("Internal error: Column index '$col' is not available.", E_USER_ERROR);
		}
		
		if ($row >= $this->rows)
		{
			trigger_error("Internal error: Row index '$row' is not available.", E_USER_ERROR);
		}
		
		return $this->data[$row][$col];
	}
	
	/// Returns the column index corresponding to a column name.
	/// @note If the colum name is not found, 'false' is returned or an error is thrown.
	/// @note If the colum name is not found more than once, an error is thrown.
	function getColumnIndex($name, $containing = false, $error_if_not_found = true)
	{
		$output = -1;
		foreach($this->headers as $index => $header)
		{
			if (!$containing && $header == $name)
			{
				if ($output!=-1)
				{
					trigger_error("Internal error: Column name '$name' is used more than once.", E_USER_ERROR);
				}
				$output = $index;
			}
			if ($containing && contains($header,$name))
			{
				if ($output!=-1)
				{
					trigger_error("Internal error: Column name '$name' is used more than once.", E_USER_ERROR);
				}
				$output = $index;
			}
		}
		
		if ($output==-1)
		{
			if ($error_if_not_found)
			{
				trigger_error("Internal error: Column name '$name' is not available.", E_USER_ERROR);
			}
			else
			{
				return false;
			}
		}
		
		return $output;
	}
	
	/// Returns a column by index.
	function getCol($col)
	{
		// Check valid indices
		if ($col >= $this->cols)
		{
			trigger_error("Internal error: Column index '$col' is not available.", E_USER_ERROR);
		}
		
		$output = array();
		for($i=0; $i<$this->rows; ++$i)
		{
			$output[] = $this->data[$i][$col];
		}
		
		return $output;
	}

	/// Returns a row by index.
	function getRow($row)
	{
		if ($row >= $this->rows)
		{
			trigger_error("Internal error: Row index '$row' is not available.", E_USER_ERROR);
		}
		
		return $this->data[$row];
	}

	/// returns Comments
	function getComments() 
	{
		return $this->comments;
	}

	/// Sets the data value with the given row and column index.
	function set($row, $col, $value)
	{
		// Check valid indices
		if ($col >= $this->cols)
		{
			trigger_error("Internal error: Column index '$col' is not available.", E_USER_ERROR);
		}
		
		if ($row >= count($this->data))
		{
			trigger_error("Internal error: Row index '$row' is not available.", E_USER_ERROR);
		}
		
		$this->data[$row][$col] = $value;
	}

	/// Sets the values of a column.
	function setCol($column, $values)
	{
		// Check row count
		if (count($values) != $this->rows)
		{
			trigger_error("Internal error: Cannot set a column with differing row count (is '".count($values)."', should be '".$this->rows."').", E_USER_ERROR);
		}
		
		for($i=0; $i<$this->rows; ++$i)
		{
			$this->data[$i][$column] = $values[$i];
		}
	}

	/// Insert a column at given index and move all other columns to the right.
	function insertCol($column_index, $values, $header = "", $desc = "")
	{
		// Check row count
		if (count($values) != $this->rows)
		{
			trigger_error("Internal error: Cannot set a column with differing row count (is '".count($values)."', should be '".$this->rows."').", E_USER_ERROR);
		}
		
		//updated header
		$tmp = array_splice($this->headers,$column_index);
		$this->headers[$column_index] = $header;
		$this->headers = array_merge($this->headers,$tmp);
		
		//update comments
		if(!empty($desc))	$this->comments[] = "#DESCRIPTION=$header=$desc";
		
		//update data
		for($i=0; $i<$this->rows; ++$i)
		{
			$tmp = array_splice($this->data[$i],$column_index);
			$this->data[$i][$column_index] = $values[$i];
			$this->data[$i] = array_merge($this->data[$i],$tmp);
		}
		
		//increase number of columns
		++$this->cols;
	}
	
	/// Moves a column to the end of the matrix
	function moveColToEnd($i_col)
	{
		$temp_header = $this->getHeaders()[$i_col];
		$temp_col = $this->getCol($i_col);
		$temp_desc = "";
		foreach($this->getComments() as $comment)
		{
			if(starts_with($comment, "#DESCRIPTION="))
			{
				list(,$name,$desc) = explode("=",$comment);
				
				if($name == $temp_header)
				{
					list(,,$temp_desc) = explode("=",$comment); //3d part of comment field contains description e.g. '#DESCRIPTION=BLA=WHAT I DESCRIBE'
					$this->removeComment($name);
					break;
				}
			}
		}
		
		$this->removeCol($i_col);
		
		$this->addCol($temp_col, $temp_header, $temp_desc);
	}
	
	/// Sets the values of a row.
	function setRow($row, $values)
	{
		// Check column count
		if (count($values) != $this->cols)
		{
			trigger_error("Internal error: Cannot set a row with differing column count (is '".count($values)."', should be '".$this->cols."').", E_USER_ERROR);
		}

		$this->data[$row] = $values;
	}
	
	/// Sets Comments
	function setComments($lines) 
	{
		$this->comments = $lines;
	}
	
	/// Adds a Comment
	function addComment($line) 
	{
		$this->comments[] = $line;
	}
	
	/// Adds a Comment (at the beginning)
	function prependComment($line)
	{
		array_unshift($this->comments, $line);
	}
	
	/// Removes a Comment
	function removeComment($comment_text, $containing = false)
	{
		
		if(!$containing)
		{
			for($i=0;$i<count($this->comments);$i++)
			{
				if(!$containing)
				{
					if($comment_text == $this->comments[$i])
					{
						array_splice($this->comments,$i,1);
						break;
					}
				}
			}
		}
		else //if containing is select, start from the end because it can affect multiple comments
		{
			for($i=(count($this->comments)-1); $i>=0; --$i)
			{
				if(strpos($this->comments[$i], $comment_text) !== false)
				{
					array_splice($this->comments,$i,1);
				}
			}
		}
	}
	
	/// Appends a column.
	function addCol($column, $name = "", $desc = "")
	{
		// Check row count
		if ($this->rows!=0 && count($column) != $this->rows)
		{
			trigger_error("Internal error: Cannot add a column with differing row count (is '".count($column)."', should be '".$this->rows."').", E_USER_ERROR);
		}
		
		if ($this->rows==0)
		{
			$this->rows = count($column);
		}
		
		$num_rows = max(count($column), $this->rows);
		for($i=0; $i<$num_rows; ++$i)
		{
			$this->data[$i][] = $column[$i];
		}
		
		++$this->cols;
		$this->headers[] = $name;
		if(!empty($desc))	$this->comments[] = "#DESCRIPTION=$name=$desc";
	}
	
	/// Appends a row.
	function addRow($row, $filename = null)
	{
		// Set column count
		if ($this->rows==0)
		{
			$this->cols = count($row);
			$this->headers = array_pad($this->headers, count($row), "");
		}
		
		// Check column count
		if (count($row) != $this->cols)
		{
			trigger_error("Internal error: Cannot add a row with differing column count (is '".count($row)."', should be '".$this->cols."')".(isset($filename) ? " in file {$filename}" : ""), E_USER_ERROR);
		}

		++$this->rows;
		$this->data[] = $row;
	}

	/// Removes a column by index.
	function removeCol($column)
	{
		if ($column >= $this->cols)
		{
			trigger_error("Internal error: Column index '$column' is not available.", E_USER_ERROR);
		}

		for($i=0; $i<$this->rows; ++$i)
		{
			array_splice($this->data[$i], $column, 1);
		}
		
		--$this->cols;
		$name = array_splice($this->headers, $column, 1)[0];
		
		$tmp_c = array();
		foreach($this->comments as $c)
		{
			if(starts_with($c,"#DESCRIPTION"))
			{
				list($t,$n,$d) = explode("=",$c);
				if($n==$name)	continue;
			}
			$tmp_c[] = $c;
		}
		$this->comments = $tmp_c;
	}
	
	//Removes all columns with "name" as header
	function removeColByName($name, $containing = false)
	{
		$col_indices = array();
		foreach($this->headers as $index => $col_name)
		{
			if(!$containing && $col_name == $name)
			{
				$col_indices[] = $index;
			}
			if($containing && strpos($col_name, $name) !== false)
			{
				$col_indices[] = $index;
			}
		}
		
		if(empty($col_indices)) return;
		
		for($i=count($col_indices)-1;$i>=0;--$i)
		{
			$this->removeCol($col_indices[$i]);
		}
	}

	/// Removes a row by index.
	function removeRow($row)
	{
		if ($row >= $this->rows)
		{
			trigger_error("Internal error: Row index '$row' is not available.", E_USER_ERROR);
		}
		
		--$this->rows;
		array_splice($this->data, $row, 1);
	}
	
	/// Returns the headers.
	function getHeaders()
	{
		return $this->headers;
	}

	/// Returns the headers.
	function getHeader($column)
	{
		if ($column >= $this->cols)
		{
			trigger_error("Internal error: Column index '$column' is not available.", E_USER_ERROR);
		}
		
		return $this->headers[$column];
	}

	/// Sets the header for a column.
	function setHeader($column, $name)
	{
		if ($column >= $this->cols)
		{
			trigger_error("Internal error: Column index '$column' is not available.", E_USER_ERROR);
		}
		
		$this->headers[$column] = $name;
	}
	
	/// Sets the header for a column.
	function setHeaders($names, $filename = null)
	{
		if ($this->cols==0)
		{
			$this->cols = count($names);
		}
		
		if (count($names) != $this->cols)
		{
			trigger_error("Internal error: Cannot set headers with differing column count (is '".count($names)."', should be '".$this->cols."')".(isset($filename) ? " in file {$filename}" : ""), E_USER_ERROR);
		}
		
		$this->headers = $names;
	}
	
	/// Loads a matrix from a TSV file and returns it.
	static function fromTSV($filename, $separator="\t", $comment="#")
	{
		$comments = array();
		$output = new Matrix();
		
		if(!file_exists($filename))	trigger_error("Could not find file '$filename'.",E_USER_ERROR);
		
		$handle = gzopen2($filename , "r");
		while (!feof($handle)) 
		{
			$line = nl_trim(fgets($handle));
			
			//skip comment line (but store them for headers)
			if (strlen($line)>0 && $line[0]==$comment)
			{
				$comments[] = substr($line, 1);
				continue;
			}
			
			//skip empty lines (sometimes at the end of the file)
			if ($line=="") continue;
			
			$output->addRow(explode($separator, $line), $filename);
		}
		gzclose($handle);
		
		//headers
		for ($i=0; $i<count($comments); ++$i)
		{
			$parts = explode($separator, trim($comments[$i]));
			if (count($parts)==$output->cols() || ($output->cols()==0 && $i==(count($comments)-1)))
			{
				$output->setHeaders($parts, $filename);
			}
			else
			{
				$output->addComment($comments[$i]);
			}
		}

		return $output;
	}
	
	/// Stores the matrix to a TSV file.
	function toTSV($filename, $separator="\t", $comment="#", $header = "#")
	{
		$handle = fopen2($filename , "w");
		
		//comments
		if (implode("", $this->comments)!="")
		{
			foreach($this->getComments() as $comment_line)
			{
				fwrite($handle , $comment.trim($comment_line)."\n");
			}
		}

		//header
		if (implode("", $this->headers)!="")
		{
			fwrite($handle , $header.implode($separator, $this->headers)."\n");
		}
		
		//content
		for($i=0; $i<$this->rows; ++$i)
		{
			fwrite($handle , implode($separator, $this->data[$i])."\n");
		}
		fclose($handle);
	}
	
	/// Prints basic infos about the matrix: Number of columns and rows, header, first and last lines.
	function printInfo($first_lines = 1, $last_lines = 0)
	{
		print "Columns: ".$this->cols()."\n";
		print "Rows   : ".$this->rows()."\n";
		print "Headers:\n";
		for ($i=0; $i<$this->cols(); ++$i)
		{
			print "  [$i]: ".$this->getHeader($i)."\n";
		}
		for ($r=0; $r<$first_lines; ++$r)
		{
			print "Row $r:\n";
			for ($c=0; $c<$this->cols(); ++$c)
			{
				print "  [$c]: ".$this->get($r, $c)."\n";
			}
		}
		for ($r=$this->rows-$last_lines; $r<$this->rows; ++$r)
		{
			print "Row $r:\n";
			for ($c=0; $c<$this->cols(); ++$c)
			{
				print "  [$c]: ".$this->get($r, $c)."\n";
			}
		}
	}
	
	///  Prints the whole matrix
	function printFull()
	{
		print implode("\t", $this->getHeaders())."\n";
		for ($i=0; $i<$this->rows(); ++$i)
		{
			print implode("\t", $this->getRow($i))."\n";
		}
	}
	
	/// Sorts the matrix according to a column
	function sort($column, $type = SORT_REGULAR, $order = SORT_ASC)
	{
		$col_data = $this->getCol($column);
		array_multisort($col_data, $type, $order, $this->data);
	}
	
	/// Resizes the matrix according to given parameters
	function resize($no_rows, $no_cols)
	{
		if($no_rows < 0 || $no_cols < 0)
		{
			trigger_error("Only positive row/col-numbers are allowed!", E_USER_ERROR);
		}
		
		if($this->rows()==0 && $this->cols()==0)
		{
			$this->addRow(array_fill(0,($no_cols-1),"."));	
		}
		
		if($this->rows() < $no_rows)
		{
			for($j=$this->rows(); $j<$no_rows; ++$j)
			{
				$this->addRow(array_fill(0,$this->cols(),"."));
			}
		}
		elseif($this->rows() > $no_rows)
		{
			for($j=($this->rows()-1); $j>=$no_rows; --$j)
			{
				$this->removeRow($j);
			}
		}
		
		if($this->cols() < $no_cols)
		{
			for($j=$this->cols(); $j<$no_cols; ++$j)
			{
				$this->addCol(array_fill(0,$this->rows(),"."));
			}
		}
		//only remove cols if there are rows
		elseif($this->cols() > $no_cols)
		{
			for($j=($this->cols()-1); $j>=$no_cols; --$j)
			{
				$this->removeCol($j);
			}
			
		}
	}
	
	
	/// Transposes the matrix (columns become rows and vice versa)
	function transpose()
	{
		// transpose data
		array_unshift($this->data, null);
		$this->data = call_user_func_array('array_map', $this->data);
		
		// swap col/row counts
		$tmp = $this->cols;
		$this->cols = $this->rows;
		$this->rows = $tmp;
		
		// clear headers
		$this->headers = array_fill(0, $this->cols, "");
	}
	
	/// Adds all columns of a second matrix based on two matching columns
	function join($column, $matrix, $column2)
	{
		// create mapping column: key => index in rhs matrix
		$col = $matrix->getCol($column2);
		$mapping = array_flip($col);
		
		// check that reference column contains unique values
		if (count($mapping)!=count($col))
		{
			trigger_error("Reference column contains non-unique values in Matrix::join function.", E_USER_ERROR);
		}
		
		// create column order array
		$order = array();
		$col = $this->getCol($column);
		foreach ($col as $value)
		{
			if (!isset($mapping[$value]))
			{
				trigger_error("The value '$value' is not contained in the reference column in Matrix::join.", E_USER_ERROR);
			}
			
			$order[] = $mapping[$value];
		}
		
		// add new columns with the determined order
		for ($i=0; $i<$matrix->cols(); ++$i)
		{
			if ($i==$column2)
			{
				continue;
			}
			
			$col = $matrix->getCol($i);
			$new_col = array();
			foreach($order as $index)
			{
				$new_col[] = $col[$index];
			}
			$this->addCol($new_col);
		}
	}
	
	///removes duplicates
	function unique()
	{
		$array_hashes = array();
		$tmp_data = array();
		
		for($i = 0; $i < $this->rows(); ++$i)
		{
			$hash = implode("_", $this->data[$i]);
			
			if(!isset($array_hashes[$hash]))
			{
				$array_hashes[$hash] = $hash;
				$tmp_data[] = $this->data[$i];				
			}			
		}
		
		$this->rows = count($tmp_data);
		$this->data = $tmp_data;
	}
	
	///returns the row data in the matrix
	function getData()
	{
		return $this->data;
	}
	
	///Reorders the columns according to the index array $columns. If $append_rest is true, all columns not mentioned in $columns are appended.
	function slice($columns, $append_rest=false)
	{
		if ($append_rest)
		{
			for ($i=0; $i<$this->cols(); ++$i)
			{
				if (!in_array($i, $columns))
				{
					$columns[] = $i;
				}
			}
		}
		
		//slice
		$this->headers = array_subset($this->headers, $columns);
		for ($i=0; $i<count($this->data); ++$i)
		{
			$this->data[$i] = array_subset($this->data[$i], $columns);
		}
		$this->cols = count($columns);
	}
	
}	
?>