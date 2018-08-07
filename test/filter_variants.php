<?php

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/all.php");

list($stdout, $stderr) = exec2(get_path("ngs-bits")."/VariantFilterAnnotations ".implode(" ", array_slice($argv, 1)));

?>