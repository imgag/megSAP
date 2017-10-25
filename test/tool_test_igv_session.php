<?php

require_once("framework.php");

$name = "igv_session";
start_test($name);

$infiles = implode(" ",
	[ data_folder() . "vc_freebayes_out1.vcf.gz",
		data_folder() . "vc_freebayes_in.bam" ]);

$out1 = output_folder() . $name . "_out1.xml";
check_exec("php " . src_folder() . "/NGS/igv_session.php -in {$infiles} -out {$out1} -relative");
check_file($out1, data_folder() . $name . "_out1.xml");

end_test();
