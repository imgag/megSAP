<?php

include("framework.php");
//require_once("all.php");

//##################################################################################
start_test("jobs_wait");

$qstat = jobsWait(array(17753,17754));
check($qstat, array(17753 => "finished", 17754 => "finished"));

end_test();

//##################################################################################
start_test("job_exit_status");

check(jobStatus(15764), 'finished');
check(jobStatus(15735), 'job error');
check(jobStatus(15768), 'finished');

end_test();
?>
