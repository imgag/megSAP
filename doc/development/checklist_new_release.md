# megSAP release

1. Check that containers are ok
	
	> php src/IMGAG/container_status.php -check_md5

1. Check for unused tools/functions (on SRV005)

	> make find_unused_tools

1. Check for missing tests

	> make find_missing_tests

1. Execute all tests

	> make test_all

1. Make sure there are no PHP warnings in the tests

	> make find_php_warnings_in_tests

1. Make a test deployment on a clean Ubuntu using WSL on out test laptop.
1. Update the release version in `doc/install_unix.md`, commit and push.
1. Compile changelog for the new release:

	> git log [last-tag]..master --oneline  
	> git diff -w [last-tag] master src/Pipelines/
 
1. Create a new release on GitHub.
1. Update Zenodo links in `README.md`.
