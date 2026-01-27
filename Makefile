help:
	@echo "Main targets:"
	@echo "  test_functions  - perform function tests." 
	@echo "  test_tools      - perform tool tests."
	@echo "  test_pipeline_a - perform DNA amplicon pipeline test." 
	@echo "  test_pipeline_x - perform DNA shotgun pipeline test."
	@echo "  test_pipeline_t - perform DNA trio pipeline test."
	@echo "  test_pipeline_m - perform DNA multi-sample pipeline test."
	@echo "  test_pipeline_s - perform DNA somatic pipeline test."
	@echo "  test_pipeline_r - perform RNA pipeline test."
	@echo "  test_pipeline_c - perform cfDNA pipeline test."
	@echo "  test_pipeline_tl - perform DNA longread trio pipeline test."
	@echo "  test_pipeline_ml - perform DNA longread multi-sample pipeline test."
	@echo "  test_all        - perform all tests in parallel (functions, tools, pipelines)."
	@echo "  test_all_status - shows the output summary of 'test_all'."
	
	@echo "" 
	@echo "Auxilary targets:"
	@echo "  pull               - Pull latest version from GitHub"
	@echo "  find_missing_tests - Checks for tools that do not have a test."
	@echo "  find_unused_tools  - Checks for tools that are not used."
	@echo "  todos              - Checks for todos in the code"
	
	
	
pull:
	git pull --recurse-submodules
	git describe --tags > megSAP_tag.txt
	git status

test_functions: dummy
	@cd test && find . -name "test_*.php" | sort | xargs -l1 php | egrep "^(FINISHED|FAILED)"
	
test_tools: dummy
	@cd test && find . -name "tool_test_*.php" | sort | xargs -l1 php | egrep "^(FINISHED|FAILED)"

test_pipeline_a: dummy
	@cd test/data_amplicon && make all

test_pipeline_x: dummy
	@cd test/data_chrx && make all

test_pipeline_s: dummy
	@cd test/data_somatic && make all

test_pipeline_t: dummy
	@cd test/data_trio && make all

test_pipeline_m: dummy
	@cd test/data_multi && make all

test_pipeline_r: dummy
	@cd test/data_rna && make all

test_pipeline_c: dummy
	@cd test/data_cfdna && make all

test_pipeline_l: dummy
	@cd test/data_longread && make all

test_pipeline_l_all: dummy
	@cd test/data_longread && make all_3_tests

test_pipeline_tl: dummy
	@cd test/data_trio_longread && make all

test_pipeline_ml: dummy
	@cd test/data_multi_longread && make all

test_all: dummy
	php src/Tools/data_setup.php > test_setup.log
	(cd test && find . -name "test_*.php" | sort | xargs -l1 php && echo "DONE") > test_function.log 2>&1
	(cd test && find . -name "tool_test_*.php" | sort | xargs -l1 php && echo "DONE") > test_tools.log 2>&1 &
	make test_pipeline_a > test_p_a.log 2>&1 &
	make test_pipeline_x > test_p_x.log 2>&1 &
	make test_pipeline_s > test_p_s.log 2>&1 &
	make test_pipeline_t > test_p_t.log 2>&1 &
	make test_pipeline_m > test_p_m.log 2>&1 &
	make test_pipeline_r > test_p_r.log 2>&1 &
	make test_pipeline_c > test_p_c.log 2>&1 &
	make test_pipeline_l > test_p_l.log 2>&1 &
	make test_pipeline_tl > test_p_tl.log 2>&1 &
	make test_pipeline_ml > test_p_ml.log 2>&1 &

test_all_status:
	@clear
	@ls test_*.log | grep -v test_setup.log | xargs tail -v -n3
	@echo ""
	@echo "### WARNINGS ###"
	@egrep -a -i "WARNING" test_*.log | grep -v "command6_exit123" || :
	@echo ""
	@echo "### ERRORS ###"
	@egrep -a -i "ERROR|FAILED" test_*.log | egrep -v "Mendelian error rate:|command6_exit123" || :

test_clear:
	@cd test/data_amplicon && make clean
	@cd test/data_chrx && make clean
	@cd test/data_somatic && make clear
	@cd test/data_trio && make clear
	@cd test/data_multi && make clear
	@cd test/data_rna && make clean
	@cd test/data_cfdna && make clean
	@cd test/data_longread && make clean
	@cd test/data_trio_longread && make clean
	@cd test/data_multi_longread && make clean
	

test_clear_check:
	git status --ignored | grep "test/data_" | grep -v data_out

find_missing_tests: dummy
	php src/IMGAG/find_missing_tests.php
	
find_unused_tools: dummy
	php src/IMGAG/find_unused_tools.php -ngsbits ../ngs-bits/ -megsap .

todos:
	find . -name "*.php" | xargs grep -i "//todo" 

find_php_warnings_in_tests:
	find test/ -type f -and -name "*.log" -or -type f -and -name "*_output" | xargs egrep -i "PHP (warning|notice)" || true

swap_settings:
	mv settings.ini settings.ini.swap
	mv settings_nightly.ini settings.ini
	mv settings.ini.swap settings_nightly.ini

clean_ignored:
	git clean -Xn | cut -f3 -d' ' | egrep -v "settings.ini|megSAP_tag.txt" | xargs rm

dummy:
