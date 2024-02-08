help:
	@echo "Main targets:"
	@echo "  test_functions  - perform function tests." 
	@echo "  test_tools      - perform tool tests."
	@echo "  test_tools_db   - perform DB tool tests." 
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
	@echo "  find_missing_tests - Checks for tools that do not have a test."
	@echo "  todos              - Checks for todos in the code"
	@echo "  pull               - Pull latest version from GitHub"
	
pull:
	git pull --recurse-submodules
	git describe --tags > megSAP_tag.txt
	git status

test_functions: dummy
	@cd test && find . -name "test_*.php" | sort | xargs -l1 php | egrep "^(FINISHED|FAILED)"
	
test_tools: dummy
	@cd test && find . -name "tool_test_*.php" | sort | xargs -l1 php | egrep "^(FINISHED|FAILED)"

test_tools_primer: dummy
	@cd test && find . -name "tool_test_*.php" | xargs grep -l "/Primer/" | sort | xargs -l1 php | egrep "^(FINISHED|FAILED)"

test_tools_ngs: dummy
	@cd test && find . -name "tool_test_*.php" | xargs grep -l "/NGS/" | sort | xargs -l1 php | egrep "^(FINISHED|FAILED)"

test_tools_tools: dummy
	@cd test && find . -name "tool_test_*.php" | xargs grep -l "/Tools/" | sort | xargs -l1 php | egrep "^(FINISHED|FAILED)"

test_tools_db:
	@cd test && find . -name "tool_test_db*.php" | sort | xargs -l1 php | egrep "^(FINISHED|FAILED)"

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

test_pipeline_tl: dummy
	@cd test/data_trio_longread && make all

test_pipeline_ml: dummy
	@cd test/data_multi_longread && make all

test_all: dummy
	(cd test && find . -name "test_*.php" | sort | xargs -l1 php && echo "DONE") > f.log 2>&1
	(cd test && find . -name "tool_test_*.php" | sort | xargs -l1 php && echo "DONE") > t.log 2>&1 &
	make test_pipeline_a > p_a.log 2>&1 &
	make test_pipeline_x > p_x.log 2>&1 &
	make test_pipeline_s > p_s.log 2>&1 &
	make test_pipeline_t > p_t.log 2>&1 &
	make test_pipeline_m > p_m.log 2>&1 &
	make test_pipeline_r > p_r.log 2>&1 &
	make test_pipeline_c > p_c.log 2>&1 &
	make test_pipeline_l > p_l.log 2>&1 &
	make test_pipeline_tl > p_tl.log 2>&1 &
	make test_pipeline_ml > p_ml.log 2>&1 &

test_all_status:
	@clear
	@tail -v -n3 *.log
	@echo ""
	@echo "### WARNINGS ###"
	@egrep -a -i "WARNING" *.log || :
	@echo ""
	@echo "### ERRORS ###"
	@egrep -a -i "ERROR|FAILED" *.log | grep -v "Medelian errors:" || :

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
	php src/Tools/find_missing_tests.php
	
find_unused_tools: dummy
	php src/Tools/find_unused_tools.php -ngsbits ../ngs-bits/ -megsap . > unused_tools.txt

todos:
	find . -name "*.php" | xargs grep -i "//todo" 

find_php_warnings_in_tests:
	find test/ -type f | xargs grep -i "PHP warning" || true
	find test/ -type f | xargs grep -i "PHP notice" || true

swap_settings:
	mv settings.ini settings.ini.swap
	mv settings_nightly.ini settings.ini
	mv settings.ini.swap settings_nightly.ini

clean_ignored:
	git clean -Xn | cut -f3 -d' ' | egrep -v "settings.ini|megSAP_tag.txt" | xargs rm

dummy:
