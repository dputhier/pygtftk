MAKEFILE=Makefile

.PHONY: help doc pylint pylintshort \
nose unittest install install_user test bats bats_cmd \
test_cmd clean_install prepare_pip


#------------------------------------------------------------------
# PYLINT arguments
#------------------------------------------------------------------

## Disable some testing in pylint
# W0123 : eval-used
# W0122 : exec-used
PYLINT_DISABLE=W0123,W0122

## Fix some issues with numpy,bx and pylint
# --extension-pkg-whitelist=numpy,bx

## Pylint arguments
PYLINT_ARGS= -d $(PYLINT_DISABLE) --extension-pkg-whitelist=numpy

#------------------------------------------------------------------
# Other variable
#------------------------------------------------------------------

CRAN_PACK_LIST=''
BIOC_PACK_LIST=''

#------------------------------------------------------------------
# Targets
#------------------------------------------------------------------


help:
	@echo ""
	@echo "- Available targets:"
	@echo ""
	@perl -ne 'if(	/^(\w+):/){print "\t",$$1,"\n"}' $(MAKEFILE)
	@echo ""
	@echo ""
	@echo "- More information:"
	@echo ""
	@echo -e "\tbats_cmd Usage: make bats_cmd CMD=select_by_key"
	@echo ""
	@echo -e "\ttest_para Usage: make test_para -j 10"

doc:
	@rm -Rf docs/manual/build/
	@rm -Rf  dist build *.egg*
	@cd docs/manual; make html; cd ../..;
	@echo ">>>> Docs is available in: docs/manual/build/html/index.html"

pylint:
	@find . -name "*.py"  -exec pylint $(PYLINT_ARGS) {} \;

pylintshort:
	@rm -Rf build dist gtftk.egg-info
	@find . -name "*.py"  -exec pylint $(PYLINT_ARGS) {} \; 2>/dev/null |perl -ne  "print if(/(Your code has been rated)|(\*\*\*\* Module)/)"

nose:
	@cd /tmp; mkdir -p gtftk_test; cd gtftk_test; a=`python -c "import os,pygtftk; print(os.path.dirname(pygtftk.__file__))"`; cd $$a ; for i in `find . -name "*.py" | perl -ne  'print unless(/(setup)|(plugin)|(libgtftk.py)|(__)/)'`; do echo "================="; echo $$i; nosetests --with-doctest $$i;   done

install:
	@rm -Rf build  dist pygtftk.egg-info ~/.gtftk; \
	python setup.py install; rm -Rf build  dist pygtftk.egg-info

install-dbg:
	@rm -Rf build  dist pygtftk.egg-info ~/.gtftk; \
	python-dbg setup.py install; rm -Rf build  dist pygtftk.egg-info

test:
	@make bats
bats:
	@gtftk -p > gtftk_test.bats
	@bats -t gtftk_test.bats
	@rm -f gtftk_test.bats


bats_cmd:
	@gtftk -p|perl -npe 's/^ +//;' |perl -ne 'BEGIN{$$/="\n}"}{print $$_ if (/\@test +"$(CMD)/)}'  > gtftk_test.bats.sub
	@bats -t gtftk_test.bats.sub
	@rm -f gtftk_test.bats.sub

test_cmd:
	make bats_cmd CMD=$(CMD)

%.bats:
	@gtftk -l > prgm_list.txt; gtftk -p > test_list.txt; for i in $$(cat prgm_list.txt); do  cat test_list.txt | grep -E "@test \"$$i" -A 3  | grep -v "^\-\-$$" > $$i.bats; done

%.completed : %.bats
	@bats -t $<
	@echo "completed" > $@

OUTPUT = $(eval OUTPUT := $$(shell gtftk -l 2>/dev/null))$(OUTPUT)

OUTPUT2 = $(addsuffix .completed, $(OUTPUT))

test_para: $(OUTPUT2)


clean:
	@make bats_cmd CMD=clean
	@rm -rf pygtftk.egg-info build airway_love.txt* ENCFF630HEX_Total_RNAseq_K562_count_mini.txt STDIN.e* STDIN.o* dist cmd_list.txt example_list.txt tmp_list.txt simple.chromInfo prgm_list.txt test_list.txt *.bats *.completed *mini_real* heatmap_* tx_classes* *~ \#* hh profile_* toto tott;  cd docs/manual/; make clean; cd ../..; find . -type f -name '*~' -exec rm -f '{}' \; git checkout docs/manual/source/conf.py pygtftk/version.py

check_cmd_has_example:
	@for i in $$(gtftk -l); do if grep -q  "^$$i" docs/manual/source/presentation.rst; then echo "" >/dev/null; else echo $$i; fi; done

check_example_has_cmd:
	@gtftk -l > cmd_list.txt
	@grep "~~" -B 1 docs/manual/source/presentation.rst | grep -v "^$$" | grep -v " " |grep -v "^~~" | grep -v  "^\-\-" > example_list.txt
	@for i in $$(cat example_list.txt); do if $$(cat cmd_list.txt  | grep -q "^$$i")  ; then echo "" >/dev/null; else echo $$i; fi; done
	@#rm -f cmd_list.txt example_list.txt tmp_list.txt

prepare_pip:
	@touch pypi_release_in_progress
	@rm -rf build/ dist/ pygtftk.egg-info/
	@python setup.py sdist
	@python setup.py bdist_wheel
	@rm -f pypi_release_in_progress

send_to_pypi:
	@twine upload --repository-url https://test.pypi.org/legacy/ dist/* --verbose

nb_test:
	@gtftk -p| perl -ne 'BEGIN{$$/="{"}{/\@test\s+"(\w+)_\d+"/; print $$1,"\n"}'| sort | uniq -c | sort -nr
