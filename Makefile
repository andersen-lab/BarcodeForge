# Description: Makefile for Nextflow pipeline development

.PHONY: clean md5sum test run-test

clean:
	rm -rf work/ results .nextflow.log* .nextflow/

md5sum:
	bash assets/test/md5.sh

run-test-newick: clean
	nextflow run main.nf -profile test,mamba -params-file assets/test/input/params_newick.json
	$(MAKE) md5sum

run-test-nexus: clean
	nextflow run main.nf -profile test,mamba -params-file assets/test/input/params_nexus.json
	$(MAKE) md5sum

test:
	$(MAKE) run-test-newick
	$(MAKE) run-test-nexus
