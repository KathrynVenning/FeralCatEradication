.PHONY: \
		check \
		clean \
		coverage \
		linter \
		mutants \
		results \
		setup \
		tests

define lint
	R -e "library(lintr)" \
	  -e "lint_dir('src', linters = with_defaults(line_length_linter(120)))"
endef

check: linter

clean:
	rm --force --recursive FeralCatEradication.Rcheck
	rm --force --recursive reports/figures
	rm --force FeralCatEradication_*.tar.gz
	rm --force FeralCatEradication/NAMESPACE
	rm --force Rplots.pdf

coverage: setup
	R -e "covr::package_coverage('FeralCatEradication')"

linter:
	$(lint)
	$(lint) | grep -e "\^" && exit 1 || exit 0

mutants: tests
	@echo "🙁🏹 No mutation testing on R 👾🎉👾"

results: src/FeralCatEradication.R
	mkdir reports/figures/ --parents
	Rscript src/FeralCatEradication.R

setup:
	R -e "devtools::document('FeralCatEradication')" && \
	R CMD build FeralCatEradication && \
	R CMD check FeralCatEradication_0.1.2.tar.gz && \
	R CMD INSTALL FeralCatEradication_0.1.2.tar.gz
	
tests:
	R -e "testthat::test_dir('tests/testthat/', report = 'summary', stop_on_failure = TRUE)"
