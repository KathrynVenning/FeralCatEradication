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
	  -e "lint_dir('R', linters = with_defaults(line_length_linter(120)))" \
	  -e "lint_dir('tests', linters = with_defaults(line_length_linter(120)))" \
	  -e "lint_dir('src', linters = with_defaults(line_length_linter(120)))"
endef

check:
	R -e "library(styler)" \
	  -e "resumen <- style_dir('R')" \
	  -e "resumen <- rbind(resumen, style_dir('src'))" \
	  -e "resumen <- rbind(resumen, style_dir('tests'))" \
	  -e "any(resumen[[2]])" \
	  | grep FALSE

clean:
	rm --force --recursive FeralCatEradication.Rcheck
	rm --force --recursive reports/figures
	rm --force FeralCatEradication_*.tar.gz
	rm --force FeralCatEradication/NAMESPACE
	rm --force Rplots.pdf

coverage: setup
	R -e "cobertura <- covr::file_coverage(c('R/feral_cat.R'), c('tests/testthat/test_feral_cat.R'))" \
	  -e "covr::codecov(covertura=cobertura, token='d40cba41-8ee3-414d-9e04-581d33a42b62')"

format:
	R -e "library(styler)" \
	  -e "style_dir('R')" \
	  -e "style_dir('src')" \
	  -e "style_dir('tests')"

linter:
	$(lint)
	if $(lint) | grep -e "\^" ; then exit 1 ; else exit 0 ; fi

mutants: tests
	@echo "🙁🏹 No mutation testing on R 👾🎉👾"

results: src/FeralCatEradication.R
	mkdir reports/figures/ --parents
	Rscript src/FeralCatEradication.R

setup:
	R -e "devtools::document()" && \
	R CMD build . && \
	R CMD check FeralCatEradication_0.1.3.tar.gz && \
	R CMD INSTALL FeralCatEradication_0.1.3.tar.gz
	
tests:
	R -e "testthat::test_dir('tests/testthat/', report = 'summary', stop_on_failure = TRUE)"
