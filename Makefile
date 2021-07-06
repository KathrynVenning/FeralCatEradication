.PHONY: results tests

define runLint
	R -e "library(lintr)" \
	  -e "lint_dir('src', linters = with_defaults(line_length_linter(100)))"
endef

clean:
	rm --force Rplots.pdf
	rm --force --recursive reports/figures

coverage:
	R -e "covr::package_coverage()"

lint:
	$(runLint)
	$(runLint) | grep -e "\^" && exit 1 || exit 0

results: src/FeralCatEradication.R src/matrixOperators.r
	mkdir reports/figures/ --parents
	Rscript src/FeralCatEradication.R

tests:
	R -e "testthat::test_dir('tests/testthat/', report = 'summary', stop_on_failure = TRUE)"

install:
	R -e "devtools::document()" &&
	R CMD build . && \
    R CMD check dimorfismo_0.1.0.tar.gz && \
    R CMD INSTALL dimorfismo_0.1.0.tar.gz
