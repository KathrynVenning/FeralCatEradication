.PHONY: \
		clean \
		coverage \
		install \
		linter \
		mutants \
		results \
		tests

define lint
	R -e "library(lintr)" \
	  -e "lint_dir('src', linters = with_defaults(line_length_linter(100)))"
endef

clean:
	rm --force Rplots.pdf
	rm --force --recursive reports/figures

coverage:
	R -e "covr::package_coverage()"

install:
	R -e "devtools::document()" &&
	R CMD build . && \
    R CMD check FeralCat_0.1.0.tar.gz && \
    R CMD INSTALL FeralCat_0.1.0.tar.gz

linter:
	$(lint)
	$(lint) | grep -e "\^" && exit 1 || exit 0

mutants:
	@echo "🙁🏹 No mutation testing on R 👾🎉👾"

results: src/FeralCatEradication.R src/matrixOperators.r
	mkdir reports/figures/ --parents
	Rscript src/FeralCatEradication.R

tests:
	R -e "testthat::test_dir('tests/testthat/', report = 'summary', stop_on_failure = TRUE)"
