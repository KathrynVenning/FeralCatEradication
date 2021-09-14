reports/figures/reduction_factor.jpg: src/plot_reduction_factor.R
	mkdir --parents $(@D)
	Rscript src/plot_reduction_factor.R

reports/figures/simulation.jpg: src/untreated_population.R
	mkdir --parents $(@D)
	Rscript src/untreated_population.R

reports/figures/constant_proportional_annual_cull.jpg: src/constant_proportional_annual_cull.R
	mkdir --parents $(@D)
	Rscript src/constant_proportional_annual_cull.R

reports/figures/monthly_time_serie_individuals.jpg: src/presentacion_210820.R
	mkdir --parents $(@D)
	Rscript src/presentacion_210820.R

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
	rm --force --recursive tests/testthat/_snaps
	rm --force FeralCatEradication_*.tar.gz
	rm --force NAMESPACE
	rm --force Rplots.pdf

coverage: setup
	Rscript tests/testthat/coverage.R

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
	R CMD check FeralCatEradication_0.2.0.tar.gz && \
	R CMD INSTALL FeralCatEradication_0.2.0.tar.gz
	
tests:
	R -e "testthat::test_dir('tests/testthat/', report = 'summary', stop_on_failure = TRUE)"
