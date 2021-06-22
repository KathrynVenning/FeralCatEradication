.PHONY: results

define runLint
	R -e "library(lintr)" \
	  -e "lint_dir('src', linters = with_defaults(line_length_linter(100)))"
endef

clean:
	rm --force Rplots.pdf
	rm --force --recursive reports/figures

lint:
	$(runLint)
	$(runLint) | grep -e "\^" && exit 1 || exit 0

results: src/FeralCatEradication.R src/matrixOperators.r
	mkdir reports/figures/ --parents
	Rscript src/FeralCatEradication.R