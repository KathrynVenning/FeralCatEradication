.PHONY: results

clean:
	rm --force Rplots.pdf
	rm --force --recursive reports/figures

results: src/FeralCatEradication.R src/matrixOperators.r
	mkdir reports/figures/ --parents
	Rscript src/FeralCatEradication.R