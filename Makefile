.PHONY: results

clean:
	rm --force Rplots.pdf
	rm --force --recursive reports/figures

results: FeralCatEradication.R matrixOperators.r
	mkdir reports/figures/ --parents
	Rscript FeralCatEradication.R