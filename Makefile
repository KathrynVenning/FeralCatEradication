.PHONY: results

clean:
	rm --force --recursive reports/figures
	
results: FeralCatEradication.R matrixOperators.r
	mkdir reports/figures/ --parents
	Rscript FeralCatEradication.R