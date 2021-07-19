#!/usr/bin/env Rscript
#
# 

library('devtools')
library('roxygen2')

package.skeleton("FeralCatEradication")

my_rpackages <- as.package("FeralCatEradication")
load_all(my_rpackages)
document(my_rpackages)