FROM islasgeci/base
RUN Rscript -e "install.packages(c('lintr', 'plotly', 'testthat'), repos='http://cran.rstudio.com')"
