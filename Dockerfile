FROM islasgeci/base
RUN Rscript -e "install.packages(c('lintr', 'plotly'), repos='http://cran.rstudio.com')"