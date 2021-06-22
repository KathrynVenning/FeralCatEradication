FROM islasgeci/base
RUN Rscript -e "install.packages(c('plotly'), repos='http://cran.rstudio.com')"