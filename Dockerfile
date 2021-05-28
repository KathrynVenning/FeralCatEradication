FROM islasgeci/base
RUN R -e "install.packages(c('plotly'), repos='http://cran.rstudio.com')"