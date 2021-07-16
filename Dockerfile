FROM islasgeci/base:0.7.0
COPY . /workdir
RUN Rscript -e "install.packages(c('plotly'), repos='http://cran.rstudio.com')"
RUN Rscript -e "install.packages(c('covr', 'devtools', 'lintr', 'roxygen2', 'styler', 'testthat'), repos='http://cran.rstudio.com')"

RUN R -e "devtools::document()" && \
	R CMD build . && \
    R CMD check FeralCatEradication_0.1.0.tar.gz && \
    R CMD INSTALL FeralCatEradication_0.1.0.tar.gz
