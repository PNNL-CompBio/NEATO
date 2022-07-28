FROM rocker/shiny:latest

RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libxml2 \
    libglpk-dev \
    libmagick++-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libgit2-dev

# RUN R -e "install.packages('devtools')"
# RUN R -e "library(devtools)"

COPY . /app
WORKDIR /app

# RUN R -e "install.packages('usethis', dependencies=TRUE, repos='http://cran.rstudio.com/')"
# RUN R -e "install.packages('roxygen2', dependencies=TRUE, repos='http://cran.rstudio.com/')"
# RUN R -e "install.packages('rversions', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('devtools', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('BiocManager', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('org.Hs.eg.db')"
RUN R -e "BiocManager::install('topGO')"
RUN R -e "devtools::install_github('Loglew12/shinyenrich')"

# RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'),  dependencies = TRUE)"
# RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"
# RUN Rscript -e "install.packages('renv', INSTALL_opts = c('--no-lock'), repos = 'https://cran.rstudio.com', dependencies = TRUE)"
# RUN Rscript -e "options(renv.consent = TRUE)"
# RUN Rscript -e "renv::consent(provided = TRUE)"
# RUN Rscript -e "renv::restore()"

# # expose port
EXPOSE 8080

# # run app on container start
CMD ["R", "-e", "shiny::runApp('../app', host = '0.0.0.0', port = 8080)"]