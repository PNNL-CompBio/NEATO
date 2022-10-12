FROM bioconductor/bioconductor_docker:latest

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
    libgit2-dev \
    apt-transport-https \
    ca-certificates \
    curl \
    software-properties-common

# ENV CONDA_DIR /opt/conda
# RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
# /bin/bash ~/miniconda.sh -b -p /opt/conda

COPY . /app
WORKDIR /app/spras

RUN INSTALL_PATH=~/anaconda \
    && wget http://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh \
    && bash Miniconda2-latest* -fbp $INSTALL_PATH
ENV PATH=/root/anaconda/bin:$PATH

# RUN apt install apt-transport-https ca-certificates curl software-properties-common
RUN curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
RUN add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu `lsb_release -cs` test"
RUN apt update
RUN apt-get install -y docker-ce

RUN conda env create -f /app/spras/environment.yml

RUN R -e "options(renv.config.connect.timeout = 600)"

# RUN R -e "install.packages('callr', dependencies=TRUE, repos='http://cran.rstudio.com/')"
# RUN R -e "install.packages('processx', dependencies=TRUE, repos='http://cran.rstudio.com/')"
# RUN R -e "install.packages('gert', dependencies=TRUE, repos='http://cran.rstudio.com/')"
# RUN R -e "install.packages('roxygen2', dependencies=TRUE, repos='http://cran.rstudio.com/')"
# RUN R -e "install.packages('sessioninfo', dependencies=TRUE, repos='http://cran.rstudio.com/')"
# RUN R -e "install.packages('htmlwidgets', dependencies=TRUE, repos='http://cran.rstudio.com/')"
# RUN R -e "install.packages('dplyr', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('devtools', dependencies=TRUE, repos='http://cran.rstudio.com/')"
# RUN R -e "install.packages('BiocManager', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('org.Hs.eg.db')"
RUN R -e "BiocManager::install('topGO')"
RUN R -e "devtools::install_github('Loglew12/shinyenrich')"
# RUN R -e "devtools::install_github('IOR-Bioinformatics/PCSF')"
# RUN R -e "install.packages('shinydashboard', dependencies=TRUE, repos='http://cran.rstudio.com/')"

# RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'),  dependencies = TRUE)"
# RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"
# RUN Rscript -e "install.packages('renv', INSTALL_opts = c('--no-lock'), repos = 'https://cran.rstudio.com', dependencies = TRUE)"
# RUN Rscript -e "options(renv.consent = TRUE)"
# RUN Rscript -e "renv::consent(provided = TRUE)"
# RUN Rscript -e "renv::restore()"

# # expose port
EXPOSE 8080

SHELL ["conda", "run", "-n", "spras", "/bin/bash", "-c"]

# # run app on container start
# CMD ["R", "-e", "shiny::runApp('../app', host = '0.0.0.0', port = 8080)"]
ENTRYPOINT ["conda", "run", "-n", "spras", "R", "-e", "shiny::runApp('/app', host = '0.0.0.0', port = 8080)"]