## Dummy Service "Filter Data"
## Last Updated: 2022_10_13

# Install latest version of rocker image
FROM bioconductor/bioconductor_docker:latest

# Install general use packages
RUN apt-get update && apt-get install -y \
    sudo \
    libssl-dev \
    libcurl4-openssl-dev \
    libhiredis-dev \
    libzmq3-dev \
    libxt-dev \
    libicu-dev \ 
    vim python3-venv 


COPY . /

# Install RWorker
RUN R -e "options(renv.config.connect.timeout = 600)"
RUN Rscript -e "install.packages(c('devtools', 'stringr', 'rzmq', 'redux', 'R6', 'processx', 'jsonlite', 'glue', 'futile.logger', 'crayon', 'base64enc','stringi'), repos =  'https://cran.rstudio.com')"
RUN Rscript -e "devtools::install_github('amwhite/rworker', repos = 'https://cran.rstudio.com')"
RUN R -e "BiocManager::install('org.Hs.eg.db')"
RUN R -e "BiocManager::install('topGO')"
RUN R -e "devtools::install_github('EMSL-Computing/mapDataAccess')"
RUN R -e "devtools::install_github('Loglew12/shinyenrich')"

# Set up python virtual env and install dependencies
COPY requirements.txt .
RUN python3 -m venv venv
RUN venv/bin/pip install --upgrade pip --trusted-host=pypi.org --trusted-host=files.pythonhosted.org 
RUN venv/bin/pip install -r requirements.txt --trusted-host=pypi.org --trusted-host=files.pythonhosted.org  

# install mapDataAccess R package from gitlab
RUN Rscript -e "install.packages(c('reticulate', 'yaml'), repos =  'https://cran.rstudio.com')"

# Add redis 
RUN Rscript -e "reticulate::use_virtualenv('/venv'); reticulate::py_install('Redis')"

# Install mapDataAccess
RUN Rscript -e 'remotes::install_github("EMSL-Computing/mapDataAccess")'

# Copy over required folders

# COPY rworker.R /
# COPY cfg /cfg/ 
CMD ["Rscript", "rworker.R", "&"]

