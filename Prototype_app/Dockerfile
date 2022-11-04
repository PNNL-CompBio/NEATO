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


COPY . /app

RUN python3 -m venv /venv
RUN /venv/bin/pip install --upgrade pip
COPY requirements.txt .
RUN /venv/bin/pip install -r requirements.txt

WORKDIR /app

# RUN INSTALL_PATH=~/anaconda \
#     && wget http://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh \
#     && bash Miniconda2-latest* -fbp $INSTALL_PATH
# ENV PATH=/root/anaconda/bin:$PATH

# RUN apt install apt-transport-https ca-certificates curl software-properties-common
# RUN curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
# RUN add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu `lsb_release -cs` test"
# RUN apt update
# RUN apt-get install -y docker-ce

# RUN conda env create -f /app/environment.yml

RUN R -e "options(renv.config.connect.timeout = 600)"
RUN R -e "install.packages('devtools', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('org.Hs.eg.db')"
RUN R -e "BiocManager::install('topGO')"
RUN R -e "devtools::install_github('EMSL-Computing/mapDataAccess')"
RUN R -e "devtools::install_github('Loglew12/shinyenrich')"

# # expose port
EXPOSE 8080

# SHELL ["conda", "run", "-n", "NEATO", "/bin/bash", "-c"]

# # run app on container start
CMD ["R", "-e", "shiny::runApp('../app', host = '0.0.0.0', port = 8080)"]
# ENTRYPOINT ["conda", "run", "-n", "NEATO", "R", "-e", "shiny::runApp('/app', host = '0.0.0.0', port = 8080)"]