NEATO - network enrichment analysis tool for Omics

An R Shiny application

## Architecture

![image](https://user-images.githubusercontent.com/65473513/171519485-dfddf6a5-8cfe-4f0d-bbfa-d5f7b55160ef.png)

The above is the planned architecture of the app. The app will take proteomics data (that had been processed by pmart) provided by the user as input. This data is used in combination with the interactions stored in the the STRING database, which will be housed locally, to contruct a network using one of the provided package. So far the app contains STRING and PCSF, but more could be added later on. These packages will output an iGraph object, which can be clustered using functions built into the iGraph package. From here, we can visualize the graph using the visNetwork package. VisNetwork displays a custimizable and interactive graph that colors the nodes based on cluster and allows the user to move around the nodes. Finally, we perform an enrichment on the data using the leapR package, as well as a local enrichment database.

## How To Install Locally

To run this package locally requires the installation of [RStudio](http://rstudio.com). Once yuou 

```
require(devtools)
devtools::install_github("PNNL-compbio/NEATO")
``` 


Click "Run app" in the top right corner of R Studio window.

## How to run via Docker


