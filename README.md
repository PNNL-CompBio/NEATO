NEATO: Network Enrichment Analysis Tool for Omics

Welcome to the NEATO home page. This tool is still heavily under development as we build out additional functionality. 

## Project Goal
The goal of this project is to facilitate the creation of networks to enhance functional enrichment for transcriptomics, proteomics, and other types of high throughput data. To achieve this goal, we plan to build a modular [R Shiny]() application that can run locally, on your own server, or be accessed via the [EMSL cloud]() depending on your use case. 


## Features

NEATO includes multiple features to enable end-to-end analysis as part of the MAP frameowrk.

### Protein/Gene filtering
We first enable filtering of the proteins/genes using p-values or log-fold change.

### Network integration
We allow for two types of network integration. The underlying network (as seen below in the architecture) is currently the STRING database. However, we plan to allow for additional networks in v2.

- Network-based enrichment: this type of enrichment maps genes/proteins that are identified as statistically significant in a biological experiment and maps them to the network. This approach is akin to what was done in the latest [STRING publication](). 
- OmicsIntegrator-based enrichment: this enrichment takes the network-based enrichment one step further to supplement the features identified by filtering step to identify proteins that might not have been captured.

### Network clustering

### Functional enrichment

## Architecture


![image](https://user-images.githubusercontent.com/65473513/171519485-dfddf6a5-8cfe-4f0d-bbfa-d5f7b55160ef.png)

The above is the planned architecture of the app. The app will take proteomics data (that had been processed by pmart) provided by the user as input. This data is used in combination with the interactions stored in the the STRING database, which will be housed locally, to contruct a network using one of the provided package. So far the app contains STRING and PCSF, but more could be added later on. These packages will output an iGraph object, which can be clustered using functions built into the iGraph package. From here, we can visualize the graph using the visNetwork package. VisNetwork displays a custimizable and interactive graph that colors the nodes based on cluster and allows the user to move around the nodes. Finally, we perform an enrichment on the data using the leapR package, as well as a local enrichment database.

## How To Install

Copy this git repository and open ui.R, server.R, or global.R from the prototype_app directory in R studio

```devtools::install_github("Loglew12/shinyenrich")``` to install dependencies

Click "Run app" in the top right corner of R Studio window.
