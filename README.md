# NEATO: Network Enrichment Analysis Tool for Omics

Welcome to the NEATO home page. This tool is still heavily under development as we build out additional functionality. 

## Project Goals
The primary goal of this project is to facilitate the creation of networks to enhance functional enrichment for transcriptomics, proteomics, and other types of high throughput data. To achieve this goal, we plan to build a modular [R Shiny](http://rstudio.com) application that can run locally, on your own server, or be accessed via the [EMSL cloud]() depending on your use case. We also hope to expand the use of network analysis tools and make them broadly interpretable.

## NEATO Features

NEATO includes multiple features to enable end-to-end analysis as part of the MAP frameowrk.

### Protein/Gene filtering
We first enable filtering of the proteins/genes using p-values or log-fold change. For proteins, you can use the output of pMART to select the proteins of interest. 

### Network integration
We allow for two types of network integration. The underlying network (as seen below in the architecture) is currently the STRING database. However, we plan to allow for additional networks in v2.

- Primary network: this type of enrichment maps genes/proteins that are identified as statistically significant in a biological experiment and maps them to the network. This approach is akin to what was done in the latest [STRING publication](). 
- Augmented network: this enrichment takes the network-based enrichment one step further to supplement the features identified by filtering step to identify proteins that might not have been captured.

### Network clustering
Once we have mapped the proteins or genes to the network (or augmented network) we employ Louvain clustering to group the proteins based on the edge connectivity. This enables us to do functional enrichment on subgraphs of the network in addition to the entire graph.

### Functional enrichment
Currently we searchf for GO biological pathways in which the selected proteins are more abundant in the pathway than we'd expect by chance. 

## Architecture

The NEATO architecture is designed so that we can configure/update the following modules:
- Network: currently we use the STRING network but can also use a kinase-substrate network or other type of interaction network
- Functional enrichment mapping: we hope to augment the enrichment to allow for different species and different types of pathways
- Network algorithm: we hope to add in additional network functionality

![image](https://user-images.githubusercontent.com/65473513/221019667-0a5a0b20-3299-4902-8e22-27471a6edfaa.png)

The above is the planned architecture of the app. The app will take proteomics data (that had been processed by pmart) provided by the user as input. Once uploaded the data is stored on MinIO, which ensures the data persists after users close the application, and also the data to be passed to the service container. This data is used in combination with the interactions stored in the STRING database, which will be housed locally, to contruct a network using one of the provided algorithms hosted in the NEATO service contianer. So far the app contains STRING and PCSF, but more could be added later on. The Redis package is used to send the request job from NEATO's frontend to the backend service container, and delivers status updates from the backend to the frontend for the user. These packages will output an iGraph object, which can be clustered using functions built into the iGraph package. From here, the output is saved to MinIO and loaded back into the frontend so we can visualize the graph using the visNetwork package. VisNetwork displays a custimizable and interactive graph that colors the nodes based on cluster and allows the user to move around the nodes. Finally, we perform an enrichment on the data using the leapR package, as well as a local enrichment database.

## How To Run

The NEATO package is not currently designed to be ran locally, and instead is being prepared for deployment on the [Multiomics Applicatoin Portal (MAP)](https://map.emsl.pnnl.gov/app/map).


