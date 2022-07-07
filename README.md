STRINGdb_exploration - contains r markdowns used for exploring different packages.

Prototype_app - Current Shiny app

## Architecture

![image](https://user-images.githubusercontent.com/65473513/171519485-dfddf6a5-8cfe-4f0d-bbfa-d5f7b55160ef.png)

The above is the planned architecture of the app. The app will take proteomics data (that had been processed by pmart) provided by the user as input. This data is used in combination with the interactions stored in the the STRING database, which will be housed locally, to contruct a network using one of the provided package. So far the app contains STRING and PCSF, but more could be added later on. These packages will output an iGraph object, which can be clustered using functions built into the iGraph package. From here, we can visualize the graph using the visNetwork package. VisNetwork displays a custimizable and interactive graph that colors the nodes based on cluster and allows the user to move around the nodes. Finally, we perform an enrichment on the data using the leapR package, as well as a local enrichment database.

## How To Install

Copy this git repository and open ui.R, server.R, or global.R from the prototype_app directory in R studio

```devtools::install_github("Loglew12/shinyenrich")``` to install dependencies

Click "Run app" in the top right corner of R Studio window.
