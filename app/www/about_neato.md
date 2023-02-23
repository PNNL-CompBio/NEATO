# NEATO: Network Enrichment Algorithm Tool for Omics

## What is NEATO

The goal of NEATO is to provide an easy-to-use application that allows users to 
perform functional enrichment analysis (FEA) on their proteomics data and visualize the
results using protein-protein interaction networks (PPIs). 

FEA provides additional insight into proteomic studies by identifying cellular pathways 
and functions most affected by a subset of proteins that are altered in a structure 
experiment or specific environmental study. In NEATO, the [LeapR](https://github.com/PNNL-CompBio/leapR) 
package is used to perform the functional enrichment, using the "enrichment_in_sets" 
method. The enrichment databases used for the app are found on the [enrichR website](https://maayanlab.cloud/Enrichr/#libraries)

PPI networks help visualize the results of FEA by using information from experiments 
that detect how proteins signal activity to one another. This helps researchers 
interpret how the differential expressed proteins could putatively be involved 
in similar cellular activity. In NEATO, the [STRING](https://string-db.org/) database is used to find the
interactions between proteins. Users have the option to view the interactions between
the top differentially expressed proteins, or can infer missing proteins, which NEATO
fills in using the [PCSF](https://github.com/IOR-Bioinformatics/PCSF) package.

## How to use NEATO

### Step 1:

**Upload the data.** NEATO accepts [Pmart](https://github.com/pmartR/pmartR) midpoint objects that contain proteomics 
data. It otherwise also accepts .csv and .tsv files containing protein identifiers,
log-fold changes, and a p-values for differencially expressed proteins. We also 
include test data based on data from [this publication by Blabby-Hass, et al.](https://www.sciencedirect.com/science/article/pii/S2211124722006076). 

### Step 2:

**Select data specific parameters.** NEATO requires the user to select the names
of the columns containing the protein identifiers, log-fold changes, and a p-values.
Users then select which of the supported species their data came from.

### Step 3:

**Select run specific parameters.** Users can select a p-value cutoff to filter 
out less significant data, change the number of proteins included in the run, and
adjust the score threshold for protein interaction given by [STRING](https://string-db.org/). Users can
also choose between gereating an explicit or inferred graph. The explicit graph shows
the interactions between proteins, while the inferred graph fills in missing proteins
using the [PCSF](https://github.com/IOR-Bioinformatics/PCSF) package.

### Step 4:

**Submit Data.** When finished setting parameters, users click the "Process Uploaded Data"
button to run the analysis. The "Check Status" button displays the status of the job
in a pop-up, which should be closed and reopenned to refresh. Once the job is done,
users can copy and paste the URL displayed in the status pop-up into a new web browser
tab, then move to NEATO's Output tab.

### Step 5:

**View Output.** Once on the output tab, users click the "Display Output" button
to load in and visualize NEATO's output. The top panel will display an interactive
PPI, with the graph clustered into sections based on an edge betweeness algroithm.
The bottom panel displays the functional enrichment table, which can be filtered 
by p-value and adjusted p-value. The background of the enrichment can be all proteins
uploaded by the user, or just those included in the network. Users can also click
a cluster to update the enrichment to include only the proteins in that cluster.

# Contact and Github

[Github](https://github.com/PNNL-CompBio/NEATO)

Sara Gosline - Sara.gosline@pnnl.gov \
Logan Lewis - logan.lewis@pnnl.gov



