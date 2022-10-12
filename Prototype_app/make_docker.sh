#!/bin/bash
#this is a file for quickly building all current docker containers for the NEATO app

cd Enrichment_db
docker build -t enrich_db .
cd ../STRING_db/
docker build -t inters_db .
cd ..
docker build -t enrich_app .
# cd /Users/lewi052/enrich_proj/MAP/Prototype_app/spras/docker-wrappers/OmicsIntegrator2
# docker build -t reedcompbio/omics-integrator-2 -f Dockerfile .
docker-compose up