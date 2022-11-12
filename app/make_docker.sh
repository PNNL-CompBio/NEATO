#!/bin/bash
#this is a file for quickly building all current docker containers for the NEATO app

# cd Enrichment_db
# docker build -t enrich_db .
# cd ../STRING_db/
# docker build -t inters_db .
#cd NEATO-Service
docker build -t neato_service NEATO-Service
#cd ..
docker build -t enrich_app .
docker-compose up
