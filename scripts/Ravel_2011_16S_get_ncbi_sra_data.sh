#! /bin/bash

## This script was used to download the metadata for the Ravel_2011_16S
## dataset from the NCBI

## NCBI e-utils is required. Check:
## https://www.ncbi.nlm.nih.gov/home/tools/
## https://www.ncbi.nlm.nih.gov/books/NBK179288/
var=$(which esearch)

if [ -z "$var" ]
then
    echo NCBI\'s e-utiliites must be installed first.
	echo Check:  https://www.ncbi.nlm.nih.gov/books/NBK179288/
    exit
fi

## Download metadata for the dataset of Ravel et. al., 2011
esearch -db sra -query SRA022855 | efetch -format runinfo > SRA022855_run_metadata.csv

