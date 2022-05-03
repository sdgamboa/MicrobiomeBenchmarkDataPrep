#! /bin/bash

var=$(which esearch)

if [ -z "$var" ]
then
    echo NCBI\'s e-utiliites must be installed first.
	echo Check:  https://www.ncbi.nlm.nih.gov/books/NBK179288/
    exit
fi

## Download metadata for the dataset of Ravel et. al., 2011
esearch -db sra -query SRA022855 | efetch -format runinfo > SRA022855_run_metadata.csv

