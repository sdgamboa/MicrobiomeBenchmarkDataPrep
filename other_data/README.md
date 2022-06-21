

## Strammler

Get metadata files

```
## EBI metadata
date=$(date -u '+%Y%m%d_%H%M%S')
curl "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB11953&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,submitted_ftp,sra_ftp&format=tsv&download=true" > Stammler_2016_16S_spikein_ebi_metadata_$date.tsv
```
