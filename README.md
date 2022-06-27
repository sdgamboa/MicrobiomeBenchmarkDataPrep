
# Prepare data for MicrobiomeBenchmarkData

## Files

+ Data files uploaded to Zenodo are in the `data` directory.
+ Other files that were necessary are in the `other_data` directory.
+ Scripts for generating the files in both `data` and `other_data` are stored in the `scripts` directory.

## Datasets components

| Dataset | count matrix | sample metadata | taxonomy table | taxa annotations | phylogenetic tree |
| ------- | ------------ | --------------- | -------------- | ---------------- | ----------------- |
| HMP_2012_16S_gingival_V13 | yes | yes | yes | yes | yes |
| HMP_2012_16S_gingival_V35 | yes | yes | yes | yes | yes |
| HMP_2012_16S_gingival_V35_subset | yes | yes | yes | yes | yes |
| HMP_2012_WMS_gingival | yes | yes | yes | yes | yes |
| Stammler_2016_16S_spikein | yes | yes | yes | | |
| Ravel_2011_16S_BV | yes | yes | yes | yes | |

## General instructions

1. Update output files (those in the data directory).
2. If sample metadata files where modified, run the sampleMetadata.R script as
well.
3. After commit and push, upload files to Zenodo.
4. Go to the MicrobiomeBenchmarkData repo and update sampleMetadata and sysdata.



# Links

+ Data on Zenodo: <<<Insert DOI here when created>>>
+ R package on Bioconductor: <<<Insert link here when created>>>
+ R package source code: https://github.com/waldronlab/MicrobiomeBenchmarkData
+ R package issues: https://github.com/waldronlab/MicrobiomeBenchmarkData/issues
+ Scripts for preparing the datasets: https://github.com/waldronlab/MicrobiomeBenchmarkDataPrep
+ Code for reproducibility of the analyses: https://waldronlab.io/MicrobiomeBenchmarkDataAnalyses
