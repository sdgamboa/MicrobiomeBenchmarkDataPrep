
## Script to generate the sample metadata, count matrix, taxonomy table, and
## taxonomy tree of the subset of the HMP 2012 data used in Calgaro 2020.

library(SummarizedExperiment)
library(TreeSummarizedExperiment)
library(dplyr)
library(tibble)
library(readr)
library(ape)
library(magrittr)
library(HMP16SData)

v35se <- V35()
colData(v35se)$library_size <- colSums(assay(v35se))

import_calgaro_2020 <- function() {
    load(url("https://github.com/mcalgaro93/sc2meta/blob/master/data/16Sdatasets_for_replicability_filtered.RData?raw=true"))
    mia::makeTreeSummarizedExperimentFromPhyloseq(
        ps_list_16S[["Subgingival_Supragingival"]]
    )
}

tse <- import_calgaro_2020()

col_names <- colnames(tse)

colData(tse)$library_size <- colData(v35se[,col_names])$library_size


# colData(tse)$library_size <- colSums(assay(tse))

# Sample metadata ---------------------------------------------------------

col_data <- colData(tse) %>%
    as.data.frame() %>%
    rownames_to_column(var = 'sample_id') %>%
    as_tibble() %>%
    set_colnames(tolower(colnames(.))) %>%
    dplyr::rename(
        subject_id = rsid,
        visit_number = visitno,
        gender = sex,
        body_site = hmp_body_site,
        body_subiste = hmp_body_subsite,
        ncbi_accession = srs_sample_id
    ) %>%
    mutate(
        sequencing_method = '16S',
        variable_region_16s = 'V3-5'
    )


# Count matrix ------------------------------------------------------------

count_matrix <- assay(tse)

# Taxonomy tree -----------------------------------------------------------

row_tree <- rowTree(tse)

# Taxonomy table ----------------------------------------------------------

## Add biosis information
fname <- "https://raw.githubusercontent.com/waldronlab/nychanesmicrobiome/master/inst/extdata/biosis.tsv"
annotations <- read.table(fname, header = TRUE, sep = '\t') %>%
    set_colnames(c('genus', 'taxon_annotation')) %>%
    mutate(
        taxon_annotation = case_when(
            taxon_annotation == 'F Anaerobic' ~ 'facultative_anaerobic',
            taxon_annotation == 'Aerobic' ~ 'aerobic',
            taxon_annotation == 'Anaerobic' ~ 'anaerobic'
        )
    )

row_data <-
    rowData(tse) %>%
    as.data.frame() %>%
    rownames_to_column(var = 'taxon_name') %>%
    as_tibble() %>%
    set_colnames(tolower(colnames(.))) %>%
    dplyr::rename(kingdom = superkingdom) %>%
    left_join(annotations, by = "genus")

# Export files ------------------------------------------------------------

count_matrix_df <- count_matrix %>%
    as.data.frame() %>%
    rownames_to_column(var = 'taxon_name') %>%
    as_tibble()

## Export sample metadata
write_tsv(col_data, 'data/HMP_2012_16S_gingival_V35_subset_sample_metadata.tsv')
## Export count matrix
write_tsv(count_matrix_df, 'data/HMP_2012_16S_gingival_V35_subset_count_matrix.tsv')
## Export taxonomy table
write_tsv(row_data, 'data/HMP_2012_16S_gingival_V35_subset_taxonomy_table.tsv')
## Export phylogenetic tree
write.tree(row_tree, 'data/HMP_2012_16S_gingival_V35_subset_taxonomy_tree.newick')

