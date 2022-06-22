
# This script is to generate the sample metadata, count matrix, taxonomy table,
# and taxonomy tree of the HMP_2012_16S_gingival_V35 and the
# HMP_2012_16S_gingival_V13 datasets.

library(SummarizedExperiment)
library(TreeSummarizedExperiment)
library(S4Vectors)
library(dplyr)
library(magrittr)
library(purrr)
library(tibble)
library(readr)
library(ape)

v35se <- HMP16SData::V35()
v13se <- HMP16SData::V13()

v35tse <- TreeSummarizedExperiment(
    assays = SimpleList(counts = assay(v35se)),
    colData = colData(v35se),
    rowData = rowData(v35se),
    rowTree = metadata(v35se)[["phylogeneticTree"]]
)

v13tse <- TreeSummarizedExperiment(
    assays = SimpleList(counts = assay(v13se)),
    colData = colData(v13se),
    rowData = rowData(v13se),
    rowTree = metadata(v13se)[["phylogeneticTree"]]
)

rm(v35se, v13se)
gc()

# Sample metadata ---------------------------------------------------------

## Sample metadata must be identical for both datasets

# v35_lgl <- colData(v35tse)$RUN_CENTER == "WUGC" & colData(v35tse)$VISITNO == 1
# v35_cols <- colnames(v35tse[, v35_lgl])
# v13_lgl <- colData(v13tse)$RUN_CENTER == "WUGC" & colData(v13tse)$VISITNO == 1
# v13_cols <- colnames(v13tse[, v13_lgl])
# intersect_samples <- intersect(v35_cols, v13_cols)

intersect_samples <- intersect(colnames(v13tse), colnames(v35tse))
select_samples <- c("Subgingival Plaque", "Supragingival Plaque")
v35_col_data <- colData(v35tse) %>%
    as.data.frame() %>%
    as_tibble(rownames = "sample_id") %>%
    filter(HMP_BODY_SUBSITE %in% select_samples,
           sample_id %in% intersect_samples,
           # RUN_CENTER == "WUGC",
           # one run center might be incorrectly annotated as "0"
           # RUN_CENTER != "0",
           !is.na(SRS_SAMPLE_ID)
        ) %>%
    group_by(HMP_BODY_SUBSITE, RSID) %>%
    slice_min(RSID) %>%
    set_colnames(tolower(colnames(.))) %>%
    rename(
        subject_id = rsid,
        visit_number = visitno,
        gender = sex,
        body_site = hmp_body_site,
        body_subsite = hmp_body_subsite,
        ncbi_accession = srs_sample_id
    ) %>%
    mutate(
        gender = tolower(gender),
        body_site = tolower(body_site),
        body_subsite = sub(" " , "_", tolower(body_subsite)),
        # sequencing_type = '16S',
        sequencing_platform = 'Roche 454',
        pmid = '22699609'
    )


v13_col_data <- colData(v13tse) %>%
    as.data.frame() %>%
    as_tibble(rownames = "sample_id") %>%
    filter(HMP_BODY_SUBSITE %in% select_samples,
           sample_id %in% intersect_samples,
           # RUN_CENTER == "WUGC",
           # one run center might be incorrectly annotated as "0"
           # RUN_CENTER != "0",
           !is.na(SRS_SAMPLE_ID)
        ) %>%
    group_by(HMP_BODY_SUBSITE, RSID) %>%
    slice_min(RSID) %>%
    set_colnames(tolower(colnames(.))) %>%
    rename(
        subject_id = rsid,
        visit_number = visitno,
        gender = sex,
        body_site = hmp_body_site,
        body_subsite = hmp_body_subsite,
        ncbi_accession = srs_sample_id
    ) %>%
    mutate(
        gender = tolower(gender),
        body_site = tolower(body_site),
        body_subsite = sub(" " , "_", tolower(body_subsite)),
        # sequencing_type = '16S',
        sequencing_platform = 'Roche 454',
        pmid = '22699609'
    )

if (nrow(v35_col_data) > nrow(v13_col_data)) {

    v35_col_data <-
        v35_col_data[v35_col_data$sample_id %in% v13_col_data$sample_id,]

} else if (nrow(v13_col_data) < nrow(v35_col_data)) {

    v13_col_data <-
        v13_col_data[v13_col_data$sample_id %in% v35_col_data$sample_id,]

}

v13_col_data[v13_col_data$sample_id == '700023388',]$run_center <- 'WUGC'

if (identical(v13_col_data, v35_col_data)) {
    hmp_gingival_samples <-
        intersect(v13_col_data$sample_id, v35_col_data$sample_id)
    # generate sample metadata
    col_data <- v13_col_data # Either v13 or v35 works
} else {
    stop("Metadata are not identical.")
}

## Add some other necessary columns
col_data <- col_data %>%
    mutate(
        study_condition = 'control',
        disease = 'healthy'
    )


## Remove taxon_name with no counts

v35_subset <- v35tse[,hmp_gingival_samples]
v35_subset <- v35_subset[rowSums(assay(v35_subset)) > 0, ]

v13_subset <- v13tse[,hmp_gingival_samples]
v13_subset <- v13_subset[rowSums(assay(v35_subset)) > 0, ]

rm(v35tse, v13tse)
gc()

# Count matrix ------------------------------------------------------------

v35_count_matrix <- assay(v35_subset)
v13_count_matrix <- assay(v13_subset)

# Taxonomy tree -----------------------------------------------------------

v35_row_tree <- rowTree(v35_subset)
v13_row_tree <- rowTree(v13_subset)

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

v35_row_data <- rowData(v35_subset) %>%
    as.data.frame() %>%
    as_tibble(rownames = "taxon_name") %>%
    set_colnames(tolower(colnames(.))) %>%
    select(-consensus_lineage) %>%
    left_join(annotations, by = 'genus')

v13_row_data <- rowData(v13_subset) %>%
    as.data.frame() %>%
    as_tibble(rownames = "taxon_name") %>%
    set_colnames(tolower(colnames(.))) %>%
    select(-consensus_lineage) %>%
    left_join(annotations, by = 'genus')

# Check that everything works ---------------------------------------------

colData <- col_data %>%
    tibble::column_to_rownames(var = "sample_id") %>%
    as.data.frame() %>%
    DataFrame()

v35rowData <- v35_row_data %>%
    tibble::column_to_rownames(var = "taxon_name") %>%
    as.data.frame() %>%
    DataFrame()

v13rowData <- v13_row_data %>%
    tibble::column_to_rownames(var = "taxon_name") %>%
    as.data.frame() %>%
    DataFrame()

tse <- TreeSummarizedExperiment(
    assays = SimpleList(counts = v13_count_matrix),
    colData = colData,
    rowData  = v13rowData,
    rowTree = v13_row_tree
)


# Export files ------------------------------------------------------------

## Add column with taxon_name as header

v35_count_matrix_df <- v35_count_matrix %>%
    as.data.frame() %>%
    rownames_to_column(var = 'taxon_name') %>%
    as_tibble()

v13_count_matrix_df <- v13_count_matrix %>%
    as.data.frame() %>%
    rownames_to_column(var = 'taxon_name') %>%
    as_tibble()

## Save sample metadata
write_tsv(col_data, "data/HMP_2012_16S_gingival_sample_metadata.tsv")
## Save count matrix
write_tsv(v35_count_matrix_df, "data/HMP_2012_16S_gingival_V35_count_matrix.tsv")
write_tsv(v13_count_matrix_df, "data/HMP_2012_16S_gingival_V13_count_matrix.tsv")
## Save phylogenetic trees
write.tree(phy = v35_row_tree, file = "data/HMP_2012_16S_gingival_V35_taxonomy_tree.newick")
write.tree(phy = v13_row_tree, file = "data/HMP_2012_16S_gingival_V13_taxonomy_tree.newick")
## Save taxonomy tables
write_tsv(v35_row_data, "data/HMP_2012_16S_gingival_V35_taxonomy_table.tsv")
write_tsv(v13_row_data, 'data/HMP_2012_16S_gingival_V13_taxonomy_table.tsv')

