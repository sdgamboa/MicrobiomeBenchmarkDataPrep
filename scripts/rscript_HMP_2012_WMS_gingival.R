# Script to generate sample metadata, count matrix, taxonomy table, and
# taxonomy tree for the HMP WMS data

library(S4Vectors)
library(SummarizedExperiment)
library(TreeSummarizedExperiment)
library(dplyr)
library(tibble)
library(curatedMetagenomicData)
library(purrr)
library(magrittr)
library(readr)
library(ape)

subjects <- as_tibble(sampleMetadata) %>%
    filter(
        study_name == 'HMP_2012',
        body_subsite %in% c('subgingival_plaque', 'supragingival_plaque')
    ) %>%
    discard(~all(is.na(.x))) %>%
    split(factor(.$body_subsite)) %>%
    map(~ pull(.x, subject_id))

subjects_in_both <- intersect(subjects[[1]], subjects[[2]])


tse <- curatedMetagenomicData(
    'HMP_2012.relative_abundance', dryrun = FALSE, counts = TRUE,
    rownames = 'long'
)[[1]]

colData(tse)$library_size <- colSums(assay(tse))

body_subsite_lgl <-
    colData(tse)$body_subsite %in% c('subgingival_plaque', 'supragingival_plaque')

subject_id_lgl <-
    colData(tse)$subject_id %in% subjects_in_both

tse <- tse[,body_subsite_lgl & subject_id_lgl]
tse <- tse[rowSums(assay(tse)) > 0,]

# Separate elements -------------------------------------------------------

## sample metadata
sample_metadata <- tse %>%
    colData() %>%
    as.data.frame() %>%
    rownames_to_column(var = 'sample_id') %>%
    as_tibble() %>%
    discard(~ all(is.na(.x))) %>%
    select(
        -study_name
    ) %>%
    mutate(
        subject_id = sub('^.+_', '', subject_id)
    )

## Taxonomy table
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

taxonomy_table <- tse %>%
    rowData() %>%
    as.data.frame() %>%
    rownames_to_column(var = 'taxon_name') %>%
    rename(kingdom = superkingdom) %>%
    left_join(annotations, by = 'genus')


## count matrix
count_matrix <- tse %>%
    assay() %>%
    as.data.frame() %>%
    rownames_to_column(var = 'taxon_name') %>%
    as_tibble()


## phylogenetic tree

tree <- rowTree(tse)
tree$tip.label <- sub('^.+\\|[a-z]__', '', tree$tip.label)
rownames(tse) <- sub('^.+\\|[a-z]__', '', rownames(tse))

# Export files ------------------------------------------------------------

write_tsv(sample_metadata, "data/HMP_2012_WMS_gingival_sample_metadata.tsv")
write_tsv(count_matrix, 'data/HMP_2012_WMS_gingival_count_matrix.tsv')
write_tsv(taxonomy_table, 'data/HMP_2012_WMS_gingival_taxonomy_table.tsv')
write.tree(phy = tree, file = 'data/HMP_2012_WMS_gingival_taxonomy_tree.newick')
