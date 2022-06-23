
## This scripts contains code to generate sample metadata, count matrix,
## taxonomy table, and taxonomy tree of the Stammler 2016 dataset (using
## spike-in bacteria)

## Execute this scritp as an independent file (not in an R package session)

library(magrittr)
library(S4Vectors)
library(SummarizedExperiment)
library(TreeSummarizedExperiment)
library(dplyr)
library(tidyr)
library(purrr)
library(taxizedb)
library(readr)
library(tibble)

# Sample metadata ---------------------------------------------------------

## The sample metadata is obtained by combining the sample information available
## in EBI and the information provided by the article.

ebi_metadata_url <- "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB11953&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,submitted_ftp,sra_ftp&format=tsv&download=true"
article_metadata_url <- "https://static-content.springer.com/esm/art%3A10.1186%2Fs40168-016-0175-0/MediaObjects/40168_2016_175_MOESM8_ESM.txt"

ebi_metadata <- readr::read_tsv(
    ebi_metadata_url, show_col_types = FALSE, progress = FALSE
    ) %>%
    filter(grepl("ASCT.MID", submitted_ftp)) %>%
    mutate(
        sample_name = sub("^.+ASCT\\.(MID[0-9]+)_.+$", "\\1", submitted_ftp)
    ) %>%
    select(-sra_ftp) %>%
    relocate(sample_name)

article_metadata <- readr::read_tsv(
    article_metadata_url, show_col_types = FALSE, progress = FALSE
)

col_data <- left_join(
    article_metadata, ebi_metadata, by = c("SampleID" = "sample_name")
) %>%
    dplyr::rename(
        sample_id = SampleID,
        subject_id = Treatment,
        study_condition = Time,
        description = Description,
        ncbi_accession = run_accession
    ) %>%
    select(sample_id, subject_id, description, study_condition, ncbi_accession) %>%
    mutate(
        body_site = 'stool',
        country = 'Germany',
        sequencing_platform = 'Roche 454',
        pmid = '27329048'
    )


# Count matrix ------------------------------------------------------------

## Import data (count_matrix/otu_table and taxonomy)
data_url <- "https://static-content.springer.com/esm/art%3A10.1186%2Fs40168-016-0175-0/MediaObjects/40168_2016_175_MOESM10_ESM.txt"
data <- utils::read.table(
    file = data_url, header = TRUE, sep = "\t", row.names = 1,
    comment.char = "#", check.names = FALSE
    )

## Count matrix / OTU table
count_matrix <- as.matrix(data[,colnames(data) != "taxonomy"])


## This needs to be adjusted with respect to one of the bacteria.
s_ruber <- count_matrix['AF323500XXXX', ]
size_factor <- s_ruber/mean(s_ruber)

# RhizPos_h<-grep("AB247615XXXX",rownames(asct.OTUs))
#Search unique ID for A. acidiphilus:
# AliPos_h<-grep("AB076660XXXX",rownames(asct.OTUs))

SCML_data <- count_matrix
for(i in seq(ncol(count_matrix))){
    SCML_data[,i] <- round(SCML_data[,i] / size_factor[i])
}


SCML_data['AF323500XXXX',]

SCML_data['AB247615XXXX',]

SCML_data['AB076660XXXX',]

# Taxonomy table ----------------------------------------------------------

taxonomy <- tibble::tibble(
    taxon_name = rownames(data), taxonomy = data[["taxonomy"]]
)

# separate(
    # taxonomy, col = 'taxonomy', into = paste0('R', 1:8),
    # sep = '; __'
# ) %>% View()

# output <- vector("list", nrow(taxonomy))
# for (i in seq_along(output)) {
#      output[[i]] <- stringr::str_split(taxonomy[i, "taxonomy"], "; __")
#      output[[i]] <- unlist(c(taxonomy[i, 'TAXA', drop = TRUE], output[[i]]))
# }


# lengths <- unique(map_int(output, length))
# divide_by_length <- vector('list', length(lengths))
# for (i in seq_along(divide_by_length)) {
#     divide_by_length[[i]] <- keep(output, ~ length(.x) == lengths[i])
#     names(divide_by_length)[i] <- as.character(lengths[i])
# }

# dddf <- divide_by_length %>%
#     modify_depth(2, ~ as.data.frame(matrix(.x, nrow = 1))) %>%
#     map(~ bind_rows(.x))

# x = dddf[[1]]




# table(vapply(output, length, integer(1)))
#
# vapply(output, \(x) tail(x, 1), character(1))
#
# taxonomy[vapply(output, length, integer(1)) == 6, ] %>% View()
# as.data.frame(table(flatten_chr(output))) %>% View()


# Test that things work ---------------------------------------------------
rowTree <- NULL

colData <- col_data %>%
    tibble::column_to_rownames(var = "sample_id") %>%
    as.data.frame() %>%
    DataFrame()

rowData <- taxonomy %>%
    tibble::column_to_rownames(var = "taxon_name") %>%
    as.data.frame() %>%
    DataFrame()

tse <- TreeSummarizedExperiment(
    assays = SimpleList(abundance = SCML_data), ## calibrated counts
    colData = colData,
    rowData = rowData,
    rowTree = rowTree
)

# Export files ------------------------------------------------------------

## Add library size to sample metadata

col_data$library_size <- colSums(SCML_data)
col_data$sequencing_method <- '16S'
col_data$variable_region_16s <- 'V3-6'

SCML_data_df <- SCML_data %>%
    as.data.frame() %>%
    rownames_to_column(var = 'taxon_name') %>%
    as_tibble()

## Sample metadata
write_tsv(col_data, "data/Stammler_2016_16S_spikein_sample_metadata.tsv")
## Count matrix
write_tsv(SCML_data_df, 'data/Stammler_2016_16S_spikein_count_matrix.tsv')
## Taxonomy table
write_tsv(taxonomy, "data/Stammler_2016_16S_spikein_taxonomy_table.tsv")

