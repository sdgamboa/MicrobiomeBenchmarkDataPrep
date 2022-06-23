library(dplyr)
library(readr)
library(purrr)

format_sample_metadata <- function(df) {

    ## One of the main formats is that everything is changed to lower case,
    ## althoug this might not be what is required needed for all cases.

    isSampleNameAndCharacter <- function(x) {
        purrr::map_lgl(x, is.character) & !grepl("(SAMPLE_NAME)|(sample_name)", colnames(x))
    }

    vct_lgl <- isSampleNameAndCharacter(df)

    df %>%
        magrittr::set_colnames(., tolower(colnames(.))) %>%
        magrittr::set_colnames(., gsub(" ", "_", colnames(.))) %>%
        purrr::map_if(.x = ., .p = vct_lgl, .f = tolower) %>%
        purrr::map_if(.x = ., .p = is.character, .f = ~gsub(" ", "_", .x)) %>%
        purrr::map_at(.at = "sample_name", .f = ~ as.character(.)) %>%
        tibble::as_tibble()
}

HMP_2012_16S_gingival_V13 <-
    read_tsv("data/HMP_2012_16S_gingival_V13_sample_metadata.tsv", col_types = cols(sequencing_method = col_character()))

HMP_2012_16S_gingival_V35 <-
    read_tsv("data/HMP_2012_16S_gingival_V35_sample_metadata.tsv", col_types = cols(sequencing_method = col_character()))

HMP_2012_16S_gingival_V35_subset <-
    read_tsv("data/HMP_2012_16S_gingival_V35_subset_sample_metadata.tsv", col_types = cols(sequencing_method = col_character()))

HMP_2012_WMS_gingival <-
    read_tsv("data/HMP_2012_WMS_gingival_sample_metadata.tsv", col_types = cols(sequencing_method = col_character()))

# Beghini_2019_16S_smoking <- readr::read_tsv(
    # file = "data-raw/Beghini_2019_16S_smoking_sample_metadata.tsv"
# ) %>%
    # format_sample_metadata()

Stammler_2016_16S_spikein <-
    read_tsv("data/Stammler_2016_16S_spikein_sample_metadata.tsv", col_types = cols(sequencing_method = col_character()))

Ravel_2011_16S_BV <-
    read_tsv("data/Ravel_2011_16S_BV_sample_metadata.tsv", col_types = cols(sequencing_method = col_character()))

# Export merged metadata --------------------------------------------------

datasets <- list(
    HMP_2012_16S_gingival_V13 = HMP_2012_16S_gingival_V13,
    HMP_2012_16S_gingival_V35 = HMP_2012_16S_gingival_V35,
    HMP_2012_16S_gingival_V35_subset = HMP_2012_16S_gingival_V35_subset,
    HMP_2012_WMS_gingival = HMP_2012_WMS_gingival,
    # Beghini_2019_16S_smoking = Beghini_2019_16S_smoking,
    Stammler_2016_16S_spikein = Stammler_2016_16S_spikein,
    Ravel_2011_16S_BV = Ravel_2011_16S_BV
)

sampleMetadata <- datasets %>%
    map(~ mutate(.x, sample_id = as.character(sample_id))) %>%
    map_if(.p = ~ 'subject_id' %in% colnames(.x),
           .f = ~ mutate(.x, subject_id = as.character(subject_id))
    ) %>%
    bind_rows(.id = 'dataset') %>%
    mutate(gender = tolower(gender))

## Export file as tsv (this for upload to Zenodo)
write_tsv(sampleMetadata, "data/sampleMetadata.tsv")

## Save data
# usethis::use_data(sampleMetadata, overwrite = TRUE)

