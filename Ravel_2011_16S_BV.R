
## This script prepares the dataset for one of the datasets for bacterial
## vaginosis (dataset 1)

## Some notes:
## Nugent scores:
## low 0 - 3
## medium 4 - 6
## high 7 - 10

library(dplyr)
library(readr)
library(SummarizedExperiment)
library(S4Vectors)
source('utils.R')

# Sample metadata ---------------------------------------------------------

## The sample metadata must be assembled from two different sources:
## 1) From the NCBI database (see the get_ncbi_sra_data.sh script).
## 2) From the supplementary data of the article.

## Data from the NCBI
ncbi_metadata <- read.csv('Ravel_2011_16S_BV_SRA022855_run_metadata.csv')
select_cols <- c(
    sample_id = 'SampleName',
    sequencing_platform = 'Platform', NCBI_accession = 'Run',
    number_bases = 'bases'

)
ncbi_metadata <- ncbi_metadata[, select_cols]
colnames(ncbi_metadata) <- names(select_cols)
ncbi_metadata$PMID <- '20534435'

## Data from the supplementary tables
## This data had to be downloaded manually from the PNAS site.
## For some reason curl and wget didn't work at this time (05/02/2022)

st04 <- readxl::read_xlsx(
    'Ravel_2011_16S_BV_st04.xlsx', sheet = 1, range = 'A3:IU397'
) %>%
    as.data.frame() %>%
    rename(sample_id = 'Subject ID')
pnas_metadata <- st04[,1:7]
new_col_names <- c(
    'sample_id', 'ethnicity', 'ph', 'nugent_score',
    'nugent_score_category', 'community_group', 'number_reads'
)
colnames(pnas_metadata) <- new_col_names

sample_metadata <- left_join(ncbi_metadata, pnas_metadata, by = 'sample_id')

# Abundance matrix --------------------------------------------------------

x <-  sample_metadata[, 'sample_id', drop = FALSE]
y <- left_join(x, st04, by = "sample_id")
relab <- as.matrix(y[,9:ncol(y)])
rownames(relab) <- x[[1]]
relab <- t(relab)
rownames(relab) <- sub("^L\\.", "Lactobacillus", rownames(relab))

# Taxonomy ----------------------------------------------------------------

taxa_names <- rownames(relab) %>%
    sub("(Incertae_sedis_5_[1-2])$", "\\1dontdelete", .) %>%
    sub("_[0-9]+$", "", .) %>%
    sub("_Incertae_Sedis", "", .) %>%
    sub("_j$", "", .) %>%
    sub("_c$", "", .) %>%
    sub("_genera_incertae_sedis", "", .) %>%
    sub("^TM7$", "Candidatus Saccharibacteria", .)

## Some taxa are Bacteria incertae sedis, so we have to add that manually
taxa_names[grepl("dontdelete", taxa_names)] <- "Bacteria incertae sedis"

taxonomy <- taxize::classification(taxa_names, db = "ncbi")
taxonomy_table <- taxize_classification_to_taxonomy_table(taxonomy)
taxonomy_table <- taxonomy_table[,-1]
rownames(taxonomy_table) <- rownames(relab)

# Microbial signatures ----------------------------------------------------

url <- 'https://docs.google.com/spreadsheets/d/e/2PACX-1vTo8BOwXFXyzqBrZYIHuOLnYPnCLPKFt3rkYG5DYEycyjv7sCbNKuVhgL1LzLaT7DjqYOnmb02I_gMv/pub?gid=0&single=true&output=tsv'

sig <- read.table(url, sep = '\t', header = TRUE)
sig$NCBI_ID <- as.character(taxize::get_uid(sig$Taxon_name))
sig$Rank <- unlist(taxize::tax_rank(sig$Taxon_name, db = 'ncbi'))

sig_taxonomy <- taxize::classification(sig$Taxon_name, db = 'ncbi')
sig_taxonomy_table <- taxize_classification_to_taxonomy_table(sig_taxonomy)

sig <- cbind(sig, sig_taxonomy_table)

sig2 <- sig %>%
    mutate(Taxon_name = sub('^(\\w+) .*$', '\\1', Taxon_name)) %>%
    select(Taxon_name, Attribute) %>%
    distinct()

signatures <-
    left_join(taxonomy_table, sig2, by = c('genus' = 'Taxon_name'))

write.table(signatures, 'signatures.tsv', sep = '\t', row.names = FALSE,
            col.names = TRUE)

# Check with SE -----------------------------------------------------------

## Check Summarized Experiment
## This is only to check that all samples and taxa are in the right order
se <- SummarizedExperiment(
    assays = SimpleList(abundance = relab),
    colData = DataFrame(sample_metadata),
    rowData = DataFrame(taxonomy_table)
)

# Exports -----------------------------------------------------------------

## Export files


## For consistency let's export abundance matrix as counts

counts <-
   t(apply(t(relab), 2, function(x) {
       round(x * sample_metadata$number_reads / 100)
       }
      )
    )

write.table(
    counts, 'Ravel_2011_16S_BV_count_matrix.tsv', sep = '\t',
    row.names = TRUE, col.names = TRUE
)

## Export sample_metadata
write.table(
    sample_metadata, 'Ravel_2011_16S_BV_sample_metadata.tsv', sep = '\t',
    row.names = TRUE, col.names = TRUE
)

## Export taxonomy_table
write.table(
    taxonomy_table, 'Ravel_2011_16S_BV_taxonomy_table.tsv', sep = '\t',
    row.names = TRUE, col.names = TRUE
)

