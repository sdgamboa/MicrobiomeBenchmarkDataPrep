
## Script to generate the sample metadata, count matrix, taxonomy table, and
## taxonomy tree for the dataset of smokers vs never smokers from
## NYCHANES and Beghini 2019

library(biomformat)
library(TreeSummarizedExperiment)
library(SummarizedExperiment)
library(S4Vectors)
library(dplyr)
library(tibble)

## The code in this script is based on two functions of the
## waldronlab/nychanesmicrobiome package:
## + loadQuiimeData: https://github.com/waldronlab/nychanesmicrobiome/blob/master/R/loadQiimeData.R
## + annotateFactors: https://github.com/waldronlab/nychanesmicrobiome/blob/master/R/annotateFactors.R

# Sample metadata ---------------------------------------------------------

## Relationship between nychanes ids and burklab ids
original_map <- readr::read_tsv(
    "https://raw.githubusercontent.com/waldronlab/nychanesmicrobiome/master/inst/extdata/original_map.tsv"
)

## Samples selected with nychanes ids and classified according to smoking status
sample_selection <- readr::read_tsv(
    "https://raw.githubusercontent.com/waldronlab/nychanesmicrobiome/master/inst/extdata/smokingsampleselection.tsv"
)

sample_selection$smokingstatus <- factor(
    x = as.matrix(sample_selection[,-1]) %*% 1:5,
    levels = 1:5,
    labels = c("alternativeonly","never","former","secondhand","cigarette")
)
colnames(sample_selection) <- toupper(colnames(sample_selection))

## Full metadata from nychanes
metadata_file <- tempfile()
download.file(
    url = "https://github.com/waldronlab/nychanesmicrobiome/blob/master/inst/extdata/public_v2_010518.sas7bdat?raw=true",
    destfile = metadata_file
)
metadata <- sas7bdat::read.sas7bdat(metadata_file) %>%
    as_tibble()

## Create sample metadata
sample_metadata <- original_map %>%
    left_join(metadata) %>%
    left_join(sample_selection) %>%
    relocate(SAMPLE_NAME = Burklab_ID) %>%
    mutate(
        SMOKER3CAT = factor(SMOKER3CAT, labels = c("Never smoker", "Current smoker", "Former smoker")),
        COTININE = as.numeric(COTININE),
        COTININECAT = factor(COTININECAT, labels = c("Active smoker", "Non-smoker + secondhand smoke", "Non-smoker no secondhand smoke", "Skipped")),
        RACE = factor(RACE, labels = c("Non-Hispanic White", "Non-Hispanic Black", "Hispanic", "Asian", "Other")),
        GENDER = factor(GENDER, labels = c("Male","Female")),
        AGEGRP5C = factor(AGEGRP5C, labels = c("20-29", "30-39", "40-49", "50-59", "60 AND OVER")),
        SR_ACTIVE = factor(SR_ACTIVE, labels = c("Very active","Somewhat active","Not very active/not active at all")),
        EDU4CAT = relevel(factor(EDU4CAT, labels = c( "Less than High school diploma", "High school graduate/GED", "Some College or associate's degree", "College graduate or more")), "College graduate or more"),
        OHQ_1 = factor(OHQ_1, levels = 1:7,labels = c("6 mos or less",">6mos, <= 1 yr",">1 yr, <= 2 yr",">2 yr, <= 3 yr",">3 yr, <=5 yr",">5 yr", "never")),
        OHQ_2 = factor(OHQ_2, levels = 1:6,labels = c("More than once a day", "Once a day", "Every few days", "Every few weeks", "Never", "No teeth or dentures")),
        OHQ_3 = factor(OHQ_3, labels = c("Yes","No",NA)),
        OHQ_5 = as.numeric(OHQ_5),
        OHQ_5_3CAT = cut(OHQ_5, breaks = c(-1,0,5,7), labels = c("0","1-5","6-7"), include.lowest = FALSE),
        DBTS_NEW = factor(DBTS_NEW, labels = c("Yes", "No",NA)),
        SMOKINGSTATUS = factor(SMOKINGSTATUS, levels = c("cigarette","never","former","alternativeonly","secondhand"),labels = c("Cigarette","Never smoker","Former smoker","Alternative smoker","Secondhand")),
        CURRENT_SMKER = factor(CURRENT_SMKER, labels = c("Yes","No")),
        CIGARETTES = factor(CIGARETTES, labels = c("Yes","No",NA)),
        BMI = as.double(BMI),
        SPAGE = as.integer(SPAGE),
        MERCURYU = as.double(MERCURYU ),
        SMOKER = factor(SMOKER, labels = c("Smoker", "Non smoker")),
        SMQ_13_1_1 = factor(SMQ_13_1_1, labels = c("Cigarettes", "Cigars/Cigarillos", "Hookah Pipe", "E-Cigsarettes", "NA")),
        # SMQ_13_1_1 = factor(SMQ_13_1_1, labels = c("Cigarettes", "Cigars/Cigarillos", #"Chewing tobacco", #"Snuff", "Hookah Pipe", "E-Cigsarettes", #"Nicotine patches/Gum/Other", "NA")),
        SMQ_7 = as.integer(SMQ_7),
        INC25KMOD = factor(INC25KMOD, labels = c("Less Than $20,000","$20,000-$49,999","$50,000-$74,999","$75,000-$99,999","$100,000 or More","NA")),
        POVGROUP6_0812CT = factor(POVGROUP6_0812CT, labels = c("0 to 5%", "5 to 10%", "10 to 20%", "20 to 30%", "30 to 40%", "40 to 100%")),
        SMQ_4UNIT = factor(SMQ_4UNIT, labels = c(1,4,52,NA)),
        US_BORN = factor(US_BORN, levels = c("1","2","NA"), labels = c("US-Born, 50 States, DC, PR and Territories","Other", NA)),
        DMQ_7YEAR = as.integer(DMQ_7YEAR),
        DMQ_2 = factor(DMQ_2,levels = 1:6,labels = c("Married","Widowed","Divorced","Separated","Never married","Living with partner")),
        AGEGRP3C = factor(AGEGRP3C, levels = 1:3, labels = c("20-34", "35-64", "65 and over")),
        EDU3CAT = relevel(factor(EDU3CAT, levels = 1:3, labels = c("High School Diploma or Less", "Some College or Associate's Degree", "College Graduate or More")),"College Graduate or More"),
        INC10K = as.numeric(as.character(INC10K)),
        INC3C = relevel(cut(INC10K, breaks = stats::quantile(INC10K, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),include.lowest = TRUE,labels = c("Less Than $30,000", "$30,000 - $60,000", "$60,000 or more")), "$60,000 or more"),
        A1C = as.integer(A1C),
        GLUCOSE = as.integer(GLUCOSE),
        DBQ_10_3CAT = cut(DBQ_10, breaks = c(-1,.999,5,999), labels = c("0-<1","1-5","6 or more"))
    )

## The following for loop is based on the function format_dbp in: https://github.com/waldronlab/nychanesmicrobiome/blob/master/R/annotateFactors.R
for (i in 3:10) {
    var <- as.numeric(sample_metadata[[paste0("DBQ_", i)]])
    varunit <- sample_metadata[[paste0("DBQ_", i, "UNIT")]]
    mult <- c(7, 1, 3/13)[varunit]
    mult[is.na(mult)] <- 0
    new_var <- var * mult
    sample_metadata[[paste0("DBQ_", i)]] <- new_var
}

## Retain metadata of only smokers and never smokers
col_data <- sample_metadata %>%
    filter(SMOKINGSTATUS %in% c("Cigarette", "Never smoker"))

# Taxonomy tree ----------------------------------------------------------------

row_tree_url <- "https://raw.githubusercontent.com/waldronlab/nychanesmicrobiome/master/inst/extdata/rep_set.tre"
row_tree_file <- tempfile()
download.file(row_tree_url, row_tree_file)
row_tree <- ape::read.tree(row_tree_file)

# Taxonomy table ----------------------------------------------------------

otu_table_file <- tempfile()
download.file(
    url = "https://github.com/waldronlab/nychanesmicrobiome/blob/master/inst/extdata/otu_table_mc10_w_tax.biom?raw=true",
    destfile = otu_table_file
)

## The otu_table_biom contains both the count matrix and OTU taxonomy table
otu_table_biom <- read_biom(otu_table_file)
row_data <- observation_metadata(otu_table_biom) %>%
    as.data.frame() %>%
    magrittr::set_colnames(
        c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    ) %>%
    as_tibble(rownames = "TAXA") %>%
    filter(
        if_any(everything(), ~ !grepl("Unassigned|Chloroplast|Mitochondria", .x))
    )

# OTU table ---------------------------------------------------------------
remove_cols <- c(
    ## controls
    "20151013CZTME1","20151020CZE3","20151020CZE4","20151020TME1","NC1","NC2",
    ## duplicates
    "NYDH0036","NYDH0051","NYDH0060","NYDH0152","NYDH0213","NYDH0487",
    "NYDH0492", "NYDH0522","NYDH0527","NYDH0545R","NYDH0649c", "NYDH0661",
    "NYDH0691", "NYDH0893","NYDH0931","NYDH0988","NYDH1042","NYDH1460",
    "NYDH1353"
)

count_matrix <- as(biom_data(otu_table_biom), "matrix")
count_matrix <- count_matrix[, !colnames(count_matrix) %in% remove_cols]
colnames(count_matrix) <- sub("c|R$", "", colnames(count_matrix))
count_matrix <- count_matrix[,colnames(count_matrix) %in% col_data$SAMPLE_NAME]
count_matrix <- count_matrix[rownames(count_matrix) %in% row_data$TAXA, ]
count_matrix <- count_matrix[rowSums(count_matrix) > 0, ]

# Intersects --------------------------------------------------------------
intersect_cols <- intersect(colnames(count_matrix), col_data$SAMPLE_NAME)
intersect_rows <- intersect(rownames(count_matrix), row_data$TAXA)

col_data <- col_data %>%
    filter(SAMPLE_NAME %in% intersect_cols)

row_data <- row_data %>%
    filter(TAXA %in% intersect_rows)

count_matrix <- count_matrix[row_data$TAXA, col_data$SAMPLE_NAME]

# Test if everything works ------------------------------------------------

colData <- col_data %>%
    column_to_rownames(var = "SAMPLE_NAME") %>%
    as.data.frame() %>%
    DataFrame()

rowData <- row_data %>%
    column_to_rownames(var = "TAXA") %>%
    as.data.frame() %>%
    DataFrame()

tse <- TreeSummarizedExperiment(
    assays = SimpleList(counts = count_matrix),
    colData = colData,
    rowData = rowData,
    rowTree = row_tree
)

# Export files ------------------------------------------------------------

## Export sample metadata
readr::write_tsv(
    x = col_data,
    file = "Beghini_2019_16S_smoking_sample_metadata.tsv"
)

## Export count matrix
write.table(
    x = count_matrix, sep = "\t", quote = TRUE, row.names = TRUE,
    col.names = TRUE,
    file = "Beghini_2019_16S_smoking_count_matrix.tsv"
)

## Export row data
readr::write_tsv(
    x = row_data,
    file = "Beghini_2019_16S_smoking_taxonomy_table.tsv"
)

## Export taxonomy tree
ape::write.tree(
    phy = row_tree,
    file = "Beghini_2019_16S_smoking_taxonomy_tree.newick"
)

