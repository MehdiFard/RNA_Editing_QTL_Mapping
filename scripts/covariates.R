#!/usr/bin/env Rscript
# Adapted from Pantry (https://github.com/PejLab/Pantry), Copyright (c) 2024 Daniel Munro, MIT License
# Rscript for making covariate file for tensorqtl

# Loading required libraries
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(impute))
suppressPackageStartupMessages(library(argparse))


# Initialize parser
parser <- ArgumentParser(description = "Process input arguments for the R script")

# Define arguments
parser$add_argument("--vcf", type = "character", required = TRUE, help = "Path to VCF file")
parser$add_argument("--bed", type = "character", required = TRUE, help = "Path to BED file")
parser$add_argument("--sex_age_info", type = "character", default = NULL, help = "Path to sex and age information file (optional, if providing, it has to have only ID, sex, age columns - order as mentioned)")
parser$add_argument("--out", type = "character", required = TRUE, help = "Path to output file")
parser$add_argument("--geno_pcs", type = "integer", required = TRUE, help = "Number of genotype PCs")
parser$add_argument("--pheno_pcs", type = "integer", required = TRUE, help = "Number of phenotype PCs")

# Parse arguments
args <- parser$parse_args()

# Assign arguments to variables
VCF_FILE <- args$vcf
BED_FILE <- args$bed
SEX_AGE_INFO <- args$sex_age_info
OUT_FILE <- args$out
N_GENO_PCS <- args$geno_pcs
N_PHENO_PCS <- args$pheno_pcs

# Defining functions
load_geno <- function(filename) {
    gt <- readGT(filename)
    geno <- apply(gt, 2, function(x) c("0|0" = 0, "0|1" = 1, "1|0" = 1, "1|1" = 2, "0/0" = 0, "0/1" = 1, "1/0" = 1, "1/1" = 2)[x])
    rownames(geno) <- rownames(gt)
    geno
}

get_PCs <- function(df, n_pcs) {
    if (sum(is.na(df)) > 0) {
        df <- impute.knn(df)$data # Expects samples as columns
    }
    df <- df[apply(df, 1, var) != 0, ]
    pca <- prcomp(t(df), center = TRUE, scale = TRUE)
    pcs <- round(pca$x[, 1:n_pcs], 6)
    pcs_df <- data.frame(ID = colnames(pcs), t(pcs))
    colnames(pcs_df) <- c("ID", rownames(pcs))  # Column names starting with digits get 'fixed' and must be changed back
    pcs_df
}

# Phenotype PCs
pheno <- read.delim(BED_FILE, check.names = FALSE, row.names = 4)[, -(1:3)]
if (ncol(pheno) < 2) stop("Computing covariate PCs requires more than 1 sample.")
pheno_pcs <- get_PCs(pheno, N_PHENO_PCS)
pheno_pcs$ID <- paste("pheno", pheno_pcs$ID, sep = "_")

# Genotype PCs
geno <- load_geno(VCF_FILE)
missing_geno <- setdiff(colnames(pheno), colnames(geno))
if (length(missing_geno) > 0) {
    stop(length(missing_geno), " phenotype sample(s) not found in the VCF (tensorQTL needs genotypes for all): ",
         paste(head(missing_geno, 10), collapse = ", "))
}
geno <- geno[, colnames(pheno)]
geno_pcs <- get_PCs(geno, N_GENO_PCS)
geno_pcs$ID <- paste("geno", geno_pcs$ID, sep = "_")

# Making covar file
stopifnot(identical(colnames(geno_pcs), colnames(pheno_pcs)))
covars <- rbind(geno_pcs, pheno_pcs)

# Adding sex and age info to covar file (if provided)
if(!is.null(SEX_AGE_INFO)){
    sex_age_df <- read.table(SEX_AGE_INFO, header = FALSE, stringsAsFactors = FALSE) %>%
    rename(ID = 1, sex = 2, age = 3) %>%  # Rename columns
    filter(ID %in% colnames(covars)[-1]) # Keep rows matching covars' column names

    # Transpose and prepare for rbind
    sex_age_covar <- sex_age_df %>%
    column_to_rownames("ID") %>% # Use ID column as rownames for transposing
    t() %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    mutate(ID = c("sex", "age")) %>%
    select(c("ID", colnames(covars)[-1])) 

# Stop if colnames aren't matched
stopifnot(identical(colnames(sex_age_covar), colnames(covars)))
covars <- rbind(covars, sex_age_covar)
}

# Save the covar file
write.table(covars, OUT_FILE, sep = "\t", quote = FALSE, row.names = FALSE)

### END