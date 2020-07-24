#from https://github.com/kelshmo/ampad-DiffExp/blob/master/gene_level_analysis/ROSMAP_geneLevel.Rmd#L117

############################
# Get and Join/Bind Inputs #
############################
Sys.setenv(R_CONFIG_ACTIVE = "rosmap")
# Join synIDs and select sample identifiers to map sequencing metadata
alignment_metadata <- get_data(config::get("alignment metadata")$synID,
                                config::get("alignment metadata")$version)

counts <- get_data(config::get("counts")$synID)[-c(1:4),]

# Filter for required mapping columns and metadata variables, keep distinct pairs
mapping <- get_data(config::get("mapping")$synID, config::get("mapping")$version) %>%
  dplyr::select(., projid, rnaseq_id) %>%
  distinct()
clinical_metadata <- get_data(config::get("clinical metadata")$synID, config::get("clinical metadata")$version) %>%
  dplyr::select(-age_first_ad_dx, -age_death, -age_at_visit_max)

# Join metadata files
metadata <- get_data(config::get("assay metadata")$synID,
                     config::get("assay metadata")$version) %>%
  dplyr::left_join(alignment_metadata,
                   by = c("Sampleid" = "specimenID")) %>%
  dplyr::select(-projid) %>%
  dplyr::left_join(mapping,
                   by = c("Sampleid" = "rnaseq_id")) %>%
  dplyr::left_join(clinical_metadata,
                   by = c("projid")) %>%
  dplyr::left_join(get_data(config::get("protected clinical metadata")$synID,
                            config::get("protected clinical metadata")$version),
                   by = c("projid"))

#######################
# Data pre-processing #
#######################

# Identify variables with missing data
removed_samples <- metadata$Sampleid[is.na(metadata$cogdx) |
                                 is.na(metadata$braaksc) |
                                 is.na(metadata$ceradsc) |
                                 is.na(metadata$RINcontinuous) |
                                 is.na(metadata$pmi) |
                                 is.na(metadata$RnaSeqMetrics__INTRONIC_BASES) |
                                 is.na(metadata$age_death)]

# Pick higher quality RIN batch
metadata <- metadata %>%
  dplyr::group_by(Sampleid) %>%
  dplyr::top_n(1, RINcontinuous)

# Remove variables with missing data
metadata <- metadata %>%
  ungroup %>%
  dplyr::filter(Sampleid %in% colnames(counts)) %>%
  dplyr::filter(!(Sampleid %in% removed_samples))

# Remove samples that do not meet criteria
counts <- counts[,metadata$Sampleid]


# Harmonize case-control status
metadata$diagnosis <- "other"
metadata$diagnosis[metadata$cogdx == 1 & metadata$braaksc <= 3 & metadata$ceradsc >= 3] <- "control"
metadata$diagnosis[metadata$cogdx == 4 & metadata$braaksc <= 4 & metadata$ceradsc >= 2] <- "AD"

# Add sex variable
metadata$sex <- "female"
metadata$sex[metadata$msex == 1] <- "male"

# Add tissue variable
metadata$tissue <- "DLPFC"

# Add APOE4 genotype = 0,1,2
metadata$APOE4 <- 0
metadata$APOE4[metadata$apoe_genotype %in% c(24,34)] <- 1
metadata$APOE4[metadata$apoe_genotype %in% c(44)] <- 2

# Compute square of RIN
metadata$RIN2 <- metadata$RINcontinuous^2

# Filter for variables to be used in downstream analysis as covariates
covariates <- c("Sampleid","projid","cogdx", "diagnosis", "APOE4", "Batch", "sex", "race", "spanish", "cogdx", "APOE4",
                "RINcontinuous", "RIN2", "age_death", "pmi", "educ", "AlignmentSummaryMetrics__PCT_PF_READS_ALIGNED",
                "RnaSeqMetrics__PCT_CODING_BASES","RnaSeqMetrics__PCT_INTERGENIC_BASES", "RnaSeqMetrics__PCT_INTRONIC_BASES",
                "RnaSeqMetrics__PCT_RIBOSOMAL_BASES")

metadata <- dplyr::select(metadata, covariates)

# Clean variable names
colnames(metadata) <- gsub("RnaSeqMetrics__", "", colnames(metadata))
colnames(metadata) <- gsub("AlignmentSummaryMetrics__", "", colnames(metadata))
colnames(metadata) <- tolower(colnames(metadata))

# Convert rownames of counts to ensemble gene id
# tmp = data.frame(Gene.ID = rownames(counts)) %>%
#   dplyr::mutate(ID = Gene.ID) %>%
#   tidyr::separate(ID, c('ensembl_gene_id', 'position'), sep = '\\.')
# rownames(tmp) = tmp$ensembl_gene_id
# rownames(counts) = tmp[rownames(counts), 'ensembl_gene_id']
