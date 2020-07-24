# Get Access to Data: https://adknowledgeportal.synapse.org/#/DataAccess/Instructions
source("../processFunctions/functions.R")

############################
# Get and Join/Bind Inputs #
############################
# Set config file
Sys.setenv(R_CONFIG_ACTIVE = "mayo")

# Login to Synapse. Make a Synapse account: https://docs.synapse.org/articles/getting_started.html
synapser::synLogin()

# Get and join metadata #
alignment_metadata <- get_data(config::get("alignment metadata")$synID,
                                config::get("alignment metadata")$version)

# Import counts
# Remove unmapped, multimapped, noFeature and ambiguous counts
counts <- get_data(config::get("counts")$synID,
                   config::get("counts")$version)[-c(1:4),]

# Harmonize clinical and alignment metadata variables
metadata <- get_data(config::get("clinical metadata cbe")$synID,
                     config::get("clinical metadata cbe")$version) %>%
  tidyr::separate(SampleID, c('Donor_ID', 'Tissue')) %>%
  dplyr::mutate(ID = paste(Donor_ID, Tissue, sep = '_')) %>%
  dplyr::rename(Gender = Sex)

metadata_tcx <- get_data(config::get("clinical metadata tcx")$synID,
                         config::get("clinical metadata tcx")$version) %>%
  tidyr::separate(SampleID, c('Donor_ID', 'Tissue')) %>%
  dplyr::mutate(ID = paste(Donor_ID, Tissue, sep = '_')) %>%
  dplyr::rename(Flowcell = FLOWCELL)

metadata <- bind_rows(metadata, metadata_tcx)

# Identify removed samples due to QC or missing data
qc <- bind_rows(
  get_data(config::get("qc cbe")$synID,
           config::get("qc cbe")$version),
  get_data(config::get("qc tcx")$synID,
           config::get("qc tcx")$version)
)

# Remove samples that have not passed QC
metadata <- metadata %>%
  dplyr::filter(!(ID %in% qc$`Sample Name`))

metadata <- metadata %>%
  dplyr::inner_join(alignment_metadata, by = c("ID" = "sample"))

# Identify variables with missing data
missing_data <- metadata$ID[is.na(metadata$RnaSeqMetrics__PCT_INTRONIC_BASES) |
                                 is.na(metadata$RIN) |
                                 is.na(metadata$AgeAtDeath) |
                                 is.na(metadata$Source)]

missing_data <- metadata[metadata$ID %in% missing_data,] %>%
  dplyr::select(ID, Tissue, Diagnosis) %>%
  dplyr::rename(`Sample ID` = ID,
                Dx = Diagnosis) %>%
  dplyr::mutate(Reason = "Missing clinical metadata")

qc <- bind_rows(qc, missing_data)

# Remove samples missing metadata
metadata <- metadata %>%
  dplyr::filter(!(ID %in% qc$`Sample Name`))

# Metadata covariate pre-processing #

# Harmonize Tissue, AgeAtDeath, Diagnosis, sex
metadata <- metadata %>%
  dplyr::rename(tissue = Tissue,
                age_death = AgeAtDeath,
                diagnosis = Diagnosis,
                source = Source,
                sex = Gender,
                sampleid = ID,
                donorid = Donor_ID,
                pmi = PMI,
                flowcell = Flowcell,
                rin = RIN,
                pct_pf_reads_aligned = `AlignmentSummaryMetrics__PCT_PF_READS_ALIGNED`,
                pct_coding_bases = `RnaSeqMetrics__CODING_BASES`,
                pct_intergenic_bases = `RnaSeqMetrics__PCT_INTERGENIC_BASES`,
                pct_intronic_bases = `RnaSeqMetrics__INTRONIC_BASES`,
                pct_ribosomal_bases = `RnaSeqMetrics__PCT_RIBOSOMAL_BASES`) %>%
  dplyr::mutate(tissue = gsub("CER", "CBE", tissue),
                age_death = gsub("_or_above", "", age_death),
                diagnosis = case_when(
                  diagnosis == "Control" ~ "control",
                  diagnosis == "Pathologic Aging" ~ "other",
                  diagnosis == "PSP" ~ "other",
                  TRUE ~ diagnosis),
                tissue_diagnosis = paste(tissue, diagnosis, sep = "_"),
                source_tissue_diagnosis = paste(source, tissue_diagnosis, sep = "_"),
                APOE4 = NA,
                pmi = NA)

# Add APOE4 genotype = 0,1,2
metadata$APOE4 <- 0
metadata$APOE4[metadata$ApoE %in% c(24,34)] <- 1
metadata$APOE4[metadata$ApoE %in% c(44)] <- 2

# Add APOE4-tissue variable
metadata$tissue_APOE4 <- paste(metadata$tissue, metadata$APOE4, sep = "_")

# Compute square of RIN
metadata <- metadata %>%
  dplyr::mutate(`rin2` = rin^2)

# Predict missing PMI
metadata$pmi <- NA
PMI <- metadata[,c("pmi", "sampleid", "donorid")]

# TODO: PMI approximation
# Get expression (convert counts to cpm)
# expr <- voom(counts, design = NULL)$E
# ind = which(rowSums(expr>=1)/dim(expr)[2] >= 0.5)
# expr = expr[ind,]

# Filter for variables to be used in downstream analysis as covariates
covariates <- c("sampleid", "donorid", "source", "diagnosis", "tissue", "tissue_diagnosis", "source_tissue_diagnosis",
                "tissue_APOE4", "APOE4", "sex", "rin", "rin2", "age_death", "pmi", "flowcell",
                "pct_pf_reads_aligned", "pct_coding_bases", "pct_intergenic_bases", "pct_intronic_bases",
                "pct_ribosomal_bases")

metadata <- dplyr::select(metadata, covariates)

# Remove samples from counts matrix
to_parse <- c("feature", colnames(counts)[colnames(counts) %in% metadata$sampleid])
counts <- counts[,to_parse]

# TODO: appropriately parse all the synIDs and version to add to provenance with config::get()
