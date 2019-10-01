#from https://github.com/kelshmo/ampad-DiffExp/blob/master/gene_level_analysis/ROSMAP_geneLevel.Rmd#L117

#######################
# Get and Join/Bind Inputs #
#######################

rosmap_inputs <- tibble(object = c("counts", "clinical metadata", "clinical metadata - uncensored ages", "assay metadata", "sequencing metadata", "mapping"), 
                        id = c("syn8691134", "syn3191087","syn7116000", "syn4300313", "syn8698240", "syn3382527"), 
                        version = c("1", "9","1", "1", "2", "7"))

data <- get_data(rosmap_inputs)

# Filter for required mapping columns and metadata variables, keep distinct pairs
data$thefile[data$object == "mapping"][[1]] <- distinct(dplyr::select(data$thefile[data$object == "mapping"][[1]], projid, rnaseq_id))
data$thefile[data$object == "clinical metadata"][[1]] <-dplyr::select(data$thefile[data$object == "clinical metadata"][[1]], -age_first_ad_dx, -age_death, -age_at_visit_max)

# Join metadata files
metadata <- data$thefile[data$object == "assay metadata"][[1]] %>%
  dplyr::left_join(data$thefile[data$object == "sequencing metadata"][[1]], by = c("Sampleid" = "sample")) %>% 
  dplyr::select(-projid) %>% 
  dplyr::left_join(data$thefile[data$object == "mapping"][[1]], by = c("Sampleid" = "rnaseq_id")) %>% 
  dplyr::left_join(data$thefile[data$object == "clinical metadata"][[1]], by = c("projid")) %>% 
  dplyr::left_join(data$thefile[data$object == "clinical metadata - uncensored ages"][[1]], by = c("projid"))
  
# Pick higher quality RIN batch
metadata <- metadata %>% 
  dplyr::group_by(Sampleid) %>%
  dplyr::top_n(1, RINcontinuous)

#######################
# Data pre-processing #
#######################

# Identify variables with missing data
missingid <- metadata$Sampleid[is.na(metadata$cogdx)|
                                 is.na(metadata$braaksc)|
                                 is.na(metadata$ceradsc)|
                                 is.na(metadata$RINcontinuous)|
                                 is.na(metadata$pmi)|
                                 is.na(metadata$RnaSeqMetrics__INTRONIC_BASES)|
                                 is.na(metadata$age_death)]

# Remove empty variables and unprocessed samples
metadata <- metadata %>% 
  ungroup %>% 
  dplyr::filter(Sampleid %in% colnames(data$thefile[data$object == "counts"][[1]])) %>%
  dplyr::filter(!is.na(cogdx), 
                !is.na(braaksc), 
                !is.na(ceradsc),
                !is.na(RINcontinuous),
                !is.na(pmi),
                !is.na(RnaSeqMetrics__PCT_INTRONIC_BASES),
                !is.na(age_death)) 

# Harmonize case-control status
metadata$diagnosis <- "other"
metadata$diagnosis[metadata$cogdx == 1 & metadata$braaksc <= 3 & metadata$ceradsc >= 3] <- "control"
metadata$diagnosis[metadata$cogdx == 4 & metadata$braaksc <= 4 & metadata$ceradsc >= 2] <- "AD"

# Add sex variable
metadata$sex <- "female"
metadata$sex[metadata$msex == 1] <- "male"

# Add APOE4 genotype = 0,1,2
metadata$APOE4 <- 0
metadata$APOE4[metadata$apoe_genotype %in% c(24,34)] <- 1
metadata$APOE4[metadata$apoe_genotype %in% c(44)] <- 2

# Get square of RIN
metadata$RINcontinuous2 <- metadata$RINcontinuous^2

# Filter for variables to be used in downstream analysis as covariates
covariates <- c("Sampleid","projid","cogdx", "diagnosis", "APOE4", "Batch", "sex", "race", "spanish", "cogdx", "APOE4",
                "RINcontinuous", "RINcontinuous2", "age_death", "pmi", "educ", "AlignmentSummaryMetrics__PCT_PF_READS_ALIGNED", 
                "RnaSeqMetrics__PCT_CODING_BASES","RnaSeqMetrics__PCT_INTERGENIC_BASES", "RnaSeqMetrics__PCT_INTRONIC_BASES", 
                "RnaSeqMetrics__PCT_RIBOSOMAL_BASES")

metadata <- dplyr::select(metadata, covariates)

# Clean variable names 
colnames(metadata) <- gsub("RnaSeqMetrics__", "", colnames(metadata))
colnames(metadata) <- gsub("AlignmentSummaryMetrics__", "", colnames(metadata))

# Convert rownames of counts from tracking id to ensemble gene id


