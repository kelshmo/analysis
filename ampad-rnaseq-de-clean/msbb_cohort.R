# Get Access to Data: https://adknowledgeportal.synapse.org/#/DataAccess/Instructions
source("../processFunctions/functions.R")

# Set config file
Sys.setenv(R_CONFIG_ACTIVE = "mssm")

# Login to Synapse. Make a Synapse account: https://docs.synapse.org/articles/getting_started.html
synapser::synLogin()

# Get and join metadata
# RNA-seq metadata
#
assay_metadata <- get_data(config::get("assay metadata")$synID, config::get("assay metadata")$version)

#Get clarity on why individual identifier differs
foo <- assay_metadata %>%
  dplyr::mutate(ext = ifelse(grepl("bam", fileName), "bam", NA),
         ext = ifelse(grepl("fastq", fileName), "fastq", ext)) %>%
  dplyr::filter(ext == "bam") %>%
  dplyr::select(individualIdentifier, individualIdentifier.inferred)
  #dplyr::select(sampleIdentifier, BrodmannArea, barcode, individualIdentifier, batch, RIN)

sequencing_metadata <- get_data(config::get("sequencing metadata")$synID)
# Join synIDs and select sample identifiers to map sequencing metadata
# sequencing_metadata <- get_data(config::get("sequencing metadata")$synID,
#                                 config::get("sequencing metadata")$version) %>%
#   full_join(get_data(config::get("synID mapping")$synID), by = c("sample" = "id")) %>%
#   dplyr::select(specimenID, everything(), -sample, -versionNumber) %>%
#   distinct()

# Import counts
# Remove unmapped, multimapped, noFeature and ambiguous counts
counts <- get_data(config::get("counts")$synID)[-c(1:4),]
