#'
#'Detect file type and download data
#'
#'This function takes synIds and version number to download any rectangular file type from Synapse.
#'
#'@param synID A character vector with a Synapse Id.
#'@param version Optional. A numeric vector with the Synapse file version number.
#'@export
#'@return A tibble.
#'@examples
#'\dontrun{
#'file <- get_data(synID = "syn1234", version = 7)
#'
#'}
get_data <- function(synID, version = NULL){
  df <- as_tibble(data.table::fread(synapser::synGet(synID, version = as.numeric(version))$path))
  df
}
