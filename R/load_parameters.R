

#' Load parameters
#'
#' @param path Location to parameters.yaml file (default: parameters.yaml)
#' @param verbose (default: TRUE)
#' @export
load_parameters <- function(
    path = "parameters.yaml",
    verbose = TRUE) {

    if (verbose) {
        cat("Loading prameters from '", path, "' ...\n", sep = "")
    }

    if (!file.exists(path)) {
        cat("ERROR: parameters path '", path, "' does not exist\n", sep = "")
    }
    yaml::read_yaml(file = path)
}
