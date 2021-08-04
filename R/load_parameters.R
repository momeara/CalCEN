

# Load parameters
#' @export
load_parameters <- function(
    path = "parameters.yaml") {
    yaml::read_yaml(file = path)
}
