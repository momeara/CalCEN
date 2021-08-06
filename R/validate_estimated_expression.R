#' Validate estimated expression runs
#'
#' After running estimate_expression(..., results_dir = <results_dir>, ...)
#' call validate_estimated_runs(results_dir = <results_dir>)
#' to get the runs that have been completed
#' 
#' @param results_dir directory to look for completed estimated expression runs
#' @param check_for_logs check for logs
#' @return data.frame with column [run_accession] of runs that have been completed
#' @export
validate_estimated_expression <- function(
    path,
    check_for_logs = TRUE) {
    if (verbose) {
        cat("Getting estimated expression runs from '", path, "' ...\n", sep = "")
    }
    done_run_results <- list.files(
        path = path,
        pattern = "*.genes.results") %>%
        stringr::str_extract("^[^.]+") %>%
        tibble::tibble(run_accession = .)
    if (check_for_logs) {
        done_run_logs <- list.files(
            path = paste0(path, "/logs"),
            pattern = "*.log") %>%
            stringr::str_extract("^[^.]+") %>%
            tibble::tibble(run_accession = .)
        done_runs <- done_run_results %>%
                dplyr::inner_join(done_run_logs, by="run_accession")
    } else{
        done_runs <- done_run_results
    }
}
