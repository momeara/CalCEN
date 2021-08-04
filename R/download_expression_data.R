
#' Download RNA-Seq expression runs from SRA
#'
#'   Use prefetch to download data from the SRA to
#'
#'     <sra_dir>/<run_accesion>/<run_accession>.sra
#'
#'   Each expression run is ~8 Gb, so roughly 1000 runs will require
#'   approximately 8 terabytes of space
#' 
#' @param runs: data.frame with columns ['study_accession', 'run_accession']
#' @param sra_dir: directory where to save the expression data
#' 
#' @export
download_expression_data <- function(
    runs,
    sra_dir,
    verbose = TRUE
    ) {

    if (verbose) {
        cat("Checking if output director '", sra_dir, "' exists ...\n", sep = "")
    }
    if (!dir.exists(paths = sra_dir)) {
        cat("Creating ", sra_dir, " ...\n", sep = "")
        dir.create(sra_dir, recursive = TRUE)
    }

    get_runs_remaining <- function(runs, sra_dir) {
        retrieved_runs <- list.files(
            path = sra_dir,
            pattern = "*.sra") %>%
            stringr::str_extract("^[^.]+") %>%
            tibble::tibble(run_accession = .)
    
        runs_remaining <- runs %>%
            dplyr::anti_join(
                retrieved_runs, by = "run_accession")
        cat("Downloading {nrow(runs_remaining)} more runs ...\n\n" %>% glue::glue())
        runs_remaining
    }        

    download_run <- function(run, sra_dir) {
        cat("Download SRA run '", run$run_accession[1], "' for study '", run$study_accession[1], "'\n", sep = "")
        
        tryCatch({
            command <- paste0("cd ", sra_dir, " && ",
                "prefetch ", run$run_accession[1], " ",
                "--progress 1")
            cat("Command: ", command, "\n", sep = "")
            system(command)
   
        },
        error = function(e) {
            cat(
                "Failed to download run '", run$run_accession[1], "' with error:\n",
                paste(e, sep = "\n"), sep = "")
            system2("rm", c("-rf", paste0(sra_dir, run$run_accession[1], ".sra")))
            Sys.sleep(5)
        })
    }
    

    done <- FALSE
    while (!done) {
        runs_remaining <- get_runs_remaining(runs, sra_dir)
        n_runs <- runs_remaining %>% nrow
        if (runs_remaining %>% nrow == 0) {
            cat("Done!\n")
            done <<- TRUE
        }
    
        # this can take ~ 2 days on a university network
        runs_remaining %>%
            plyr::a_ply(1, function(run) {
                status <- download_run(run, sra_dir)
                if (status == "success") {
                    n_runs <<- n_runs - 1
                    cat(paste0("SUCCESS: ", n_runs, " runs to go ...\n", sep = ""))
                    
                }
            })
    }
}

#' Check that the downloaded runs are not corrupted
#'
#'    use vdb-validate from the prefetch package to validate all
#'    all .sra files under <sra_dir>.
#'
#'        <sra_dir>/<run_accession>/<run_accession>.sra
#'
#' @param sra_dir: path to find downloaded SRA runs
#' @export
validate_downloaded_runs <- function(sra_dir) {
    # check retrieved runs for download integrety
    retrieved_runs <- list.files(
        path = sra_dir,
        pattern = "*.sra$",
        recursive = TRUE) %>%
        tibble::tibble(run_fname = .) %>%
        dplyr::mutate(
            run_accession = run_fname %>% stringr::str_extract("^[^/.]+"))

    cat("Found ", nrow(retrieved_runs), " runs to validate ...\n", sep = "")
    
    validated_runs <- retrieved_runs %>%
        plyr::adply(1, function(run) {
            sra_fname <- paste0(sra_dir, "/", run$run_fname[1])
            cat("checking SRA run '", run$run_accession[1], "': ", sra_fname, "\n", sep="")
            is_consistent <- system2("vdb-validate", sra_fname, stdout = TRUE, stderr = TRUE) %>%
                stringr::str_detect("is consistent") %>%
                any()
            tibble::tibble(
                sra_fname = sra_fname,
                is_consistent = is_consistent)
        })
    
    # remove broken runs
    "Found {n_good} good and {n_bad} bad runs...\n\n" %>%
        glue::glue(
            n_good = validated_runs %>% dplyr::filter(is_consistent) %>% nrow,
            n_bad = validated_runs %>% dplyr::filter(!is_consistent) %>% nrow) %>%
        cat()
        
    validated_runs
}
