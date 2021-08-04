

#' Get connection to the sequence read archive database
#'
#' if it doesn't exist download it to
#'
#'   <parameters$data_path$scratch_dir>/sra_meta/SRAmetadb.sqlite
#'
#' @export
get_sra_connection <- function(parameters, verbose=TRUE) {
    sra_fname <- paste0(
        parameters$data_paths$scratch_dir, "/sra_meta/SRAmetadb.sqlite")

    if (verbose) {
        cat("Getting SRA connection to database ", sra_fname, "\n", sep = "")
    }
    
    if (!file.exists(sra_fname)) {
        if (verbose) {
            cat("SRA database does not exist, downloading...\n")
        }
        if (!dir.exists(paste0(parameters$data_paths$scratch_dir, "/sra_meta"))) {
            if (verbose) {
                cat("Creating directory ", parameters$data_paths$scratch_dir, "/sra_meta ...\n", sep = "")
            }
            dir.create(paste0(parameters$data_paths$scratch_dir, "/sra_meta"))
        }
	SRAdb::getSRAdbFile(
            destdir = paste0(parameters$data_paths$scratch_dir, "/sra_meta"),
            destfile = "SRAmetadb.sqlite.gz")
    }
    sra_con <- DBI::dbConnect(RSQLite::SQLite(), sra_fname)
}

#'Retrieve SRA runs from the SRAmeta database
#'
#' @param sra_con a DBI connection to the SRA database. E.g. the
#'     result of calling get_sra_connection()
#' @param curated_studies a character vector with study_accession values
#'     or a data.frame with a column study_accession
#' @param verbose (default: TRUE)
#'
#' @return data.frame with one run per row with columns in
#'     curated_studies an additionally:
#' 
#'         study_alias
#'         study_title
#'         experiment_accession
#'         sample_accession
#'         library_strategy
#'         library_layout
#'         sample_alias
#'         taxon_id
#'         scientific_name
#'         common_name
#'         run_accession
#'         run_alias
#' 
#' @export
retrieve_sra_runs <- function (
    sra_con,
    curated_studies,
    verbose = TRUE) {

    if (is.character(curated_studies)) {
        curated_studies <- data.frame(study_accession = curated_studies)
    }
    if (!("study_accession" %in% names(curated_studies))) {
        cat("ERROR: curated_studies must be a data.frame with a column column 'study_accession'\n")
    }

    if (verbose) {
        cat("Copying curated studies to sra database temporary ...\n")
    }
    sra_con %>%
        dplyr::copy_to(
            df = curated_studies,
            overwrite = TRUE)

    if (verbose) {
        cat("Checking if the studies exist in SRA ...\n")
    }
    # Curated studies that can't be found in SRA
    missing_studies <- sra_con %>%
        dplyr::tbl("curated_studies") %>%
        dplyr::anti_join(
            sra_con %>%
            dplyr::tbl("study") %>%
            dplyr::select(study_accession),
            by = "study_accession") %>%
        dplyr::collect(n = Inf)
    if (nrow(missing_studies) > 0) {
        cat("WARNING: ", nrow(missing_studies), " studies could not be found in the SRA database:\n", sep = "")
        print(missing_studies)
    }

    if (verbose) {
        cat("Retrieving runs for ", nrow(curated_studies), " curated studies ...\n", sep = "")
    }
    runs <- sra_con %>%
        dplyr::tbl("curated_studies") %>%
        dplyr::left_join(
            sra_con %>%
            dplyr::tbl("study") %>%
            dplyr::select(
                study_accession,
                study_alias,
                study_title,
                submission_accession),
            by = "study_accession") %>%
        dplyr::left_join(
            sra_con %>%
            dplyr::tbl("experiment") %>%
            dplyr::filter(library_strategy == "RNA-Seq") %>%
            dplyr::select(
                study_accession,
                experiment_accession,
                sample_accession,
                library_strategy,
                library_layout),
            by = "study_accession") %>%
        dplyr::left_join(
            sra_con %>%
            dplyr::tbl("sample") %>%
            dplyr::select(
                sample_accession,
                sample_alias,
                taxon_id,
                scientific_name,
                common_name),
            by = "sample_accession") %>%
        dplyr::left_join(
            sra_con %>%
            dplyr::tbl("run") %>%
            dplyr::select(
                experiment_accession,
                run_accession,
                run_alias),
            by = "experiment_accession") %>%
        dplyr::collect(n = Inf)
}
