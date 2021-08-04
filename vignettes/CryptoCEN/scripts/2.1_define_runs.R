
# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(CalCEN)
library(SRAdb)
library(plyr)
library(dplyr)
library(stringr)
library(magrittr)
library(readxl)




parameters <- CalCEN::load_parameters()


curated_studies <- readxl::read_xlsx("raw_data/studies_20210726.xlsx", 1) %>%
    dplyr::select(
        -taxon_id,
        -study_title)

# strategy for looking up runs from study_accession
#   1) use the SRAdb
#       a) experiment(study_accession) -> run(experiment_accession)
#       b) study -> run(submission_accession)
#   2) hand curated runs

runs_hand_curated <- readr::read_tsv("raw_data/runs_hand_curated_20210802.tsv")

sra_con <- CalCEN::get_sra_connection(
    parameters = parameters,
    verbose = TRUE)

runs_via_experiment <- sra_con %>%
    CalCEN::retrieve_sra_runs(
        curated_studies = curated_studies)

runs_via_submission <- sra_con %>%
    dplyr::copy_to(
        df = curated_studies,
        overwrite = TRUE) %>%
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
        dplyr::tbl("run") %>%
        dplyr::select(
            submission_accession,
            run_accession,
            run_alias),
        by = "submission_accession") %>%
    dplyr::collect(n = Inf)


runs <- curated_studies %>%
    dplyr::left_join(
        dplyr::bind_rows(
            runs_hand_curated %>%
            dplyr::select(study_accession, run_accession),
            runs_via_experiment %>%
            dplyr::select(study_accession, run_accession),
            runs_via_submission %>%
            dplyr::select(study_accession, run_accession)) %>%
        dplyr::distinct(study_accession, run_accession),
        by = "study_accession")


$ study_accession                        <chr> "SRP253688"
$ submission_accession                   <chr> "SRA1058183"


sra_con %>% dplyr::tbl("study") %>% dplyr::filter(study_accession == "SRP253688")
$ study_accession      <chr> "SRP253688"
$ submission_accession <chr> "SRA1058183"

sra_con %>% dplyr::tbl("run") %>% dplyr::filter(run_accession == "SRR11391643")
$ submission_accession <chr> "SRA1058183"
$ run_accession        <chr> "SRR11391643"


runs <- %>% dplyr::left_join(
    dplyr::


runs_via_submission <- sra_con %>%
    CalCEN::


studies_missing_experimental_details <-
    curated_studies %>%
    dplyr::semi_join(
        runs %>% dplyr::filter(is.na(run_accession)),
        by = "study_accession")

# filter to relevant runs
runs <- runs %>% dplyr::filter(!is.na(run_accession))
runs <- runs %>%
    dplyr::filter(
        !(study_accession == "SRP015456") | (scientific_name %>% stringr::str_detect("H99")),
        !(study_accession == "SRP065179") | (taxon_id == 178876),
        !(study_accession == "SRP073615") | (taxon_id == 178876),
        !(study_accession == "SRP165294") | (run_accession %in% c(
            "SRR8042971", "SRR8042973", "SRR8042977", "SRR8042978",
            "SRR8042984", "SRR8042990", "SRR8042991", "SRR8042972",
            "SRR8042979", "SRR8042983", "SRR8042985", "SRR8042989")),
        !(study_accession == "SRP212774") | (run_accession %in% hand_curated_runs$run_accession))

# do the reported number of runs match the number of runs collected from SRA?
runs %>%
    dplyr::distinct(study_accession, .keep_all = TRUE) %>%
    dplyr::left_join(
        runs %>% dplyr::count(study_accession, sort = T),
        by = "study_accession") %>%
    dplyr::filter(n != n_runs) %>%
    dplyr::select(study_accession, n, n_runs)

# get see if we can get missing studies
sra_con %>%
    dplyr::copy_to(
        df = studies_missing_experimental_details,
        overwrite = TRUE)
runs2 <- sra_con %>%
    dplyr::tbl("studies_missing_experimental_details") %>%
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
        dplyr::tbl("run") %>%
        dplyr::select(
            submission_accession,
            run_accession,
            run_alias),
        by = "submission_accession") %>%
    dplyr::collect(n = Inf)


studies_missing_experimental_details2 <- runs2 %>%
    dplyr::filter(is.na(run_accession))

runs2 <- runs2 %>%
    dplyr::filter(!is.na(run_accession))


runs2 <- runs2 %>%
    dplyr::inner_join(hand_curated_runs, by = c("study_accession", "run_accession"))


# do the reported number of runs match the number of runs collected from SRA?
runs2 %>%
    dplyr::distinct(study_accession, .keep_all = TRUE) %>%
    dplyr::left_join(
        runs2 %>% dplyr::count(study_accession),
    by = "study_accession") %>%
    dplyr::filter(n != n_runs) %>%
    dplyr::select(study_accession, n, n_runs)
    

runs2 %>%
    dplyr::distinct(study_accession, .keep_all = TRUE) %>%
    dplyr::left_join(
        runs2 %>% dplyr::count(study_accession, sort = T),
        by = "study_accession") %>%
    dplyr::filter(n != n_runs)

