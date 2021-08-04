# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(CalCEN)


parameters <- CalCEN::load_parameters()
sra_dir <- paste0(parameters$data_paths$scratch_dir, "/sra")
runs <- readr::read_tsv("intermediate_data/runs_20210803.tsv")

CalCEN::download_expression_data(
    runs = runs,
    sra_dir = sra_dir)

validated_runs <- sra_dir %>%
		CalCEN::validate_downloaded_runs()

cat("Downloaded but in consinstent runs:\n")
validated_runs %>%
		dlpyr::filter(!is_consistent)

cat("Missing runs:\n")
runs %>%
		dplyr::anti_join(
				validated_runs,
				by = "run_accession")

# remove broken runs
validated_runs %>%
    dplyr::filter(!is_consistent) %>%
    plyr::a_ply(1,
				function(ca_run) { system2("rm", ca_run$sra_fname[1]) })
