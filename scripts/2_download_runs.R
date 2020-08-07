# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:


library(SRAdb)
library(plyr)
library(dplyr)
library(stringr)
library(magrittr)
library(glue)
source("scripts/parameters.R")

ca_runs <- readr::read_tsv("product/ca_runs_180608.tsv")
SRAmetadb_fname <- paste0(base_dir, "/sra_meta/SRAmetadb.sqlite")
sra_dir <- paste0(base_dir, "/sra_80601")


if(!file.exists(SRAmetadb_fname)){
	SRAdb::getSRAdbFile(destfile = SRAmetadb_fname)
}

sra_con <- DBI::dbConnect(RSQLite::SQLite(), SRAmetadb_fname)


done <- FALSE
while(!done){
	retrieved_runs <- list.files(
		path=sra_dir,
		pattern="*.sra") %>%
		stringr::str_extract("^[^.]+") %>%
		tibble::data_frame(run_accession = .)

	to_get_runs <- ca_runs %>%
		dplyr::anti_join(
			retrieved_runs, by="run_accession")

	n_runs <- to_get_runs %>% nrow
	if(to_get_runs %>% nrow == 0){
		cat("Done!\n")
		done <<- TRUE
	}

	cat("Downloading {n_runs} more runs ...\n\n" %>% glue::glue())

	# this can take ~ 2 days on a university network
	to_get_runs %>%
		plyr::a_ply(1, function(ca_run){
			cat("Download SRA run '", ca_run$run_accession[1], "' for study '", ca_run$study_accession[1], "' (", n_runs, " runs to go)\n", sep="")

		tryCatch({
			SRAdb::getSRAfile(
				in_acc=ca_run$run_accession[1],
				sra_con=sra_con,
				destDir=sra_dir,
				fileType='sra',
				makeDirectory=TRUE)

			n_runs <<- n_runs - 1

		},
		error=function(e){
			cat("Failed to download run '", ca_run$run_accession[1], "' with error:\n", paste(e, sep="\n"), sep="")
			system2("rm", c("-rf", paste0(sra_dir, ca_run$run_accession[1], ".sra")))
			Sys.sleep(5)
		})
	})
}

# check retrieved runs for download integrety
retrieved_runs <- list.files(
	path=sra_dir,
	pattern="*.sra") %>%
	stringr::str_extract("^[^.]+") %>%
	tibble::data_frame(run_accession = .)

validated_runs <- retrieved_runs %>%
	plyr::adply(1, function(ca_run){
		sra_fname <- paste0(sra_dir, "/", ca_run$run_accession, ".sra")
		cat("checking SRA run '", ca_run$run_accession[1], "': ", sra_fname, "\n", sep="")

		is_consistent <- system2("vdb-validate", sra_fname, stdout=TRUE, stderr=TRUE) %>%
			stringr::str_detect("is consistent") %>%
			any
		tibble::data_frame(
			sra_fname=sra_fname,
			is_consistent = is_consistent)
	})

# remove broken runs
"Found {n_good} good and {n_bad} bad downloaded runs. Deleting the bad ones...\n" %>%
	glue::glue(
		n_good=validated_runs %>% dplyr::filter(is_consistent) %>% nrow,
		n_bad=validated_runs %>% dplyr::filter(!is_consistent) %>% nrow) %>%
	cat

validated_runs %>%
	dplyr::filter(!is_consistent) %>%
	plyr::a_ply(1, function(ca_run){ system2("rm", ca_run$sra_fname[1]) })
