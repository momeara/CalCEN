# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(magrittr)
library(stringr)
library(readr)

ca_runs <- dplyr::bind_rows(
		readr::read_tsv("product/ca_runs_180608.tsv") %>% dplyr::mutate(year_collected = 2018),
		readr::read_tsv("product/ca_runs_200928.tsv") %>% dplyr::mutate(year_collected = 2020))

log_files <- tibble::tibble(
		log_dir = c(
				"intermediate_data/estimated_expression/logs",
				"intermediate_data/estimated_expression_180611/logs",
				"intermediate_data/estimated_expression_20201007/logs")) %>%
		dplyr::rowwise() %>%
		dplyr::do({
				list.files(path = .$log_dir[1], pattern="*log", full.names=TRUE) %>%
						tibble::tibble(log_fname= .)
		}) %>%
		dplyr::mutate(
				run_accession = log_fname %>%
						stringr::str_match("([^/]+)[.]log") %>%
		        magrittr::extract(,2)) %>%
		dplyr::inner_join(
				ca_runs %>% dplyr::select(run_accession, is_paired),
				by="run_accession")

estimated_expression_meta <- log_files  %>%
	plyr::adply(1, function(run){
		log_fname <- run$log_fname[1]
		cat("processing log file: ", log_fname, " ...\n", sep="")
		lines <- readr::read_lines(file=log_fname)
		n_reads <- lines %>%
			stringr::str_detect(" reads; of these:") %>%
			magrittr::extract(lines, .) %>%
			stringr::str_extract("^[0-9]+") %>%
			as.numeric()
		if(!run$is_paired[1]){
			n_reads_unpaired <- lines %>%
				stringr::str_detect(" were unpaired; of these:") %>%
				magrittr::extract(lines, .) %>%
				stringr::str_match("^  ([0-9]+)") %>%
				unlist() %>%
				magrittr::extract2(2) %>%
				as.numeric()
			if(n_reads_unpaired != n_reads){
				cat("The run is unpaired and there were ", n_reads, " reads but only ", n_unpaired, " unpaired reads.\n", sep="")
			}
		} else {
			n_reads_paired <- lines %>%
				stringr::str_detect(" were paired; of these:") %>%
				magrittr::extract(lines, .) %>%
				stringr::str_match("^  ([0-9]+)") %>%
				unlist() %>%
				magrittr::extract2(2) %>%
				as.numeric()
			if(n_reads_paired != n_reads){
				cat("The run is paired and there were ", n_reads, " reads but only ", n_unpaired, " paired reads.\n", sep="")
			}
		}
		n_reads_aligned_0_times = lines %>%
				stringr::str_detect(" aligned[ a-z]* 0 times") %>%
				magrittr::extract(lines, .) %>%
				stringr::str_match("^ +([0-9]+)") %>%
				unlist() %>%
				magrittr::extract2(2) %>%
				as.numeric()
		n_reads_aligned_1_time <- lines %>%
				stringr::str_detect(" aligned[ a-z]* exactly 1 time") %>%
				magrittr::extract(lines, .) %>%
				stringr::str_match("^ +([0-9]+)") %>%
				unlist() %>%
				magrittr::extract2(2) %>%
				as.numeric()
		n_reads_aligned_mutiple_times <- lines %>%
				stringr::str_detect(" aligned[ a-z]* >1 times") %>%
				magrittr::extract(lines, .) %>%
				stringr::str_match("^ +([0-9]+)") %>%
				unlist() %>%
				magrittr::extract2(2) %>%
				as.numeric()
		run_time <- lines %>%
				stringr::str_detect("# Runtime: ") %>%
				magrittr::extract(lines, .) %>%
				stringr::str_extract("[0-9.]+$") %>%
				as.numeric()
		run_time <- ifelse(length(run_time) == 0, NA, run_time)
		tibble::tibble(
				n_reads = n_reads,
				n_reads_aligned_0_times = n_reads_aligned_0_times,
				n_reads_aligned_1_time = n_reads_aligned_1_time,
				n_reads_aligned_mutiple_times =	n_reads_aligned_mutiple_times,
				run_time = run_time)
	})


estimated_expression_meta <- estimated_expression_meta %>%
		dplyr::left_join(
				ca_runs %>%
				dplyr::distinct(
						run_accession,
						study_accession),
				by=c("run_accession"))

save(estimated_expression_meta, file="intermediate_data/estimated_expression_meta.Rdata")

estimated_expression_meta %>%
	readr::write_tsv("product/estimated_expression_meta_20201007.tsv")



