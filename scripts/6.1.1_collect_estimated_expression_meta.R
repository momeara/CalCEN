# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(magrittr)
library(stringr)
library(readr)

ca_runs <- readr::read_tsv("product/ca_runs_180608.tsv")


estimated_expression_meta <- list.files(
	path="intermediate_data/estimated_expression/logs",
	pattern="*log") %>%
	data.frame(log_fname=., stringsAsFactors=FALSE) %>%
	dplyr::mutate(
		run_accession = log_fname %>% stringr::str_extract("^[^.]+")) %>%
	dplyr::left_join(
		ca_runs %>% dplyr::select(run_accession, is_paired),
		by="run_accession") %>%
	plyr::adply(1, function(run){
		log_fname <- run$log_fname[1]
		# cat("processing log file: ", log_fname, " ...\n", sep="")
		lines <- readr::read_lines(
			file=paste0("intermediate_data/estimated_expression/logs/", log_fname))

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


		tibble::data_frame(
			n_reads = n_reads,
			n_reads_aligned_0_times = lines %>%
				stringr::str_detect(" aligned[ a-z]* 0 times") %>%
				magrittr::extract(lines, .) %>%
				stringr::str_match("^ +([0-9]+)") %>%
				unlist() %>%
				magrittr::extract2(2) %>%
				as.numeric(),
			n_reads_aligned_1_time = lines %>%
				stringr::str_detect(" aligned[ a-z]* exactly 1 time") %>%
				magrittr::extract(lines, .) %>%
				stringr::str_match("^ +([0-9]+)") %>%
				unlist() %>%
				magrittr::extract2(2) %>%
				as.numeric(),
			n_reads_aligned_mutiple_times = lines %>%
				stringr::str_detect(" aligned[ a-z]* >1 times") %>%
				magrittr::extract(lines, .) %>%
				stringr::str_match("^ +([0-9]+)") %>%
				unlist() %>%
				magrittr::extract2(2) %>%
				as.numeric(),
			run_time = lines %>%
				stringr::str_detect("# Runtime: ") %>%
				magrittr::extract(lines, .) %>%
				stringr::str_extract("[0-9.]+$") %>%
				as.numeric())
	})
save(estimated_expression_meta, file="intermediate_data/estimated_expression_meta.Rdata")
estimated_expression_meta %>%
	readr::write_tsv("product/estimated_expression_meta_180621.tsv")



