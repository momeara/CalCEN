# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(stringr)
library(readr)
library(BioChemPantry)

ca_runs <- readr::read_tsv("product/ca_runs_180608.tsv")

gather_estimated_expression <- function(path){
	cat("Gathering estimated expression for runs in '", path, "' ...\n", sep="")

	estimated_expression <- list.files(
		path=path,
		pattern="*genes.results") %>%
		plyr::ldply(function(results_fname){
			results <- readr::read_tsv(
				file=paste0(path, "/", results_fname),
				col_types=readr::cols(
					gene_id = readr::col_character(),
					`transcript_id(s)` = readr::col_character(),
					length = readr::col_double(),
					effective_length = readr::col_double(),
					expected_count = readr::col_double(),
					TPM = readr::col_double(),
					FPKM = readr::col_double())) %>%
				dplyr::rename(transcript_ids = `transcript_id(s)`) %>%
				dplyr::mutate(
					run_accession = results_fname %>% stringr::str_extract("^[^.]+"))
		}) %>%
		dplyr::inner_join(
			ca_runs %>% dplyr::select(study_accession, run_accession, is_paired),
			by="run_accession")
}

estimated_expression <- rbind(
	gather_estimated_expression("intermediate_data/estimated_expression"),
	gather_estimated_expression("intermediate_data/estimated_expression_180611"))

# gene_id and transcript_ids are 1-1
#estimated_expression %>% BioChemPantry::summarize_map("gene_id", "transcript_ids")
estiamted_expression <- estimated_expression %>% dplyr::select(-transcript_ids)

save(estimated_expression, file="intermediate_data/estimated_expression.Rdata")
estimated_expression %>%
	readr::write_tsv("product/estimated_expression_180621.tsv")
