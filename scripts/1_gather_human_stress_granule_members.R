# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(stringr)
library(readr)
library(tidyr)
library(magrittr)
library(googledrive)
googledrive::drive_auth("~/.R/googledrive_auth.rds")

library(googlesheets)
googlesheets::gs_auth("~/.R/googledrive_auth.rds")


# retrieved from https://en.wikipedia.org/wiki/Stress_granule 18/07/12
human_stress_granules <- googledrive::as_dribble(
	x="~/Collaborations/Candida Functional Genomics/HSP90 Physical Interactors/Datasets/Stress Granules") %>%
	magrittr::extract2("id") %>%
	googlesheets::gs_key() %>%
	googlesheets::gs_read(
		ws="Human Wikipedia",
		col_types = readr::cols(
			`Gene ID` = col_character(),
			`Protein Name` = col_character(),
			Description = col_character(),
			References = col_character(),
			`Also found in processing bodies?` = col_character())) %>%
	dplyr::rename(
		gene_name = `Gene ID`,
		protein_name = `Protein Name`,
		description = Description,
		references = References,
		in_p_bodies = `Also found in processing bodies?`) %>%
	dplyr::mutate(
		in_p_bodies = !is.na(in_p_bodies))

save(human_stress_granules, file="intermediate_data/human_stress_granules.Rdata")
