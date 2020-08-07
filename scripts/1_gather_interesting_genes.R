# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(stringr)
library(readr)
library(tidyr)
library(magrittr)
library(googledrive)
googledrive::drive_auth(email='mattjomeara@gmail.com')

library(googlesheets)
# z <- googlesheets::gs_auth(); saverRDS(z, "~/.R/googlesheets_auth.rds")
googlesheets::gs_auth("~/.R/googlesheets_auth.rds")

gdrive_path <- "~/Collaborations/Candida Functional Genomics/HSP90 Physical Interactors/Case Studies/"
raw_path <- "raw_data"

load("intermediate_data/chromosome_features.Rdata")
load("intermediate_data/sac_stress_granule_annotations.Rdata")
load("intermediate_data/sac_pbody_annotations.Rdata")

# retrieved from genes fo interest from GoogleDrive
genes_of_interest <- googledrive::as_dribble(
	x="Gene Sets of Interest") %>%
	magrittr::extract2("id") %>%
	googlesheets::gs_key() %>%
	googlesheets::gs_read(
		ws="Sheet1",
		col_types = readr::cols(
			Include = readr::col_character(),
			Set = readr::col_character(),
			Order = readr::col_integer(),
			feature_name = readr::col_character(),
			Gene = readr::col_character(),
			Notes = readr::col_character(),
			Reference = readr::col_character())) %>%
	dplyr::rename(
		include = Include,
		set = Set,
		set_order = Order,
		gene_name = Gene,
		notes = Notes,
		reference = Reference) %>%
	dplyr::mutate(include = (include == "TRUE"))


# For the genes without feature_name given, look them up in the CGD
genes_of_interest <- rbind(
	genes_of_interest %>%
		dplyr::filter(!is.na(feature_name)),
	genes_of_interest %>%
		dplyr::filter(is.na(feature_name)) %>%
		dplyr::select(-feature_name) %>%
		dplyr::left_join(
			chromosome_features %>%
				dplyr::filter(feature_name %>% stringr::str_detect("A$")) %>%
				dplyr::select(gene_name, feature_name),
			by="gene_name")) %>%
	dplyr::rename(gene_name_orig = gene_name) %>%
	dplyr::left_join(
		chromosome_features %>%
			dplyr::select(
				feature_name,
				gene_name,
				sac_ortholog,
				description),
		by="feature_name") %>%
	dplyr::arrange(set, set_order) %>%
	dplyr::select(
		include,
		set,
		set_order,
		feature_name,
		gene_name_orig,
		gene_name,
		sac_ortholog,
		description,
		notes,
		reference)

# check if any genes weren't located
genes_of_interest %>% dplyr::filter(is.na(feature_name))


genes_of_interest <- genes_of_interest %>%
	rbind(
		sac_stress_granule_annotations %>%
			dplyr::filter(!is.na(feature_name)) %>%
			dplyr::transmute(
				include=FALSE,
				set="Sac Ortholog Stress Granule",
				set_order=1:n(),
				feature_name,
				gene_name_orig=gene_name,
				gene_name,
				sac_ortholog,
				description,
				notes="SGD go_id=GO:0010494; go_term=cytoplasmic stress granule; evidence=IDA; date_collected=18/08/22",
				reference),
		sac_pbody_annotations %>%
			dplyr::filter(!is.na(feature_name)) %>%
			dplyr::transmute(
				include=FALSE,
				set="Sac Ortholog P-body",
				set_order=1:n(),
				feature_name,
				gene_name_orig=gene_name,
				gene_name,
				sac_ortholog,
				description,
				notes="SGD go_id=GO:0000932; go_term=P-body; evidence=IDA; date_collected=18/08/22",
				reference))

save(genes_of_interest, file="intermediate_data/genes_of_interest.Rdata")
