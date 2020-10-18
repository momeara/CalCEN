# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(magrittr)
library(stringr)
library(readr)

load("intermediate_data/ca_coexp.Rdata")

source("scripts/gene_report.R")

top_ca_coexp <- ca_coexp %>%
	dplyr::filter(feature_name_1 < feature_name_2) %>%
	dplyr::arrange(desc(score)) %>%
	head(50) %>%
	gene_gene_report()


PGA52 <- ca_coexp %>%
	dplyr::filter(feature_name_1 == "C2_00100C_A") %>%
	dplyr::filter(feature_name_2 != feature_name_1) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::select(-score) %>%
	gene_gene_report()

C7_00310C_A <- ca_coexp %>%
	dplyr::filter(feature_name_1 == "C7_00310C_A") %>%
	dplyr::filter(feature_name_2 != feature_name_1) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::select(-score) %>%
	gene_gene_report() %>%
	dplyr::select(
		-description_1,
		-gene_name_1, -sac_ortholog_1, -blastp_EValue, -physical_interaction, -genetic_interaction, -sac_genetic_interaction, -sac_physical_interaction
	) %>%
	dplyr::mutate(description_2 = description_2 %>% stringr::str_sub(1,50)) %>%
	head(50)



FOX2 <- ca_coexp %>%
	dplyr::filter(feature_name_1 == "C3_00810C_A") %>%
	dplyr::filter(feature_name_2 != feature_name_1) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::select(-score) %>%
	gene_gene_report() %>%
	dplyr::select(
		-description_1, -gene_name_1, -sac_ortholog_1,
#		-blastp_EValue, -physical_interaction, -genetic_interaction, -sac_genetic_interaction, -sac_physical_interaction
	) %>%
	dplyr::mutate(description_2 = description_2 %>% stringr::str_sub(1,80)) %>%
	head(50)


#
seeds <- readr::read_tsv("raw_data/deORFinizing_seeds_20201006.tsv")

# all present
missing_seeds <- seeds %>% dplyr::anti_join(
		chromosome_features,
		by = c("Seed" = "feature_name"))



# and in the top 30 associations to a seed
# and itself not a seed
seed_to_neighbor <- ca_coexp %>%
		dplyr::filter(feature_name_1 != feature_name_2) %>%
		dplyr::semi_join(
				seeds,
				by = c("feature_name_1" = "Seed")) %>%
		dplyr::group_by(feature_name_1) %>%
		dplyr::arrange(desc(score)) %>%
		dplyr::slice(1:50) %>%
		dplyr::ungroup() %>%
		dplyr::transmute(
				seed = feature_name_1,
				source_type = "seed",
				target_type = "neighbor",
				feature_name_1,
				feature_name_2)

neighbor_to_neighbor <- seed_to_neighbor %>%
		dplyr::group_by(feature_name_1) %>%
		dplyr::do({
				data <- .
				tidyr::expand_grid(
						feature_name_1 = data$feature_name_2,
						feature_name_2 = data$feature_name_2) %>%
						dplyr::mutate(
								seed = data$feature_name_1[1])}) %>%
		dplyr::ungroup() %>%
		dplyr::filter(feature_name_1 < feature_name_2) %>%
		dplyr::transmute(
				seed,
				source_type = "neighbor",
				target_type = "neighbor",
				feature_name_1,
				feature_name_2)


interactions <- dplyr::bind_rows(
		seed_to_neighbor,
		neighbor_to_neighbor)

interactions <- interactions %>%
		dplyr::left_join(
				interactions %>%
						dplyr::distinct(feature_name_1, feature_name_2) %>%
						gene_gene_report() %>%
						dplyr::distinct(feature_name_1, feature_name_2, .keep_all=T),
				by=c("feature_name_1", "feature_name_2"))

interactions <- interactions %>%
		dplyr::filter(!(
				coexp_score < .9 &
				is.na(blastp_EValue) &
				is.na(genetic_interaction) &
				is.na(physical_interaction) &
				is.na(sac_genetic_interaction) &
				is.na(sac_physical_interaction)))

# the sac intractions are all na to remove them
interactions <- interactions %>%
		dplyr::select(
				-sac_genetic_interaction,
				-sac_physical_interaction)

interactions <- interactions %>%
		dplyr::mutate(
				source_label = ifelse(!is.na(gene_name_1), gene_name_1, feature_name_1),
				target_label = ifelse(!is.na(gene_name_2), gene_name_2, feature_name_2))

interactions %>%
		readr::write_tsv("product/deORFinizing_seed_to_neighborhood_20201007.tsv")


### seeds for kyla
ca_coexp_for_kyla <- ca_coexp %>%
	dplyr::semi_join(
		readr::read_tsv("raw_data/seeds_selvig_20201013.tsv"),
		by = c("feature_name_1" = "feature_name")) %>%
	dplyr::filter(feature_name_1 != feature_name_2) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::group_by(feature_name_1) %>%
	dplyr::slice(1:50) %>%
	dplyr::ungroup() %>%
	gene_gene_report()

ca_coexp_for_kyla %>%
		readr::write_tsv("product/seeds_for_selvid_coexp_20201013.tsv")

##########
# ERG11

report <- ca_coexp %>%
	dplyr::filter(feature_name_1 == "C5_00660C_A") %>%
	dplyr::filter(feature_name_2 != feature_name_1) %>%
	dplyr::arrange(desc(score)) %>%
	dplyr::slice(1:50) %>%
	dplyr::select(-score) %>%
	gene_gene_report()


report %>%
	readr::write_tsv("product/ERG11_gene_gene_report_20201015.tsv")
