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









