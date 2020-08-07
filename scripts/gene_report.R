# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(magrittr)
library(stringr)
library(readr)

load("intermediate_data/chromosome_features.Rdata")
load("intermediate_data/ca_coexp.Rdata")
load("intermediate_data/ca_blastp.Rdata")
load("intermediate_data/ca_biogrid.Rdata")
load("intermediate_data/sac_biogrid.Rdata")
#load("intermediate_data/go_pred.Rdata")


gene_gene_report <- function(data){
	report <- data %>%
		dplyr::select(feature_name_1, feature_name_2)

  # Gene Annotations
	report <- report %>%
		dplyr::left_join(
			chromosome_features %>%
				dplyr::select(
					 feature_name_1 = feature_name,
					 gene_name_1 = gene_name,
					 description_1 = description,
					 sac_ortholog_1 = sac_ortholog),
			by=c("feature_name_1")) %>%
		dplyr::left_join(
			chromosome_features %>%
				dplyr::select(
					 feature_name_2 = feature_name,
					 gene_name_2 = gene_name,
					 description_2 = description,
					 sac_ortholog_2 = sac_ortholog),
			by=c("feature_name_2"))

	# Co-Expression
	report <- report %>%
		dplyr::left_join(
			ca_coexp %>% dplyr::rename(coexp_score = score),
			by=c("feature_name_1", "feature_name_2"))

	# Sequence Similarity
	report <- report %>%
		dplyr::left_join(
			ca_blastp %>%
				dplyr::select(
					feature_name_1,
					feature_name_2,
					blastp_EValue = EValue),
			by=c("feature_name_1", "feature_name_2"))

	# Genetic Interaction
	report <- report %>%
		dplyr::left_join(
			ca_biogrid %>%
				dplyr::semi_join(data, by=c("feature_name_1", "feature_name_2")) %>%
				dplyr::filter(experimental_system_type == "genetic") %>%
				dplyr::group_by(feature_name_1, feature_name_2) %>%
				dplyr::summarize(
					genetic_interaction = paste0(
						experimental_system_abbreviation, "_", throughput_abbreviation, collapse="|")),
			by=c("feature_name_1", "feature_name_2"))

	# Physical Interaction
	report <- report %>%
		dplyr::left_join(
			ca_biogrid %>%
				dplyr::semi_join(data, by=c("feature_name_1", "feature_name_2")) %>%
				dplyr::filter(experimental_system_type == "physical") %>%
				dplyr::group_by(feature_name_1, feature_name_2) %>%
				dplyr::summarize(
					physical_interaction = paste0(
						experimental_system_abbreviation, "_", throughput_abbreviation, collapse="|")),
			by=c("feature_name_1", "feature_name_2"))

	# Sac Ortholog Genetic Interaction
	report <- report %>%
		dplyr::left_join(
			sac_biogrid %>%
				dplyr::semi_join(data, by=c("feature_name_1", "feature_name_2")) %>%
				dplyr::filter(experimental_system_type == "genetic") %>%
				dplyr::group_by(gene_symbol_1, gene_symbol_2) %>%
				dplyr::summarize(
					sac_genetic_interaction = paste0(
						experimental_system_abbreviation, "_", throughput_abbreviation, collapse="|")),
			by=c("sac_ortholog_1" = "gene_symbol_1", "sac_ortholog_2" = "gene_symbol_2"))

	# Sac Ortholog Physical Interaction
	report <- report %>%
		dplyr::left_join(
			sac_biogrid %>%
				dplyr::semi_join(data, by=c("feature_name_1", "feature_name_2")) %>%
				dplyr::filter(experimental_system_type == "physical") %>%
				dplyr::group_by(gene_symbol_1, gene_symbol_2) %>%
				dplyr::arrange(desc(throughput_abbreviation), experimental_system_abbreviation) %>%
				dplyr::summarize(
					sac_physical_interaction = paste0(
						experimental_system_abbreviation, "_", throughput_abbreviation, collapse="|")) %>%
				dplyr::ungroup(),
			by=c("sac_ortholog_1" = "gene_symbol_1", "sac_ortholog_2" = "gene_symbol_2"))

}

