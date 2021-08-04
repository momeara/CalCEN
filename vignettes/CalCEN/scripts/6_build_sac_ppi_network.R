# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(EGAD)

load("intermediate_data/ca_genes.Rdata")
load("intermediate_data/sac_biogrid.Rdata")
load("intermediate_data/sac_sge_interactions.Rdata")

# e.g.
#  feature_name_1 feature_name_2 interactions
#  <chr>          <chr>          <chr>
#1 C2_09900C_A    C2_04680W_A    3*AP-MS_lt|23*AP-W_lt|CoF_lt|2*CoP_lt|3*ReC_lt
ca_sac_ortholog_physical_ppi_summary <- sac_biogrid %>%
	dplyr::filter(experimental_system_type == "physical") %>%
		dplyr::transmute(
			feature_name_1 = ca_feature_name_1,
			feature_name_2 = ca_feature_name_2,
			sac_phys_ppi = paste(
				experimental_system_abbreviation,
				throughput_abbreviation, sep="_")) %>%
	dplyr::filter(!is.na(feature_name_1), !is.na(feature_name_2)) %>%
	dplyr::group_by(feature_name_1, feature_name_2, sac_phys_ppi) %>%
	dplyr::summarize(
			n_ppis = paste0(ifelse(n() == 1, "", paste0(n(), "*")), sac_phys_ppi[1])) %>%
	dplyr::summarize(
			interactions = n_ppis %>% paste0(collapse="|")) %>%
	dplyr::ungroup()
save(ca_sac_ortholog_physical_ppi_summary, file="intermediate_data/ca_sac_ortholog_physical_ppi_summary.Rdata")


ca_sac_ortholog_genetic_ppi_network <- sac_biogrid %>%
	dplyr::filter(experimental_system_type == "genetic") %>%
	dplyr::distinct(ca_feature_name_1, ca_feature_name_2) %>%
	as.data.frame() %>%
	EGAD::build_binary_network(ca_genes) %>%
	EGAD::build_coexp_network(ca_genes)

ca_sac_ortholog_physical_ppi_network <- sac_biogrid %>%
	dplyr::filter(experimental_system_type == "physical") %>%
	dplyr::distinct(ca_feature_name_1, ca_feature_name_2) %>%
	as.data.frame() %>%
	EGAD::build_binary_network(ca_genes) %>%
	EGAD::extend_network()
ca_sac_ortholog_physical_ppi_network[is.na(ca_sac_ortholog_physical_ppi_network)] <- 0

save(
	ca_sac_ortholog_genetic_ppi_network,
	file="intermediate_data/ca_sac_ortholog_genetic_ppi_network.Rdata")
save(
	ca_sac_ortholog_physical_ppi_network,
	file="intermediate_data/ca_sac_ortholog_physical_ppi_network.Rdata")

ca_sac_ortholog_genetic_ppi_network %>%
	as.data.frame() %>%
	readr::write_tsv("product/ca_sac_ortholog_genetic_ppi_network_181010.tsv")
ca_sac_ortholog_physical_ppi_network %>%
	as.data.frame() %>%
	readr::write_tsv("product/ca_sac_ortholog_physical_ppi_network_181010.tsv")


ca_sac_ortholog_genetic_ppi <- ca_sac_ortholog_genetic_ppi_network %>%
	data.frame() %>%
	tibble::rownames_to_column("feature_name_1") %>%
	tidyr::gather(key="feature_name_2", value="score", -feature_name_1)
save(ca_sac_ortholog_genetic_ppi, file="intermediate_data/ca_sac_ortholog_genetic_ppi.Rdata")

ca_sac_ortholog_physical_ppi <- ca_sac_ortholog_physical_ppi_network %>%
	data.frame() %>%
	tibble::rownames_to_column("feature_name_1") %>%
	tidyr::gather(key="feature_name_2", value="score", -feature_name_1)
save(ca_sac_ortholog_physical_ppi, file="intermediate_data/ca_sac_ortholog_physical_ppi.Rdata")


#############################################################################33


ca_sac_sge_ortholog_physical_ppi_summary <- sac_sge_interactions %>%
	dplyr::filter(interaction_type == "Physical") %>%
		dplyr::transmute(
			feature_name_1 = ca_feature_name_1,
			feature_name_2 = ca_feature_name_2,
			sac_phys_ppi = paste(
				experimental_system_abbreviation,
				annotation_abbreviation, sep="_")) %>%
	dplyr::filter(!is.na(feature_name_1), !is.na(feature_name_2)) %>%
	dplyr::group_by(feature_name_1, feature_name_2, sac_phys_ppi) %>%
	dplyr::summarize(
			n_ppis = paste0(ifelse(n() == 1, "", paste0(n(), "*")), sac_phys_ppi[1])) %>%
	dplyr::summarize(
			interactions = n_ppis %>% paste0(collapse="|")) %>%
	dplyr::ungroup()
save(ca_sac_sge_ortholog_physical_ppi_summary, file="intermediate_data/ca_sac_sge_ortholog_physical_ppi_summary.Rdata")

