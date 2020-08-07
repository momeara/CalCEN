# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(readr)
library(reshape2)
library(EGAD)

source("scripts/embed_network.R")
load("intermediate_data/ca_silac_hsp90_intensities.Rdata")
load("intermediate_data/ca_genes.Rdata")

ca_silac_network <- ca_silac_hsp90_intensities %>%
	dplyr::filter(
		condition=="wildtype",
		!is.na(feature_name)) %>%
	dplyr::select(feature_name, fraction_index, intensity) %>%
	reshape2::acast(feature_name ~ fraction_index, value.var="intensity") %>%
	EGAD::build_coexp_network(ca_genes)
ca_silac_network[is.na(ca_silac_network)] <- 0

ca_silac_network <- ca_silac_network %>%
	embed_network(ca_genes)

save(
	ca_silac_network,
	file="intermediate_data/ca_silac_network.Rdata")

ca_silac <- ca_silac_network %>%
	data.frame() %>%
	tibble::rownames_to_column("feature_name_1") %>%
	tidyr::gather(key="feature_name_2", value="score", -feature_name_1)
save(
	ca_silac,
	file="intermediate_data/ca_silac.Rdata")
