# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(EGAD)

load("intermediate_data/ca_biogrid.Rdata")
load("intermediate_data/ca_genes.Rdata")

ca_ppi_network <- ca_biogrid %>%
	dplyr::distinct(feature_name_1, feature_name_2) %>%
	as.data.frame() %>%
	EGAD::build_binary_network(ca_genes) %>%
	EGAD::extend_network()

save(ca_ppi_network, file="intermediate_data/ca_ppi_network.Rdata")
ca_ppi_network %>%
	as.data.frame() %>%
	readr::write_tsv("product/ca_ppi_network_180621.tsv")
