# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(dplyr)
library(EGAD)

load("intermediate_data/ca_blastp.Rdata")
load("intermediate_data/ca_genes.Rdata")

ca_blastp_network <- ca_blastp %>%
	dplyr::select(feature_name_1, feature_name_2, bit_score) %>%
	as.data.frame() %>%
	EGAD::build_weighted_network(ca_genes)

save(ca_blastp_network, file="intermediate_data/ca_blastp_network.Rdata")

ca_blastp_network %>%
	as.data.frame() %>%
	readr::write_tsv("product/ca_blastp_network_20201113.Rdata")

